"""Persistent AlphaFold 3 process: loads and compiles the model once per job.

Runs with the AF3 module's Python (started by ``trial.run_af3``), so it
imports nothing from ``hpo``. Each line on stdin is one request, a JSON object
with ``input_dir``, ``output_dir`` and ``log``. The worker predicts every input
of ``input_dir`` like ``run_alphafold.py --input_dir --norun_data_pipeline``,
or, with ``--msa_database``, with AF3's data pipeline searching that FASTA.
It writes AF3's own output (including ``Done processing fold input <name>.``) to
``log`` and answers ``done`` on its original stdout. It exits at the end of
stdin and on any error, so a failed prediction never leaves a broken model
state behind; the traceback is in the request's log.
"""
import argparse
import concurrent.futures
import datetime
import json
import os
import pathlib
import shutil
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--model_dir', required=True)
parser.add_argument('--buckets', help='comma-separated; default: those of run_alphafold.py')
parser.add_argument('--msa_database', help='protein FASTA; AF3 searches it with jackhmmer instead of using the input MSA')
args = parser.parse_args()

# Answers keep the original stdout; everything else AF3 prints goes to stderr
# (the log of the trial that started the worker) and later to each request's log.
reply = os.fdopen(os.dup(1), 'w', buffering=1)
os.dup2(2, 1)

sys.path.insert(0, os.environ['ALPHAFOLD_BIN_DIR'])
import jax  # noqa: E402
import run_alphafold  # noqa: E402
from alphafold3.common import folding_input  # noqa: E402
from alphafold3.constants import chemical_components  # noqa: E402
from alphafold3.data import featurisation  # noqa: E402
from alphafold3.data import pipeline  # noqa: E402

# The same model and defaults as run_alphafold.py without further flags.
runner = run_alphafold.ModelRunner(
    config=run_alphafold.make_model_config(),
    device=jax.local_devices(backend='gpu')[0],
    model_dir=pathlib.Path(args.model_dir),
)
buckets = tuple(int(bucket) for bucket in (args.buckets.split(',') if args.buckets
                                           else run_alphafold._BUCKETS.default))

# Featurisation runs on the CPU (1.4 s for 36 nt), inference on the GPU (3.1 s).
# The next input is featurised in a thread while the GPU predicts the current
# one, and process_fold_input receives it through this wrapper: 27% less AF3
# time, identical models (A100, 2026-11-10).
featurise = featurisation.featurise_input
prefetched = {}


def featurise_input(*, fold_input, **kwargs):
    future = prefetched.pop(fold_input.name, None)
    return future.result() if future else featurise(fold_input=fold_input, **kwargs)


# With a data pipeline, AF3 searches the MSA before featurising, so the raw
# input cannot be featurised ahead.
if not args.msa_database:
    featurisation.featurise_input = featurise_input
prefetcher = concurrent.futures.ThreadPoolExecutor(1)

data_pipeline_config = None
if args.msa_database:
    # Only UniRef90 is searched (with the given FASTA). The other protein
    # databases get one unrelated sequence that never passes the E-value, so
    # their MSAs hold only the query. RNA and template databases stay unused:
    # the inputs are protein chains with `templates: []`.
    placeholder = pathlib.Path(args.msa_database).with_name('placeholder.fasta')
    placeholder.write_text('>placeholder\nACDEFGHIKLMNPQRSTVWY\n')
    jackhmmer = shutil.which('jackhmmer')
    if jackhmmer is None:
        raise RuntimeError('jackhmmer not on PATH; the AF3 module should provide it')
    data_pipeline_config = pipeline.DataPipelineConfig(
        jackhmmer_binary_path=jackhmmer,
        nhmmer_binary_path='', hmmalign_binary_path='',
        hmmsearch_binary_path='', hmmbuild_binary_path='',
        uniref90_database_path=args.msa_database,
        small_bfd_database_path=str(placeholder),
        mgnify_database_path=str(placeholder),
        uniprot_cluster_annot_database_path=str(placeholder),
        ntrna_database_path='', rfam_database_path='', rna_central_database_path='',
        seqres_database_path='', pdb_database_path='',
        # The four protein searches run in parallel.
        jackhmmer_n_cpu=max(1, int(os.environ.get('SLURM_CPUS_PER_TASK', 8)) // 4),
        max_template_date=datetime.date(2021, 9, 30),  # run_alphafold.py's default; unused
    )

for line in sys.stdin:
    request = json.loads(line)
    with open(request['log'], 'a') as log:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
    fold_inputs = list(folding_input.load_fold_inputs_from_dir(pathlib.Path(request['input_dir'])))
    for fold_input, upcoming in zip(fold_inputs, fold_inputs[1:] + [None]):
        if upcoming is not None and data_pipeline_config is None:
            # The arguments of run_alphafold.predict_structure.
            prefetched[upcoming.name] = prefetcher.submit(
                featurise, fold_input=upcoming, buckets=buckets,
                ccd=chemical_components.cached_ccd(user_ccd=upcoming.user_ccd),
                verbose=True, conformer_max_iterations=None,
            )
        run_alphafold.process_fold_input(
            fold_input=fold_input,
            data_pipeline_config=data_pipeline_config,
            model_runner=runner,
            output_dir=os.path.join(request['output_dir'], fold_input.sanitised_name()),
            buckets=buckets,
        )
    sys.stdout.flush()
    sys.stderr.flush()
    reply.write('done\n')
