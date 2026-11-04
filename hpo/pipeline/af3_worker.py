"""Persistent AlphaFold 3 process: loads and compiles the model once per job.

Runs with the AF3 module's Python (started by ``trial.run_af3``), so it
imports nothing from ``hpo``. Each line on stdin is one request, a JSON object
with ``input_dir``, ``output_dir`` and ``log``. The worker predicts every input
of ``input_dir`` like ``run_alphafold.py --input_dir --norun_data_pipeline``,
writes AF3's own output (including ``Done processing fold input <name>.``) to
``log`` and answers ``done`` on its original stdout. It exits at the end of
stdin and on any error, so a failed prediction never leaves a broken model
state behind; the traceback is in the request's log.
"""
import argparse
import json
import os
import pathlib
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--model_dir', required=True)
parser.add_argument('--buckets', help='comma-separated; default: those of run_alphafold.py')
args = parser.parse_args()

# Answers keep the original stdout; everything else AF3 prints goes to stderr
# (the log of the trial that started the worker) and later to each request's log.
reply = os.fdopen(os.dup(1), 'w', buffering=1)
os.dup2(2, 1)

sys.path.insert(0, os.environ['ALPHAFOLD_BIN_DIR'])
import jax  # noqa: E402
import run_alphafold  # noqa: E402
from alphafold3.common import folding_input  # noqa: E402

# The same model and defaults as run_alphafold.py without further flags.
runner = run_alphafold.ModelRunner(
    config=run_alphafold.make_model_config(),
    device=jax.local_devices(backend='gpu')[0],
    model_dir=pathlib.Path(args.model_dir),
)
buckets = tuple(int(bucket) for bucket in (args.buckets.split(',') if args.buckets
                                           else run_alphafold._BUCKETS.default))

for line in sys.stdin:
    request = json.loads(line)
    with open(request['log'], 'a') as log:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
    for fold_input in folding_input.load_fold_inputs_from_dir(pathlib.Path(request['input_dir'])):
        run_alphafold.process_fold_input(
            fold_input=fold_input,
            data_pipeline_config=None,
            model_runner=runner,
            output_dir=os.path.join(request['output_dir'], fold_input.sanitised_name()),
            buckets=buckets,
        )
    sys.stdout.flush()
    sys.stderr.flush()
    reply.write('done\n')
