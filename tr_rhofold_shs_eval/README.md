# trRosettaRNA / RhoFold+ SHS evaluation

Feeds this repository's synthetic homologous sequences (SHS) into two *alignment-based*
RNA 3D predictors other than AlphaFold 3 — [trRosettaRNA v1.1](https://yanglab.qd.sdu.edu.cn/trRosettaRNA/)
and [RhoFold+](https://github.com/ml4bio/RhoFold) — and scores the resulting structures
with the evaluation stack of this repo (`../run_evaluation.sh`). It answers the question
"do synthetic homologs help predictors that are not AlphaFold 3?", and supports the
SHS-vs-natural-database-MSA ablation for both predictors.

```
 SHS A3M ─┬─► trRosettaRNA (predict.py + fold.py) ─► *_model1.pdb ─┐
          │                                                        ├─► CIF ─► ../run_evaluation.sh ─► results/*.csv ─► comparison vs AF3
 DB  A3M ─┴─► RhoFold+ (inference.py) ────────────► relaxed_*.pdb ─┘
```

All commands below are run **from this folder** (`tr_rhofold_shs_eval/`) unless stated
otherwise.

---

## 1. What you need

### 1.1 Conda environments

| env | used for | create with |
| --- | --- | --- |
| `synEvo` | SHS generation (imports RnaBench) and the evaluator | `./install.sh` in the repo root, see [Install Guidelines](../docs/INSTALL.md) |
| `rhofold` | RhoFold+ inference | `conda env create -f envs/rhofold_env.yml` |
| `trRNA` | trRosettaRNA inference | `conda env create -f envs/trRosettaRNA_env.yml` |

The two env files are exact exports from the machine the paper runs were produced on;
if `conda` refuses them on your platform, install RhoFold+/trRosettaRNA per their own
instructions and point the scripts at the resulting env with `CONDA_ENV=<name>`.

The SLURM scripts activate `rhofold` or `trRNA` and call the `synEvo` interpreter
directly (`SYNEVO_PY`, default `~/.conda/envs/synEvo/bin/python`) wherever SHS have to be
generated, so no env switch happens inside a job.

Additionally, `convert_prediction_pdbs_to_eval_cif.py` needs **Biopython**
(`pip install biopython`) — it is already in `envs/rhofold_env.yml`.

### 1.2 External predictors (not part of this repo)

Both are expected as **siblings of the SyntheticEvolution checkout**, i.e. `../../RhoFold`
and `../../trRosettaRNA_v1.1` relative to this folder:

```
<parent>/
├── SyntheticEvolution/      <- this repo (this folder lives in here)
├── RhoFold/                 <- inference.py + pretrained checkpoint
└── trRosettaRNA_v1.1/       <- predict.py, fold.py, params/model_1/
```

Every script takes an override, so any other location works:

- RhoFold+: `--rhofold-bin /path/to/RhoFold/inference.py`, `--checkpoint /path/to/model.pt`
  (job scripts: `RHOFOLD_BIN`, `RHOFOLD_CKPT`, `RHOFOLD_ROOT`).
  Only `--msa-mode auto_db` additionally needs RhoFold's `database/` (rnacentral.fasta, nt)
  and `rhofold/data/bin/` (blastn + perl helpers) — `DATABASE_DPATH`, `BINARY_DPATH`.
- trRosettaRNA: `--trrosetta-root /path/to/trRosettaRNA_v1.1` (job scripts: `TRROSETTA_ROOT`).
  The root must contain `predict.py`, `fold.py` and `params/model_1/`. Without an explicit
  `-ss` input trRosettaRNA runs its bundled SPOT-RNA, which is slow; see §4 for the
  base-pair-CSV route that skips it.

### 1.3 Data from this repo

Run `./download.sh` in the repo root (see the [main README](../README.md)) to get:

- `../data/datafiles/alphafold3/` — the AF3 input JSONs whose `unpairedMsa` holds the SHS.
  Needed only if you want to re-derive FASTA/A3M inputs yourself; the prepared folders
  shipped here (§2) already contain them.
- `../evaluation/predictions/gt/` — ground-truth mmCIF structures, required for evaluation.
- `../results/csvs/All_alphafold.csv` — AlphaFold 3 baseline, required for
  `compare_algorithms.py`.

### 1.4 Evaluation tooling

`../run_evaluation.sh` shells out to **US-align** as `./USalign/USalign` *relative to the
repository root*, while the repo ships the source under `evaluation/USalign/`. Before the
first evaluation, make sure the binary is reachable there, e.g. from the repo root:

```bash
cd evaluation/USalign && make && cd ../..
ln -s evaluation/USalign USalign      # so ./USalign/USalign resolves
```

DSSR is **not** needed: `run_evaluation.sh` invokes the evaluator with `-dssr`
(see [docs/EVALUATION.md](../docs/EVALUATION.md)).

### 1.5 Site-specific settings you will want to change

The SLURM scripts are the ones that produced the paper runs, on a
[bwForCluster Helix](https://wiki.bwhpc.de/e/Helix)-style setup. Nothing below is
required for the Python drivers — those run anywhere — but check these before your first
`sbatch`:

| assumption | where | how to change |
| --- | --- | --- |
| `--partition=gpu-single` / `cpu-single` | `#SBATCH` lines in every job script | `sbatch --partition=<yours> ...` overrides the directive, or edit the file |
| `--gres=gpu:A100:1` (some jobs `gpu:A40:1`) | `#SBATCH` lines | `sbatch --gres=gpu:<type>:1 ...`, or drop the type if your site does not use one |
| `module purge && module load devel/cuda` | top of every GPU job script | replace with your site's CUDA module, or delete both lines if CUDA comes from the conda env |
| conda env names `rhofold`, `trRNA` | `CONDA_ENV="${CONDA_ENV:-...}"` | `sbatch --export=ALL,CONDA_ENV=<name> ...` |
| conda base fallback `$HOME/miniconda3` | `CONDA_BASE=$(conda info --base ...)` | only used if `conda info --base` fails; install path is otherwise auto-detected |
| `synEvo` interpreter at `~/.conda/envs/synEvo/bin/python` | `SYNEVO_PY` in `run_rhofold_shs_seeds_job.sh`, `DEFAULT_GEN_PYTHON` in `generate_shs_seed_inputs.py` | `SYNEVO_PY=/path/to/python` / `--gen-python /path/to/python` |
| walltimes, `--mem`, `--cpus-per-task` | `#SBATCH` lines | tuned for ~90 RNA monomers; scale to your queue limits |
| RhoFold+ / trRosettaRNA at `../../` | script defaults | see §1.2 |

No account or project directives (`--account`, `--qos`) are set anywhere, so add yours if
your scheduler requires them. Three scripts carry a disabled mail block
(`##SBATCH --mail-user=you@example.org`) — put your address in and drop one `#` to enable
it.

---

## 2. What is in this folder

### Input preparation

| file | purpose |
| --- | --- |
| `prepare_shs_inputs_from_af3_json.py` | AF3 JSON (`../data/datafiles/alphafold3`) → `fasta/` + `a3m/` + `metadata/`. Single-RNA targets only by default. |
| `generate_shs_seed_inputs.py` | Regenerates one SHS variant under a chosen **generator seed** from `RNA_Monomers_base_pairs.csv`, using `../rna_msa_generator_base_pair.py`. |
| `prepare_standard_msa_from_blastn.py` | Natural database MSAs via RhoFold's BLASTN helper (the non-SHS arm of the ablation). |
| `prepare_msa_plus_shs_inputs.py` | Stacks a natural-DB A3M and an SHS A3M into one combined alignment. |
| `prepare_trrosetta_ss_from_base_pairs_csv.py` | Base-pair lists → per-target `spot_prob` matrices for trRosetta's `-ss`. |
| `RNA_Monomers_base_pairs.csv` | Per-target base pairs for every structure-prediction method (`pairs`, `wc_pairs`, `wobble_pairs`, `nc_pairs`, `lone_pairs`, `multiplets`). |

### Folding

| file | purpose |
| --- | --- |
| `run_rhofold_shs_batch.py` | Batch RhoFold+ over a `fasta/`+`a3m/` folder; `--msa-mode provided \| auto_db \| single_seq`. |
| `run_trrosetta_shs_batch.py` | Batch trRosettaRNA (geometry `predict.py` → `fold.py`); `--msa-mode provided \| single_seq`, optional custom `--ss-dir`. |
| `run_rhofold.py`, `run_trRosettaRNA.py` | Standalone single-folder runners, independent of the SHS layout (useful for smoke tests). |

### Evaluation and comparison

| file | purpose |
| --- | --- |
| `convert_prediction_pdbs_to_eval_cif.py` | Predicted PDBs → mmCIF named `<PDBID>_<algorithm>.cif`, with the `_entity_poly` rows the evaluator's parser needs. |
| `evaluate_with_synthevolution.py` | Runs `../run_evaluation.sh` on those CIFs and files the resulting CSVs under a per-arm label. |
| `compare_algorithms.py` | Merges trRosetta / RhoFold / AlphaFold 3 CSVs into per-target and mean tables (+ bar plot). |
| `compare_msa_sources.py` | SHS vs standard-MSA comparison within one algorithm. |

### SLURM job scripts

Every `*_job.sh` is a thin wrapper: it activates a conda env, sets defaults, and calls one
of the Python drivers. All paths are `VAR="${VAR:-default}"`, so anything can be replaced
via `sbatch --export=ALL,VAR=value`.

| script | env / partition | what it runs |
| --- | --- | --- |
| `run_rhofold_shs_from_{rnafold,spotrna,rnaformerN100}_job.sh` | rhofold / gpu | RhoFold+ on the matching `prepared_shs_inputs_from_*_folder` |
| `run_rhofold_shs_from_rna_monomers_*_job.sh` | rhofold / gpu | RhoFold+ on the four DSSR-denoised (second-round) arms |
| `run_rhofold_msa_plus_shs_from_rnaformerN100_job.sh` | rhofold / gpu | RhoFold+ on natural-MSA + SHS combined inputs |
| `run_rhofold_shs_single_seq_job.sh` | rhofold / gpu | RhoFold+ single-sequence control |
| `run_rhofold_stdmsa_job.sh` | rhofold / gpu | RhoFold+ with its own BLASTN database search (`auto_db`) |
| `run_rhofold_shs_from_spotrna_9CXF_retry_job.sh` | rhofold / gpu | single-target rerun (9CXF has no SPOT-RNA row) |
| `run_rhofold_shs_seeds_job.sh`, `submit_rhofold_shs_seed_jobs.sh` | rhofold + synEvo / gpu array | multi-seed SHS sweep, see §5 |
| `run_trrosetta_shs_job.sh` | trRNA / gpu | trRosettaRNA on `prepared_shs_inputs` |
| `run_trrosetta_shs_job_CPU_{rnafold,spotrna,rnaformerN100}.sh` | trRNA / cpu | trRosettaRNA per SHS arm, CPU-only |
| `run_trrosetta_shs_job_CPU_csv_ss.sh` | trRNA / cpu | converts CSV base pairs to `spot_prob` SS, then folds with them (§4) |
| `run_trrosetta_stdmsa_job.sh` | trRNA / gpu | trRosettaRNA on standard-DB MSAs |
| `run_trrosetta_DB_MSA_job_CPU.sh` | trRNA / cpu | trRosettaRNA on `prepared_db_msa_inputs`, CPU-only |
| `run_rhofold_job.sh`, `run_trRosettaRNA_job.sh`, `run_trRosettaRNA_job_CPU.sh` | — | standalone runners for the two scripts above |

The trRosetta CPU jobs default to `TRROSETTA_USE_JOB_LOCAL_COPY=1`: they copy the
trRosettaRNA tree into the job's scratch first, because concurrent jobs sharing one
installation race on its temporary files.

### Prepared inputs (shipped)

Each folder holds `fasta/`, `a3m/` and `metadata/` (manifest, exported IDs) for the same
91 RNA monomer targets:

| folder | MSA content |
| --- | --- |
| `prepared_shs_inputs_from_rnafold_folder` | SHS from RNAfold secondary structures |
| `prepared_shs_inputs_from_spotrna_folder` | SHS from SPOT-RNA secondary structures |
| `prepared_shs_inputs_from_rnaformerN100_folder` | SHS from RNAformer, N=100 (plus `custom_ss_rnaformer_pred_pairs/`) |
| `prepared_shs_inputs_from_rna_monomers_dssrN100_folder` | SHS from DSSR annotations of the N100 predictions |
| `prepared_shs_inputs_from_rna_monomers_{rnafold,rnaformerN100,spotrna}_denoise_folder` | second-round ("denoised") SHS arms |
| `prepared_db_msa_inputs` | natural BLASTN database MSAs (the non-SHS baseline) |
| `prepared_msa_plus_shs_from_rnaformerN100_folder` | natural MSA + SHS, concatenated |

Prediction outputs and run logs are **not** in git (see `.gitignore`); every
`predictions_*` folder is reproduced by the corresponding job script.

---

## 3. Standard workflow

### Step 1 — inputs

Use a shipped folder (e.g. `prepared_shs_inputs_from_rnaformerN100_folder`), or rebuild
them from the AF3 data files:

```bash
python prepare_shs_inputs_from_af3_json.py \
  --af3-json-dir ../data/datafiles/alphafold3 \
  --out-dir ./prepared_shs_inputs
```

For the natural-MSA arm (needs RhoFold's BLASTN database, §1.2):

```bash
python prepare_standard_msa_from_blastn.py \
  --fasta-dir ./prepared_shs_inputs/fasta \
  --output-a3m-dir ./prepared_standard_msa/a3m \
  --rhofold-root ../../RhoFold
```

### Step 2 — fold

```bash
# RhoFold+  (env: rhofold)
python run_rhofold_shs_batch.py \
  --fasta-dir ./prepared_shs_inputs_from_rnaformerN100_folder/fasta \
  --a3m-dir   ./prepared_shs_inputs_from_rnaformerN100_folder/a3m \
  --msa-mode provided \
  --output-dir ./predictions_rhofold_on_shs_from_rnaformerN100_folder \
  --rhofold-bin ../../RhoFold/inference.py \
  --device cuda --skip-existing

# trRosettaRNA  (env: trRNA)
python run_trrosetta_shs_batch.py \
  --fasta-dir ./prepared_shs_inputs_from_rnaformerN100_folder/fasta \
  --a3m-dir   ./prepared_shs_inputs_from_rnaformerN100_folder/a3m \
  --output-dir ./predictions_trrosetta_on_shs_from_rnaformerN100_folder \
  --trrosetta-root ../../trRosettaRNA_v1.1 \
  --gpu 0 --skip-existing
```

Both take `--ids-file` (one PDB ID per line) and `--skip-existing`, so an interrupted run
resumes where it stopped.

### Step 3 — convert predictions to mmCIF

```bash
python convert_prediction_pdbs_to_eval_cif.py \
  --algorithm rhofold \
  --input-pred-dir ./predictions_rhofold_on_shs_from_rnaformerN100_folder \
  --output-cif-dir ./converted_cif/rhofold_shs_rnaformerN100

python convert_prediction_pdbs_to_eval_cif.py \
  --algorithm trrosetta \
  --input-pred-dir ./predictions_trrosetta_on_shs_from_rnaformerN100_folder \
  --output-cif-dir ./converted_cif/trrosetta_shs_rnaformerN100
```

`--algorithm` only selects which PDB inside a target folder is the model
(`relaxed_*_model.pdb` / `unrelaxed_model.pdb` for RhoFold, `structures/*_model1.pdb` for
trRosetta).

### Step 4 — evaluate

```bash
python evaluate_with_synthevolution.py \
  --algorithm rhofold_shs_rnaformerN100 \
  --pred-cif-dir ./converted_cif/rhofold_shs_rnaformerN100 \
  --eval-type rna_rna \
  --output-dir ./evaluation_results
```

This runs `../run_evaluation.sh <cifs> ../evaluation/predictions/gt rna_rna` in the repo
root and copies `results/pred_rna_rna.csv` / `results/exp_rna_rna.csv` to
`evaluation_results/<label>_*.csv` plus a small mean-metric summary. Because
`run_evaluation.sh` always writes to the same `results/` filenames and rebuilds
`database/rbpDatabase.db`, run one arm at a time.

Equivalent manual call from the repo root:

```bash
./run_evaluation.sh tr_rhofold_shs_eval/converted_cif/rhofold_shs_rnaformerN100 \
                    evaluation/predictions/gt rna_rna
```

### Step 5 — compare

```bash
python compare_algorithms.py \
  --trrosetta-csv ./evaluation_results/trrosetta_shs_rnaformerN100_pred_rna_rna.csv \
  --rhofold-csv   ./evaluation_results/rhofold_shs_rnaformerN100_pred_rna_rna.csv \
  --alphafold-csv ../results/csvs/All_alphafold.csv \
  --output-dir ./comparison

python compare_msa_sources.py \
  --left-csv  ./evaluation_results/rhofold_shs_pred_rna_rna.csv \
  --right-csv ./evaluation_results/rhofold_stdmsa_pred_rna_rna.csv \
  --left-label rhofold_shs --right-label rhofold_stdmsa \
  --output-dir ./comparison_msa_sources/rhofold
```

---

## 4. trRosetta with secondary structure from the base-pair CSV

Instead of letting trRosettaRNA run SPOT-RNA, feed it the base pairs already in
`RNA_Monomers_base_pairs.csv`:

```bash
python prepare_trrosetta_ss_from_base_pairs_csv.py \
  --pairs-csv ./RNA_Monomers_base_pairs.csv \
  --fasta-dir ./prepared_shs_inputs_from_rnaformerN100_folder/fasta \
  --method rnaformer_pred --pairs-column pairs \
  --output-dir ./prepared_custom_ss/rnaformer_pred --output-suffix .ssprob.txt

python run_trrosetta_shs_batch.py \
  --fasta-dir ./prepared_shs_inputs_from_rnaformerN100_folder/fasta \
  --a3m-dir   ./prepared_shs_inputs_from_rnaformerN100_folder/a3m \
  --output-dir ./predictions_trrosetta_csv_ss \
  --trrosetta-root ../../trRosettaRNA_v1.1 --gpu -1 \
  --ss-dir ./prepared_custom_ss/rnaformer_pred \
  --ss-fmt spot_prob --ss-suffix .ssprob.txt --require-ss --skip-existing
```

Notes:

- The converter writes only the upper triangle — trRosetta's `spot_prob` loader
  symmetrises internally (`ss += ss.T`).
- Keep `--ss-fmt spot_prob` for matrices built from pair lists.
- On SLURM the same two stages are one job:
  `sbatch --export=ALL,CSV_METHOD=rnaformer_pred,PAIRS_COLUMN=pairs run_trrosetta_shs_job_CPU_csv_ss.sh`.

---

## 5. Multi-seed SHS sweep

RhoFold+ inference is deterministic (one forward pass plus Amber minimisation), so
re-folding a fixed A3M under different RNG seeds reproduces the same structure. The only
place a seed changes anything is the **synthetic MSA**, and that is what this sweep varies:
each variant's MSAs are regenerated from its base-pair set under seeds 0/42/137 and folded
again.

| variant | CSV method | sequence source |
| --- | --- | --- |
| rnafold | rnafold_pred | prepared_shs_inputs_from_rnafold_folder |
| spotrna | spotrna_pred | prepared_shs_inputs_from_spotrna_folder |
| rnaformerN100 | rnaformer_pred | prepared_shs_inputs_from_rnaformerN100_folder |
| rna_monomers_dssrN100 | dssrN100_dssr | prepared_shs_inputs_from_rna_monomers_dssrN100_folder |
| rna_monomers_rnafold_denoise | rnafoldN100_dssr | prepared_shs_inputs_from_rna_monomers_rnafold_denoise_folder |
| rna_monomers_rnaformerN100_denoise | rnaformerN100_dssr | prepared_shs_inputs_from_rna_monomers_rnaformerN100_denoise_folder |
| rna_monomers_spotrna_denoise | spotrnaN100_dssr | prepared_shs_inputs_from_rna_monomers_spotrna_denoise_folder |

(`python generate_shs_seed_inputs.py --list-variants` prints the same table.)

```bash
bash submit_rhofold_shs_seed_jobs.sh                                        # all 7 variants x 3 seeds
sbatch --export=ALL,VARIANT=spotrna run_rhofold_shs_seeds_job.sh            # one variant
sbatch --export=ALL,VARIANT=spotrna --array=2 run_rhofold_shs_seeds_job.sh  # seed 137 only
```

Each array task writes `prepared_shs_inputs_from_<variant>_seed<SEED>_folder/` (with
`fasta/`, `a3m/`, the generator's `af3_json/` — the same MSAs in AF3 job form — and
`metadata/run_info.json`) and folds it into
`predictions_rhofold_on_shs_from_<variant>_seed<SEED>_folder/`. The shipped single-run
folders are never touched. Array index 0/1/2 selects seed 0/42/137; both stages pass
`--skip-existing`, so a task that hits the wall clock resumes.

Job knobs (`sbatch --export=ALL,VARIANT=...,KNOB=value`): `SEEDS` (default `"0 42 137"`,
keep `--array` in range), `RUN_PREP` / `RUN_FOLD`, `N_MSA` (default 100), `SHS_PARAMS`
(JSON dict of extra generator flags, e.g. `'{"wobble-prob": 0.2}'`), `PREP_DIR` / `OUT_DIR`,
`SYNEVO_PY`, `SHS_GENERATOR`, `RHOFOLD_BIN`, `RHOFOLD_CKPT`, `CONDA_ENV`.

Generation is CPU-only but runs inside the GPU job for convenience. To keep GPUs free,
pre-generate with `EXTRA_EXPORT=RUN_FOLD=0 bash submit_rhofold_shs_seed_jobs.sh` and submit
the folds afterwards with `RUN_PREP=0`.

Locally, one variant/seed without SLURM:

```bash
python generate_shs_seed_inputs.py --variant rnafold --seed 137
python run_rhofold_shs_batch.py \
  --fasta-dir ./prepared_shs_inputs_from_rnafold_seed137_folder/fasta \
  --a3m-dir   ./prepared_shs_inputs_from_rnafold_seed137_folder/a3m \
  --msa-mode provided \
  --output-dir ./predictions_rhofold_on_shs_from_rnafold_seed137_folder \
  --rhofold-bin ../../RhoFold/inference.py --device cuda --skip-existing
```

Targets without a row for the variant's method are skipped and logged in
`metadata/manifest.tsv` (SPOT-RNA has none for 8TKK, DSSR none for 7EFG/7VFT, and 9CXF is
absent from the pairs CSV entirely). Seeded folders are gitignored: ~1.6 MB each and
exactly reproducible from the pairs CSV plus the seed.

**Before spending GPU time**, regenerate one variant under seed 42 and diff it against the
shipped folder — that tells you whether your generator settings match the ones a given
input folder was built with:

```bash
python generate_shs_seed_inputs.py --variant rnafold --seed 42 \
  --out-dir /tmp/seed42check
diff -rq /tmp/seed42check/a3m prepared_shs_inputs_from_rnafold_folder/a3m
```

Note that the shipped `prepared_shs_inputs_from_*_folder` were exported from the paper's
AF3 job files, not produced by this seeding script; a spot check of three targets shows
the same query, the same depth (100) and the same style of alignment, but not byte-identical
rows. Treat the shipped folders as the paper arm and a seed sweep as its own set of
arms — and keep every seed of one comparison generated by the same script and knobs.

---

## 6. Migration notes

This folder came from the separate `shs-eval-pipeline` repository. Changes made while
moving it in:

- Paths that pointed at a sibling SyntheticEvolution checkout now point one level up:
  the SHS generator is `../rna_msa_generator_base_pair.py`, AF3 data files are
  `../data/datafiles/alphafold3`, the evaluator root is `..` and the AF3 baseline is
  `../results/csvs/All_alphafold.csv`. These in-repo defaults resolve from the script's
  own location, so they hold no matter which directory you call the scripts from.
  RhoFold+ and trRosettaRNA are still expected next to the repo (`../../RhoFold`,
  `../../trRosettaRNA_v1.1`), which is unchanged.
- `../rna_msa_generator_base_pair.py` gained an optional `--pairs` flag (a JSON list of
  0-based `[i, j]` pairs) that takes precedence over `--structure` /
  `--structure_predictor`. `generate_shs_seed_inputs.py` needs it to drive the generator
  from `RNA_Monomers_base_pairs.csv`, and unlike dot-bracket it carries pseudoknots and
  multiplets. Nothing else in the generator changed and no existing default moved.
- `convert_prediction_pdbs_to_eval_cif.py`, `evaluate_with_synthevolution.py`,
  `compare_algorithms.py` and `compare_msa_sources.py` were restored from the source
  repository's history — they had been deleted there, which left the pipeline without its
  evaluation stage.
- Prediction outputs (~270 MB) and run logs were not migrated; they are reproduced by the
  job scripts, and this repo distributes predictions through `download.sh` instead.
- Machine-specific strings were scrubbed: the `prefix:` line of both `envs/*.yml` exports
  (conda ignores it when creating an environment) and the absolute cluster workspace paths
  in the `json_file` column of every `prepared_*/metadata/manifest.tsv`, which are now
  written relative to the repo (`../data/datafiles/...`) or to this folder (`./...`).
