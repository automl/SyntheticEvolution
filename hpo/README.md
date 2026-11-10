# SHS hyperparameter search on Helix

Searches the parameters of the synthetic MSA generator (SHS) with NePS. Each
trial generates an MSA for every RNA of a dataset, predicts the structure with
AlphaFold 3, extracts the base pairs with DSSR and scores them against the
target pairs. The trial loss is the mean per-RNA loss.

```text
NePS samples a config → SHS generator → AlphaFold 3 → x3dna-dssr → F1 / WL loss → NePS
```

## Entry points

Run from the repository root with the controller environment
(`~/.conda/envs/controller/bin/python`).

| Command | What it does |
| --- | --- |
| `python -m hpo.submit_standalone_trial --config <yaml>` | Submit one controller job that evaluates the `fixed` parameters once, or every trial in `trials`. |
| `python -m hpo.submit_neps_hpo --config <yaml>` | Submit `controller.count` NePS controllers (a Slurm array). |
| `python -m hpo.run_standalone_trial --config <yaml> [--validate-only]` | Run a standalone run inside an existing allocation; `--validate-only` only prints the parsed dataset. |
| `python -m hpo.run_neps_hpo --config <yaml>` | Run one NePS controller inside an existing allocation. |
| `python -m hpo.export_results --config <yaml> [--include-archives]` | Copy the trajectory configs (check run: every finished trial) into `results/` next to the config, for `analysis/`. |
| `python -m hpo.export_results --config <yaml> --prune [--dry-run]` | After a run: delete the archives of configs that are neither on the trajectory nor in the top `keep_top_k`. |
| `python -m hpo.make_checks <results> [--heldout <csv>]` | Locally, after an export: write a check run for the best configs (see below). |

The submit scripts check the config and the environment (interpreters and
their packages, DSSR, AF3 models and module) and report every problem before
anything waits in the queue. They ask before resuming an existing run;
`submit_neps_hpo` also refuses a run stopped by the failure limit. Jobs run
`hpo/slurm/controller.slurm`, which loads the AF3 module and starts the
controller. Analysis of an export: `python -m analysis.analyze_results <results>`.

## Configuration

Start from [examples/helix.example.yaml](examples/helix.example.yaml); every
setting is commented there.

- **`run_name`:** the run directory `<workspace>/shs_hpo/<run_name>/`. The
  fingerprint covers the dataset, scientific settings and the pipeline and
  generator code; if one of them changes, a run cannot be resumed and needs a
  new name. Resources, timeouts, `evaluations`, `ignore_errors`, `storage`
  and `failures` may change.
- **`evaluations`:** NePS counts it **per controller**: a run with
  `controller.count: 4` evaluates up to 4 × `evaluations` configs, and every
  resubmission adds that again. Leave it out to evaluate as many configs as fit
  into `controller.time`.
- **`fixed` / `search`:** generator parameters; NePS samples `search`,
  standalone trials use `fixed` and the generator defaults.
- **`trials`** (standalone only, optional): named trials that run one after
  another in one job with one AF3 worker, each with `parameters` over `fixed`
  and optionally its own `shs_seed` / `af3_seed`. Other keys are kept for
  analysis. Not part of the fingerprint: trials can be added and the run
  resubmitted. A failed trial does not stop the others.
- **`ignore_errors: true`:** NePS keeps sampling when a trial fails. A failed
  trial gets the worst loss (1.0), so the optimizer avoids its region.
- **`failures`:** the failure limit (default: stop when more than 20% of the
  last 50 finished trials failed, after at least 5).
- **`storage`:** `scratch` (default `/dev/shm`), `keep_artifacts`
  (`all`, the default, or `trajectory`) and `keep_top_k` (default 5).
- **`scoring_objective`:** `f1` (default), `f1_soft` (F1 with partial credit
  for pairs slipped by 1–3 nucleotides: 0.75 / 0.5 / 0.25) or `wl`
  (Weisfeiler-Lehman graph kernel, `wl_iterations`).
- **`af3_timeout_seconds`:** per RNA. One AF3 worker per controller loads and
  compiles the model once and predicts the RNAs of every trial.
- **`af3_buckets`:** optional token sizes AF3 pads to (default: AF3's, from
  256). Small buckets speed up short RNAs but change the predictions like a
  different seed, so they are part of the fingerprint; every bucket the
  dataset uses is compiled once per controller (~25 s).

## Check runs

A check run tests whether the best configs of an HPO run hold up. Its loss is
the minimum of many single-seed losses and therefore optimistic.
`make_checks` reads an export and writes `checks/config.yaml` (a standalone
config with `trials`) and `checks/dataset/` next to `results/`:

- `c<id>_hpo`: each of the top `--top-n` configs (default 3) with the HPO seeds;
  reproduces the HPO loss with the same GPU type, AF3 module and buckets,
- `c<id>_s<seed>`: the same configs with `--seeds` new seed pairs (default 5),
- `nomsa_s<seed>`: baseline `N: 1` (only the query sequence, i.e. AF3 without MSA).

The dataset holds the HPO training rows (`1_train.csv`, from
`run_metadata.json`) and the optional held-out RNAs (`2_heldout.csv`). The
`input` written into the config is the Helix path `checks/dataset` next to the
HPO input (`--input` overrides it). Commit the check directory to the thesis
repo, pull it on Helix and submit the config with `submit_standalone_trial`.

## Dataset

A CSV file, or a directory whose `*.csv` files are combined in name order.
Columns: `id` (unique), `sequence` (A, C, G, U) and at least one of

- `pairs`: the scoring targets as a JSON list of zero-based position pairs,
  e.g. `"[[0,14],[1,13]]"`; `[]` for no pairs. Reversed pairs and duplicates
  are normalized.
- `interactions`: `[[i, j, strength], ...]` with strengths in [0, 1], used by
  the generator. Without `pairs`, the pairs with a positive strength are the
  scoring targets.

Optional `mutation_rates`: one rate per position. Examples:
[examples/example_input.csv](examples/example_input.csv),
[examples/weighted_input.csv](examples/weighted_input.csv).

## Scoring

DSSR's full `pairs` list is used, including noncanonical pairs, mapped through
its nucleotide records (single chain A, residues 1..N, matching sequence).
F1 loss: `1 - 2TP / (2TP + FP + FN)`; both empty gives 0, exactly one empty
gives 1. Every RNA has equal weight. Tool or parser failures fail the trial;
they never become an empty prediction.

## Storage

The workspace gets one archive per RNA (complete AF3 output incl. all samples,
the MSA, the generator's effective parameters and `dssr.json`) and one trial
record; logs and temporary files stay in a per-trial scratch directory in RAM.
Details: [docs/storage_layout.md](docs/storage_layout.md).

- **`keep_artifacts: all`** keeps the archives of every config. Use it for
  scientific runs.
- **`keep_artifacts: trajectory`** keeps them only for configs on NePS's
  best-config trajectory or in the top `keep_top_k`; the others are pruned
  during the run and keep their `trial_result.json`. Failed trials are always
  kept.
- `/dev/shm` counts towards `controller.memory`.

## Failures, stops and restarts

- **Failed trial:** `artifacts/failure_report.json` holds the stage, RNA,
  error, signal and the end of each log. The controller logs one `ERROR` line
  per failure and a run summary at the end; the export copies all reports to
  `results/errors/` and warns.
- **Failure limit:** the first controller that sees the limit exceeded writes
  `STOPPED.json`; every controller reports its next trial as stopped and
  exits with code 3. To restart, rename the file as printed (`STOPPED-<n>.json`);
  the limit then counts only trials that finish later.
- **Time limit:** a controller does not start a trial that would not finish
  before its time limit (judged by its longest trial so far); it reports the
  sampled trial as stopped and exits with code 0.
- **Stop signal and `scancel`:** jobs get `SIGUSR1` shortly before their time
  limit (`signals.STOP_MARGIN_SECONDS`, 2 minutes). A trial still running then
  is stopped and reported as failed, and the controller exits with code 4.
  Finished AF3 predictions are kept.
- **Resume:** submitting the same config again resumes the run. A standalone
  run skips finished trials, reuses archives and only reruns what is missing; NePS never reruns
  a finished or failed trial.

## Tests

```bash
~/.conda/envs/synEvo/bin/python -m pytest hpo/tests
```

The tests need no cluster: Slurm, AF3 and DSSR are replaced. The generator CLI
tests run the real generator.
