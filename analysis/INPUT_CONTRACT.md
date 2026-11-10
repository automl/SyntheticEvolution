# Analysis input contract

Analysis reads the `results/` directory created by
`hpo/export_results.py`. The exporter copies a selected HPO run from the
cluster workspace; analysis must not read the cluster workspace directly.

## Current layout

```text
results/
├── run_metadata.json
├── manifest.json
├── summary/
│   ├── best_config_trajectory.txt
│   ├── best_config.txt
│   ├── short.csv
│   └── full.csv
├── errors/                               # optional: only if trials failed
│   ├── shared_errors.jsonl
│   └── reports/config_<id>.json
└── configs/
    └── config_<id>/
        ├── trial_result.json
        ├── archives/rna_<index>.tar.gz   # optional: export --include-archives
        └── shs/
            └── rna_<index>/
                ├── af3_input.json
                └── shs_generation_request.json
```

## File ownership

- `run_metadata.json` is created by the HPO pipeline and copied unchanged.
- The four files in `summary/` are created by NePS and copied unchanged.
- Each `trial_result.json` is created by the HPO pipeline and copied unchanged.
- The SHS files are derived by the exporter for each exported trajectory
  configuration: `af3_input.json` is AF3's `<name>_data.json` from the RNA's
  archive (the AF3 input with the synthetic MSA), and
  `shs_generation_request.json` is rebuilt from `run_metadata.json` exactly as
  the pipeline built it for the generator.
- The optional `errors/` files are NePS's crash log and the pipeline's failure
  reports, copied unchanged. The optional archives hold the complete AF3
  output (all samples), the generator's effective parameters and `dssr.json`;
  see `hpo/docs/storage_layout.md`.
- `manifest.json` is created by `hpo/export_results.py`; it is not produced by
  the HPO pipeline.

## Required relationships

- `best_config_trajectory.txt` contains configuration IDs on lines matching
  `Config ID: <id>`.
- Every ID in that file has one corresponding
  `configs/config_<id>/trial_result.json`.
- Configuration IDs are integers, and their order in the trajectory is
  meaningful.
- Only configurations listed in the trajectory are exported.

## Analysis rules

- Treat `results/` as read-only input.
- Write generated reports, tables, and plots to a separate `analysis/` output
  directory next to `results/` by default. An optional output path may override
  this default.
- Missing required files or trajectory trial results are errors.
- Additional files may be ignored unless an analysis explicitly uses them.
- Do not assume fields beyond those needed by a specific analysis. The
  structure of `trial_result.json` and the summary files remains intentionally
  open for now.

## Analysis orchestration

- One main analysis entry point runs all registered analysis modules.
- There is no per-module selection interface in the initial implementation.
- Each module receives shared run data and may write only inside its assigned
  output directory.
- Modules should be independent so new optimization, trial, and later sample
  analyses can be added without redesigning the input loader.
- Trial results are exposed as plain dictionaries. Analyses must tolerate
  additional keys and metrics.
- Individual samples are not part of the initial analysis scope. The structure
  should allow later modules to inspect selected configurations or the top `N`
  configurations when the exporter provides the required artifacts.
- The trajectory dot-bracket analysis uses each exported
  `configs/config_<id>/trial_result.json` and writes its combined output under
  the module's analysis directory.
- The trial comparison analysis compares `parameters` and the per-RNA
  `scores[].predicted_pairs` of all trajectory configurations and writes
  `trial_comparison.txt`. It is skipped when the trajectory has only one
  configuration.
- The optimization history analysis reads `summary/full.csv` (every trial)
  and writes `loss_over_time.png` (loss and best-so-far over wall-clock time,
  crashed trials as ticks), `loss_histogram.png` and `parameters_vs_loss.png`.
  It is skipped with a warning when `full.csv` is missing (older exports).
- The per-RNA analysis plots each trajectory configuration's per-RNA
  `scores[].loss` as `per_rna_heatmap.png` and `per_rna_lines.png` (one line
  per RNA chain, plus the trial loss).
- SHS analysis uses the exported `shs/rna_<index>/` files and never reads the
  original HPO workspace.

## Output directory behavior

- By default, analysis output is the sibling directory `analysis/` beside
  `results/`.
- An optional output argument may select another destination.
- An optional JSON analysis configuration can be passed with `--config`.
- The SHS module reads `shs.top_n` from that file; when omitted, it analyzes
  only the authoritative best configuration.
- If the output directory already exists, the program asks whether it should be
  recreated.
- An empty response, `y`, or `yes` recreates the directory.
- Any other response stops the program without running analyses.
- Generated output may be recreated; the exported `results/` input is never
  modified.

## Failure behavior

- The orchestrator logs each major step and each module start and completion.
- If a module fails, the error is logged and the remaining modules continue.
- After all modules finish, the orchestrator raises an error if any module
  failed.

## Best configuration

`summary/best_config.txt` is the authoritative NePS best-configuration summary.
The initial analysis code may use it directly rather than reconstructing the
best configuration from trial losses or other summary files.

## Manifest

The manifest records the source config path, source workspace, export time,
exported configuration IDs, a `failures` block, and the paths copied into the
snapshot. The `failures` block is always present: `status` is `none` or
`found`, with the count, rate, groups by error type, signal, stage and RNA,
and the copied error files. It is
bookkeeping for the export and is not required to reproduce the pipeline
results.

## Check runs

`python -m analysis.analyze_checks <checks/results> [--hpo-results <results>]`
reads the `results/` of a check run:
a standalone run with `trials`, written by `hpo/make_checks.py` and exported
by `hpo/export_results.py`.

```text
checks/results/
├── run_metadata.json                 # config (incl. trials) and dataset
├── manifest.json                     # `trials`: the config's trial list; failures block
├── errors/reports/<name>.json        # optional: failed trials
└── trials/<name>/
    ├── trial_result.json             # parameters, shs_seed, af3_seed, loss, scores
    ├── archives/rna_<index>.tar.gz   # optional: export --include-archives
    └── shs/rna_<index>/{af3_input.json, shs_generation_request.json}
```

- The trial list is `manifest.json`'s `trials`, not `run_metadata.json`'s
  config: trials may be added after the run started. Each entry has `name`
  and `parameters`; trials written by `make_checks` also have `check`
  (`reproduction`, `seed` or `baseline`), `base` (the HPO config id) and the seeds.
- Trials listed without a `trials/<name>/` directory have failed or not run;
  failed ones have an error report.
- Training vs held-out RNAs: an RNA is *train* if it is in the HPO run's
  dataset (the HPO export's `run_metadata.json`, matched by id and sequence),
  otherwise *held-out*. The analysis therefore also reads the HPO export.
- `--hpo-results` defaults to `results/` two levels up
  (`<version>/checks/results` → `<version>/results`); it needs
  `run_metadata.json` and `summary/full.csv` (HPO loss per config id).
- Groups: `c<base>` for trials with a `base`, `no MSA` for `check: baseline`,
  otherwise the trial name. Trials with `check: reproduction` are not part of
  the seed statistics.
- Output (default `analysis/` next to `results/`, recreated after asking):
  `trials.csv` (loss on all, train and held-out RNAs per trial),
  `summary.csv` (per group: HPO loss, reproduction and its difference to the
  HPO loss, mean and SD over seeds on train and held-out RNAs),
  `stability.csv` (per group and RNA: mean pairwise F1 between the seeds'
  predicted pairs, loss mean and SD), `predictions.txt` (per RNA: target and
  every trial's prediction as dot-bracket with its loss; pairs dot-bracket
  cannot show are counted), `seed_spread.png`, `per_rna_loss.png`
  (per RNA: each group's seed mean and range, baseline as ring) and `stability.png`.
