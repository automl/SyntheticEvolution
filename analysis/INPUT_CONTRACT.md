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
│   └── short.csv
└── configs/
    └── config_<id>/
        └── trial_result.json
```

## File ownership

- `run_metadata.json` is created by the HPO pipeline and copied unchanged.
- The three files in `summary/` are created by NePS and copied unchanged.
- Each `trial_result.json` is created by the HPO pipeline and copied unchanged.
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

## Output directory behavior

- By default, analysis output is the sibling directory `analysis/` beside
  `results/`.
- An optional output argument may select another destination.
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
exported configuration IDs, and the paths copied into the snapshot. It is
bookkeeping for the export and is not required to reproduce the pipeline
results.
