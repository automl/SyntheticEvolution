# Storage layout

- **Run:** one experiment with shared configuration and dataset. Start with `submit_standalone_trial.py` or `submit_neps_hpo.py`; inside an existing CPU allocation, use `run_standalone_trial.py` or `run_neps_hpo.py`.
- **Trial:** one parameter configuration evaluated across the entire dataset, producing one aggregate loss. Standalone runs contain one trial; NePS runs evaluate configurations across multiple trials, subject to the evaluation budget.
- **Task:** processing one RNA within a trial, from SHS generation through AF3 and DSSR to scoring. Each task has one AF3 array element.

## Run

```text
<workspace>/shs_hpo/<run_name>/
├── run_metadata.json                 # Configuration, dataset, fingerprint
├── controller_logs/<job_id>.out/.err  # When submitted through a helper
│
│   # Standalone mode:
├── standalone_trial/                 # Pipeline trial files directly here
│
│   # NePS mode:
├── controller.lock                   # NePS controller lock
└── neps/
    ├── <NePS-managed search files>
    └── <NePS-managed trial path>/    # Repeated per trial; see below
```

Standalone and NePS are alternative modes. NePS placeholders indicate ownership, not literal filenames or a verified internal directory structure. `.out/.err` denotes two files.

## One trial

NePS supplies the trial directory; `artifacts/` separates our pipeline outputs from NePS-managed files. Standalone needs no extra level because its directory contains only pipeline outputs.

```text
Standalone:                               NePS:
standalone_trial/                         <NePS-managed trial path>/
└── <pipeline trial files>                ├── <NePS-managed files>
                                          └── artifacts/
                                              └── <pipeline trial files>
```

Both locations contain the same pipeline trial files:

```text
<pipeline trial directory>/
├── trial_metadata.json
├── trial_result.json
├── af3_tasks.json
├── af3_submission_intent.json
├── af3_job_id.json
├── logs/<array_id>_<task_index>.out/.err
└── rna_<index>/                      # Repeated per dataset RNA; see below
```

## One task

Indices follow dataset order: `rna_00000`, `rna_00001`, etc. Original RNA IDs are retained in request and evaluation files.

```text
rna_<index>/
├── shs_generation_request.json
├── shs_generation.log
├── af3_input.base.json
├── af3_input.json
├── base_pair_evaluation.json
├── af3/
│   ├── af3_success.json
│   └── rna_<index>/
│       ├── rna_<index>_model.cif
│       ├── rna_<index>_data.json
│       ├── rna_<index>_confidences.json
│       ├── rna_<index>_summary_confidences.json
│       ├── ranking_scores.csv
│       ├── TERMS_OF_USE.md
│       └── seed-<seed>_sample-<sample>/  # Repeated per AF3 sample
│           ├── model.cif
│           ├── confidences.json
│           └── summary_confidences.json
└── dssr/
    ├── dssr.json
    ├── dssr.log
    └── dssr-*                        # Auxiliary outputs; vary by structure
```

## File purposes

| File | Purpose |
| --- | --- |
| `trial_metadata.json` | Fingerprint and parameter overrides; checked before reuse. |
| `trial_result.json` | Aggregate loss, all task evaluations, parameters, and job ID. |
| `af3_tasks.json` | Ordered input/output/model-directory paths for AF3 array tasks. |
| `af3_submission_intent.json` | Submission command saved before calling `sbatch`. |
| `af3_job_id.json` | Submitted array ID used to reattach. |
| `shs_generation_request.json` | RNA, target pairs, parameter overrides, and seeds. |
| `af3_input.base.json` / `af3_input.json` | Base generator input / final AF3 input with synthetic MSA. |
| `af3_success.json` | Verified top-level model path used for DSSR scoring. |
| `dssr.json` | Nucleotide annotations and predicted pairs used for scoring. |
| `base_pair_evaluation.json` | RNA ID, TP/FP/FN, F1, loss, target/predicted pairs, and model path. Also included in the trial result. |

Logs capture process output. AF3 sample/confidence files and DSSR auxiliary files support inspection. Saved parameters contain overrides, not expanded generator defaults.

## Recovery and export

After metadata checks, an existing trial result returns immediately. Otherwise, a saved job ID skips generation/submission; submission intent without a job ID stops for inspection. Without either, a new array is submitted. Task evaluations are not checkpoints: without a trial result, DSSR/scoring repeats. This does not itself recover the NePS search.

`export_results.py` targets completed **NePS trials only**, writing to `hpo/reports/<run_name>/`:

- `scores.csv`: per-RNA metrics for every completed trial.
- `best.json`: lowest-loss trial result and its source path.
- `run_metadata.json`: copied run metadata.

Raw predictions remain in the workspace; exported model paths still refer there.
