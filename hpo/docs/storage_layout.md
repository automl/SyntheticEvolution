# Storage layout

- **Run:** one experiment with one configuration and dataset, in
  `<workspace>/shs_hpo/<run_name>/`. Standalone runs contain one trial; NePS
  runs contain one trial per sampled configuration.
- **Trial:** one parameter configuration evaluated on every RNA of the
  dataset, producing one aggregate loss.
- **Scratch:** a private directory per trial in RAM (`/dev/shm`, setting
  `storage.scratch`) for everything temporary. It is removed when the trial
  ends. The workspace only receives finished files.

The workspace holds as few files as possible: per trial one record and one
archive per RNA. Logs and tool side files stay in scratch.

## Run

```text
<workspace>/shs_hpo/<run_name>/
├── run_metadata.json              # config, dataset, fingerprint (the only absolute paths)
├── controller_logs/               # Slurm output of the submitted controllers
├── STOPPED.json                   # only if the failure limit stopped the run
├── STOPPED-<n>.json               # earlier stops, renamed to restart the run
│
│   # Standalone mode
├── standalone_trial/              # one trial, see below; plus logs/
│
│   # NePS mode
└── neps/
    ├── optimizer_info.yaml, optimizer_state.pkl, pipeline_space.pkl, .trial_cache.pkl
    ├── shared_errors.jsonl        # every crashed trial (NePS), only if one crashed
    ├── summary/                   # full.csv, short.csv, best_config.txt, best_config_trajectory.txt
    └── configs/config_<id>/
        ├── config.yaml, metadata.json, report.yaml   # NePS
        └── artifacts/             # one trial, see below
```

## One trial

```text
Successful trial                     Failed trial
artifacts/                           artifacts/
├── trial_result.json                ├── trial_metadata.json
└── rna_<index>.tar.gz               ├── failure_report.json
    (one per RNA)                    └── rna_<index>.tar.gz
                                         (RNAs whose AF3 prediction finished)
```

A running trial has `trial_metadata.json` (fingerprint and parameters). On
success, `trial_result.json` replaces it. Standalone trials also keep the
tool logs in `standalone_trial/logs/` (`af3.log`, `rna_<index>/shs_generation.log`,
`rna_<index>/dssr.log`).

With `storage.keep_artifacts: trajectory`, successful configs that are
neither on NePS's best-config trajectory nor among the `keep_top_k` lowest
losses keep only `trial_result.json`. Failed trials are never pruned.

### Per-RNA archive

`rna_<index>.tar.gz` holds everything kept from one RNA: the complete AF3
output (without AF3's `<name>/` folder), the generator's effective parameters
and the DSSR output.

```text
rna_<index>.tar.gz
├── rna_<index>_model.cif               # top-ranked model, scored with DSSR
├── rna_<index>_data.json               # AF3 input incl. the synthetic MSA
├── rna_<index>_confidences.json
├── rna_<index>_summary_confidences.json
├── ranking_scores.csv
├── TERMS_OF_USE.md
├── seed-<seed>_sample-<n>/             # one per AF3 sample (5)
│   ├── model.cif
│   ├── confidences.json
│   └── summary_confidences.json
├── shs_effective_parameters.json       # generator parameters after defaults
└── dssr.json                           # added after DSSR; missing if DSSR crashed
```

The archive is published right after AF3 and published again once DSSR has
written `dssr.json`, before its pairs are validated. Its existence marks a
finished AF3 prediction: a resumed trial reuses it instead of running the
generator and AF3 again.

## File purposes

| File | Purpose |
| --- | --- |
| `run_metadata.json` | Run config, dataset and fingerprint; the exporter rebuilds generator requests from it. |
| `trial_metadata.json` | Fingerprint and parameters of a running or failed trial; checked before reuse. |
| `trial_result.json` | Fingerprint, parameters, aggregate loss, per-RNA scores (`model` points into the archive), files and bytes published (`storage`) and, in `trajectory` mode, the `retention` decision. |
| `failure_report.json` | Stage, RNA, error, command, signal, timeout flag, host, Slurm job and the last 64 KB of each log of the failing RNA. |
| `STOPPED.json` | Reason and rates of a failure-limit stop; the run refuses to start until it is renamed. |
| `rna_<index>.tar.gz` | See above. |

## Flow of one trial

1. Generate the AF3 inputs of all RNAs without a readable archive (scratch).
2. Predict all of them with the controller's AF3 worker (scratch), which
   loads and compiles the model once per controller. Publish the archive of
   every finished prediction, also after an AF3 failure or stop.
3. Per RNA: run DSSR in scratch, add `dssr.json` to the archive and score.
4. Write `trial_result.json`, remove `trial_metadata.json`, apply retention.

On a failure, `failure_report.json` is written and the archives already
published stay.

## Export

`export_results.py --config <yaml>` writes a snapshot of the best-config
trajectory into `results/` next to the config:

```text
results/
├── run_metadata.json
├── manifest.json                     # export bookkeeping, incl. a failures block
├── summary/                          # best_config_trajectory.txt, best_config.txt, short.csv, full.csv
├── errors/                           # only if trials failed
│   ├── shared_errors.jsonl
│   └── reports/config_<id>.json      # every failure_report.json
└── configs/config_<id>/
    ├── trial_result.json
    ├── archives/rna_<index>.tar.gz   # only with --include-archives
    └── shs/rna_<index>/
        ├── af3_input.json            # the archive's rna_<index>_data.json
        └── shs_generation_request.json   # rebuilt from run_metadata.json
```

`--prune [--dry-run]` instead deletes the archives of configs that are
neither on the trajectory nor in the top `keep_top_k`, also for a run that
used `keep_artifacts: all`.
