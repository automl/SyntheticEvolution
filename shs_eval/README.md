# SHS evaluation pipeline

Compares SHS generation approaches by folding their alignments with AlphaFold 3 and reading the
resulting base pairs back out with DSSR. Four resumable stages; re-running any stage only does the
work that is still missing.

```
plan  ->  fold  ->  annotate  ->  score
```

## Setup

Copy `config.example.yaml` and edit the two commands to match your cluster. Those templates are
the only site-specific part; nothing else needs changing.

```bash
cp shs_eval/config.example.yaml shs_eval/config.yaml
$EDITOR shs_eval/config.yaml     # set the AF3 and DSSR invocations
```

Check the rendering before spending GPU time:

```bash
python -m shs_eval.run plan --config shs_eval/config.yaml --limit 4
python -m shs_eval.run fold --config shs_eval/config.yaml --limit 4 --dry-run
```

`--dry-run` prints the exact commands (or writes the sbatch script) without executing anything.

## Running

```bash
python -m shs_eval.run all \
  --config shs_eval/config.yaml \
  --work-dir results/shs_eval \
  --approaches potts watson_crick_cov \
  --shs-seeds 1 2 \
  --input-json-dir data/datafiles/alphafold3 \
  --gt-dssr-dir evaluation/dssr/gt
```

`--input-json-dir` reuses the original AF3 data files so protein chains and modifications survive
into the run; without it, RNA-only inputs are synthesised from the target pickle.
`--gt-dssr-dir` enables the accuracy columns; without it only loop closure is reported.

Anything after `--generator-flags` is passed straight through to `shs_generator.py`:

```bash
  --generator-flags --potts-coupling 4.0 --potts-wobble 0.9 -N 100
```

## On SLURM

Set `scheduler.kind: slurm` and fill in `sbatch_args` and `setup`. The fold stage then writes
`af3_tasks.txt` plus one or more `af3_array*.sbatch` and submits them, instead of running
serially.

Three things the config must get right, because the compute node shares almost nothing with your
login shell:

- **`setup`** — the module loads and environment activation the node needs. Without it the AF3
  command runs against the node's bare default environment.
- **`sbatch_args`** — partition, GPU, memory and walltime. The defaults in the example are
  placeholders, not recommendations.
- **`max_array_size`** — arrays larger than the site's `MaxArraySize` are rejected outright.
  Check with `scontrol show config | grep MaxArraySize`; work above the limit is split across
  several submissions automatically.

Submission returns immediately, so `all` stops after `fold` and prints the two commands to run
once the array drains:

```bash
python -m shs_eval.run annotate --config shs_eval/config.yaml --work-dir results/shs_eval
python -m shs_eval.run score    --config shs_eval/config.yaml --work-dir results/shs_eval
```

Both are resumable, so partial arrays are fine — run them, look at the counts, resubmit `fold`
for whatever is still missing (it re-plans only the jobs with no model CIF).

Note that `annotate` runs DSSR serially in the calling process. That is seconds per job and fine
for a few thousand, but it is still work on a login node; run it inside a small allocation if your
site is strict about that.

## Output

`<work-dir>/scores.csv`, one row per (target, approach, SHS seed, AF3 seed), plus a summary table
per approach. Two metric families:

- `closure_*` — predicted pairs against the **seeded** pairs. *Did the alignment transmit the
  intended constraint into AF3?* Available for every target, no ground truth needed.
- `accuracy_*` — predicted pairs against **ground-truth** DSSR, plus `input_f1` for how good the
  seeded structure was in the first place. Requires `--gt-dssr-dir`.

Report accuracy as the result. Loop closure is a diagnostic and a training signal, not a headline
number — it is high whenever AF3 obeys the input, including when the input is wrong.

## Layout

```
<work-dir>/
  manifest.jsonl     one line per job; edit or filter it to resubmit a subset
  inputs/            AF3 input JSONs with the generated SHS
  af3/<job>/         AF3 output
  dssr/<job>.json    DSSR annotations
  logs/              per-job stdout/stderr
  scores.csv
```

## Smoke test

`tests/mock_cluster.yaml` wires the pipeline to stand-ins that fabricate AF3 and DSSR outputs, so
the whole flow can be exercised without a cluster:

```bash
python -m shs_eval.run all --config shs_eval/tests/mock_cluster.yaml \
  --work-dir /tmp/shs_eval_smoke --limit 8
```

The mocks test plumbing only. `mock_dssr` derives its pairs from the query sequence alone, so
**every approach scores identically under the mock** — that is the expected result there, not a
finding about the approaches.

## Two things that will bite

**Skip the genetic search.** Pass `--norun_data_pipeline` (or your build's equivalent) in the AF3
command. The SHS is already in the input JSON, so without that flag AF3 spends most of its wall
clock on an HMMER search whose RNA output it will not use.

**Use at least two AF3 seeds.** AF3's diffusion sampling is stochastic. A single-seed comparison
cannot separate a real difference between approaches from sampling noise, and the seed variance
has not yet been measured for this benchmark.
