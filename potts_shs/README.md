# Potts-model SHS: generator, calibration, and evaluation

Everything on the `learning-shs` branch in one place: the new Potts substitution model for
synthetic homolog sets, how it was calibrated against Rfam, the evaluation pipeline built to
test it, and what that evaluation found.

```
potts_shs/
  README.md              this file
  potts.py               the Potts substitution model
  calibrate_rfam.py      fits coupling and wobble to Rfam alignments
  validate_potts.py      offline rate-matching and covariation checks
  eval/                  the four-stage evaluation pipeline
    run.py config.py dssr.py metrics.py
    tests/               AF3 and DSSR stand-ins for the smoke test
  configs/               cluster configs (example, Helix array, Helix serial)
  jobs/                  SLURM scripts for the full evaluation
  tests/                 synthetic Rfam fixture for calibrate_rfam
  runs/
    rna_only/            the 91-target evaluation (manifest, inputs, af3, dssr, scores.csv)
    calibration/         potts_calibration.json, potts_calibration_refined.json
```

The pair representation the model reads (`pair_map.py`) and the generator it plugs into
(`shs_generator.py`) stay in `SHS-Generator/`, since every approach uses them; `potts.py` is
reached from there via `--pair-mutation-approach potts`.

**Bottom line: the Potts model transmits a seeded structure into AlphaFold 3 more faithfully
than the existing approach, but that does not improve structure recovery. Covariation itself is
worth ~+0.058 MCC; which covariation model produces it makes no measurable difference.**

---

## 1. What was built

| component | file | what it does |
|---|---|---|
| Potts substitution model | `potts_shs/potts.py` | samples a row's substitutions jointly from a Potts model |
| Generator integration | `SHS-Generator/shs_generator.py` | `--pair-mutation-approach potts` plus its parameters |
| Pair representation | `SHS-Generator/pair_map.py` | matrix form; interaction strengths and per-position rates |
| Offline validation | `potts_shs/validate_potts.py` | rate matching and covariation readout on toy structures |
| Rfam calibration | `potts_shs/calibrate_rfam.py` | fits coupling and wobble to real alignment statistics |
| Evaluation pipeline | `potts_shs/eval/` | plan → fold (AF3) → annotate (DSSR) → score |

### The model

For a coupled component of positions, the distribution over substitutions is

```
P(x) ∝ ∏ f_i(x_i) · exp( γ · Σ p_ij · s(x_i, x_j) )
```

`f_i` puts weight on the query base (from that position's mutation rate), `s` scores canonical
pairs 1 and wobbles `--potts-wobble`, `p_ij` is the interaction strength from the pair map, and
γ is `--potts-coupling`. Indels are untouched — the generator keeps its own insertion/deletion
layer and the Potts draw only replaces the substitution layer.

Components up to `--potts-max-enumerate` (7) are tabulated exactly; larger ones fall back to
Gibbs sampling. Fields are solved so realised mutation rates match the requested ones for
components up to `--potts-max-rate-match` (4).

On the 919-row protein-RNA set, 71% of targets have all components ≤4 (exact + rate-matched),
14% are 5–7 (exact, rates undershoot), 8% need Gibbs, largest component seen 29. Measured
realised rates are within 0.005 of each other across all three strata, so the fallbacks are
**not** a confound.

---

## 2. Calibration against Rfam

`calibrate_rfam.py` matches covariation statistics against Rfam seed alignments at the depth the
generator is actually used at.

```bash
python potts_shs/calibrate_rfam.py \
  --seed-file evaluation/data/rfam/Rfam.seed.gz --depth 100 \
  --gamma-grid 2.0 2.25 2.5 2.75 3.0 3.5 --wobble-grid 0.3 0.4 0.5 0.6 0.7 \
  --out potts_shs/runs/calibration/potts_calibration_refined.json
```

Timing: ~5 s to scan the archive, then ~15 s per grid point over all qualifying families. The
default 9×6 grid is ~14 minutes on one core. Only 164 of 4227 families qualify at depth 100
(706 at depth 20, 314 at 50), so `--max-families 300` is never reached.

**The default grid is too coarse.** It steps γ 1→2→3, and the Rfam MI target falls between 2
and 3. The coarse grid returned γ=3.0/w=0.70 (score 0.1895); refining to a 0.25 grid returned
**γ=2.25, w=0.40 (score 0.0161)** — a 12× better fit, with both parameters interior.

| | MI@pairs | wobble fraction | canonical |
|---|---|---|---|
| Rfam target | 0.5519 | 0.0896 | 0.8518 |
| γ=2.25, w=0.40 | 0.5369 | 0.0890 | 0.7447 |
| γ=6.0, w=0.85 *(old default)* | ~1.26 | 0.0905 | 0.8772 |

Calibrated values, applied to **every** arm so mutation rate is never confounded with the model:

```
--potts-coupling 2.25 --potts-wobble 0.40
--mutation-rate-paired 0.345 --mutation-rate-unpaired 0.366
```

### Two caveats on the calibration

**MI and pairing fidelity cannot both be matched.** At the MI optimum, 16% of seeded pairs are
non-pairing combinations against 4.8% in Rfam; at γ≈6 composition matches but MI is twice too
high. `discrepancy()` scores only MI and wobble, so it optimises into that corner. The evaluation
tested both settings and found no difference, which closes the question empirically.

**The Rfam MI target is inflated by ~20%.** 22.3% of seeded-pair columns in real alignments are
gaps or `N`, and the MI is computed over a 6-symbol alphabet including them, while the Potts
sampler emits only ACGU. Restricting to gap-free rows drops the target from 0.5519 to 0.4392.

**Leakage:** the 164 calibration families were not filtered against the benchmark. Harmless for
ranking arms (all share γ), but not for a standalone "calibrated to Rfam" claim. Excluding them
needs a cmscan of the queries against the already-`cmpress`ed `evaluation/data/rfam/Rfam.cm`.

---

## 3. Evaluation pipeline (`potts_shs/eval/`)

Four resumable stages; re-running any stage only does the work still missing.

```
plan  →  fold (AlphaFold 3)  →  annotate (DSSR)  →  score
```

AF3 and DSSR are drop-in commands rendered from templates, so no site paths are baked into the
code. `potts_shs/configs/config.helix.yaml` (arrays) and `config.helix.single.yaml` (serial) are filled
in and verified for Helix.

### Running it

```bash
python -m potts_shs.eval.run plan --config potts_shs/configs/config.helix.single.yaml \
  --work-dir potts_shs/runs/rna_only \
  --targets-pickle data/rnaformer_predictions.pkl \
  --target-ids potts_shs/runs/rna_only/rna_only_targets.txt \
  --input-json-dir data/datafiles/alphafold3 --shs-seeds 1 2 \
  --approaches potts --arm-label potts_g225 \
  --generator-flags --potts-coupling 2.25 --potts-wobble 0.40 \
                    --mutation-rate-paired 0.345 --mutation-rate-unpaired 0.366 -N 100

python -m potts_shs.eval.run annotate --config <cfg> --work-dir <wd> --gt-cif-dir evaluation/predictions/gt
python -m potts_shs.eval.run score    --config <cfg> --work-dir <wd> --canonical-only
```

One `plan` pass per arm; `--arm-label` names the arm and the manifest merges across passes.
`--arms` narrows fold/annotate/score to one arm so arms can run as parallel jobs against a
shared work-dir.

### Layout

```
<work-dir>/
  manifest.jsonl     one line per job; filter it to resubmit a subset
  inputs/            AF3 input JSONs with the generated SHS
  af3/<job>/         AF3 output
  dssr/<job>.json    DSSR annotations
  dssr_gt/           ground-truth annotations, from --gt-cif-dir
  dssr_scratch/      per-job working directories for DSSR's aux files
  logs/              per-job stdout/stderr
  scores.csv
```

### Metrics

- `closure_*` — predicted pairs against the **seeded** pairs. Did the alignment transmit the
  intended constraint into AF3? Available for every target.
- `accuracy_*` — predicted pairs against **ground-truth** DSSR, plus `input_f1` for how good the
  seeded structure was. Requires `--gt-cif-dir`.

`--canonical-only` restricts both sides to Watson-Crick and wobble, matching the protocol behind
the existing AF3 and RhoFold MCC numbers.

---

## 4. Benchmark

**91 RNA-only targets** — the entries in `data/datafiles/alphafold3` with no protein chain
(`potts_shs/runs/rna_only/rna_only_targets.txt`). All 91 have ground-truth structures in
`evaluation/predictions/gt/`, and RNA-only inputs carry **no templates**, so there is no
self-templating leakage.

Pairs come from `data/rnaformer_predictions.pkl`, one row per chain; identical homodimer copies
are deduplicated to one row per (target, sequence), giving exactly 91.

**The seeded pairs are RNAformer predictions, not ground truth.** `get_secondary_structure_prediction.py`
uses this table as a prediction cache. So the pipeline under test is

```
RNAformer 2D prediction → SHS MSA → AF3 3D → DSSR 2D → compare to deposited structure
```

`closure_*` measures obedience to a *prediction*; `input_f1` ≈ 0.85 is RNAformer's own accuracy.

### Why not the 919-row protein-RNA set

It was tried first and abandoned: ground truth covered only 358 of 742 targets, and **~15% of
targets carried their own PDB entry as a protein template**, handing AF3 the deposited structure.
Both problems vanish with the RNA-only set.

---

## 5. Results

91 targets × 5 arms × 2 SHS seeds × 2 AF3 seeds = 1,820 folds, 5.4–6 h per arm on one A100,
all completed.

### Arms

```
nomsa  →  none  →  watson_crick_cov / potts_g225  →  potts_g6
no MSA    depth      real-level covariation          over-coupled
```

`none` still builds a full depth-100 MSA — only covariation is removed. `nomsa` is depth-1
(`-N 1`), AF3 from its own prior, with no protein MSAs either on this target set.

### Headline: covariation works; the Potts model does not beat the incumbent

Paired per-target, n=83, accuracy MCC, canonical-only:

| arm | acc MCC | median | vs `none` | p |
|---|---|---|---|---|
| `watson_crick_cov` | **0.9472** | 0.9658 | +0.0590 | 1.2e-06 |
| `potts_g225` | 0.9459 | 0.9670 | +0.0577 | 1.4e-06 |
| `potts_g6` | 0.9433 | 0.9670 | +0.0555 | — |
| `none` | 0.8882 | 0.9442 | — | — |
| `nomsa` | 0.8861 | 0.9541 | — | — |

**Covariation is worth ~+0.058 MCC** over a depth-matched control (63 wins / 11 losses),
reproducing the earlier ablation's +0.075–0.080 on a cleaner benchmark.

**The three covariation arms are statistically indistinguishable.** `watson_crick_cov − potts_g225`
= +0.0013, p=0.74, with **44 of 83 targets exact ties**. A genuine null result on the main
question.

### Two findings worth more than the null

**1. MSA depth alone is worthless.** `none` vs `nomsa` is +0.0021 mean, and by rank `nomsa`
*wins* on 47 of 83 targets (p=0.03). 100 sequences without covariation are no better than
nothing. The entire benefit of the SHS is covariation, not depth.

**2. The pipeline improves on its own input.** `input_f1` = 0.8476 (RNAformer vs deposited), but
output `accuracy_f1` = 0.9462 — **+0.10 F1**. AF3 corrects RNAformer's errors rather than merely
transmitting them.

### Where Potts wins — and why it doesn't matter

```
potts_g6   closure_mcc 0.8760   vs watson_crick_cov  +0.0044  p=0.014
potts_g225 closure_mcc 0.8741   vs watson_crick_cov  +0.0025  p=0.093
```

Potts MSAs transmit the seeded structure **more faithfully** — significant for γ=6 — but that
fidelity does not convert into accuracy, because the input is only ~0.85 accurate. Following an
imperfect structure more closely buys nothing.

### Reliability

SHS-seed sensitivity is low across all arms (|Δ| ≤ 0.010; ≤0.003 for the covariation arms), so
the null is not noise-driven. Both γ settings behave the same, which closes the calibration
tension empirically.

### Caveats

- 83 of 91 targets; 4 have ground truth that cannot be offset-aligned, 4 lack canonical pairs.
- Seeded pairs are RNAformer predictions, so closure measures obedience to a prediction.
- Potts calibration used Rfam families overlapping the benchmark — harmless here, since all arms tie.

Results: `potts_shs/runs/rna_only/scores.csv` (canonical-only) and `scores_all_pairs.csv`.

---

## 6. Things that will bite

Every one of these cost a run or produced silently wrong numbers.

**DSSR chain keys.** DSSR keys `chains` as `m1_chain_A` but reports `chain_name: "A"` on every
nucleotide. Comparing them directly matches nothing and yields an **empty pair list with no
error** — every metric reads 0. The mocks reproduce the real key shape so this cannot hide again.

**Ground truth rarely has the query's length.** Deposited structures carry an extra 5′ G, a
modified residue, a disordered terminus. A strict equal-length guard dropped 51 of 91 targets;
87 of 91 align by a single offset. `align_reference()` remaps indices; without it more than half
the benchmark is silently unscored.

**Chain selection.** In a complex holding two different RNAs, "the longest chain" is the wrong
chain. Predictions are matched by the AF3 chain ids the SHS was written into; ground truth is
matched by **sequence**, because deposited chains use the depositor's labelling (`m1_chain_Q`).

**One pickle row is one RNA chain, not one PDB entry.** 919 rows collapse onto 742 four-character
ids; colliding rows overwrite each other's inputs. Rows after the first are keyed `1FEU_c2`.

**Relative paths break containers.** Singularity requires bind *destinations* to be absolute. A
relative `--work-dir` renders `--bind results/...:results/...` and every container aborts. The
work-dir is now resolved; `stage_fold` also returns non-zero when *every* fold fails, which it
previously did not — all five arms once reported COMPLETED having produced nothing.

**DSSR litters its working directory** with fixed-name aux files (`dssr-pairs.pdb`, …), so
concurrent runs collide and drop targets. Each job gets a private cwd — which means paths in the
command templates must be absolute; use `{root}` for anything in the repository. DSSR also
mangles paths beyond ~400 characters and can abort, so keep `--work-dir` short.

**The QOS submit cap, not `MaxArraySize`.** QOS `normal` allows 1500 submitted jobs; 9,190
one-fold array tasks are rejected outright (`QOSMaxSubmitJobPerUserLimit`). `jobs_per_task`
batches folds per array task. `gres/gpu=100` is the concurrency allowance, so `array_throttle: 8`
wastes most of it.

**Empty references are undefined, not zero.** Targets with no seeded pairs would score a hard 0
and drag every mean down; `pair_metrics` returns NaN so they are excluded from summaries.

**Skip the genetic search.** `--norun_data_pipeline` is what makes this cheap — the SHS is
already in the input JSON.

**Use at least two AF3 seeds.** Diffusion sampling is stochastic, and a single seed cannot
separate a real difference from noise.

---

## 7. Smoke test

```bash
python -m potts_shs.eval.run all --config potts_shs/eval/tests/mock_cluster.yaml \
  --work-dir /tmp/potts_smoke --limit 8 --input-json-dir data/datafiles/alphafold3
```

Stand-ins fabricate AF3 and DSSR output so the whole flow runs without a cluster. They reproduce
the two shapes that matter — several RNA chains per prediction, and DSSR's `m1_chain_A` keys.
`mock_dssr` derives pairs from the query sequence alone, so **every approach scores identically
under the mock**; that is the expected result there, not a finding. `--gt-cif-dir` is the one
thing the mocks cannot stand in for: annotating a deposited structure needs the real DSSR.
