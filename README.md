# Reproducing the NeurIPS 2026 SHS paper (camera-ready)

Code and small data to regenerate every AlphaFold 3 table and the rebuttal analyses of
*"Encoding RNA Topology into Synthetic Alignments for 3D Structure Prediction"*.
Every AF3 condition is evaluated with **4 AF3 seeds (1, 0, 42, 137) × 5 diffusion samples = 20 predictions per target**;
a target's score is the mean over its 20 predictions.

Large files are separate downloads from
<https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/>
(`bash download_data.sh` fetches, verifies and unpacks them into `downloads/`):

| Bundle | Size | Contents |
|---|---|---|
| [`2026_neurips_shs_inputs.zip`](https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/2026_neurips_shs_inputs.zip) | 26.8 GB | AF3 input JSONs (natural MSA / SHS / controls embedded), reference mmCIFs |
| [`2026_neurips_shs_outputs.zip`](https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/2026_neurips_shs_outputs.zip) | 20.1 GB | AF3 structures (`model.cif`) + summary confidences for all seeds/samples, RhoFold+ outputs, metric DBs (`metric_dbs/`); full PAE/pLDDT matrices omitted |

sha256: inputs `2131fa721fa0c38dae15a957d94d08386b0eef9fd62fec381a1310ecb1463aca`,
outputs `47ca85aacba9814b532d73675dab0ddff130d5dd076e45c175835176a98e276a` (also published as `<zip>.sha256`).
Unpacked: ~218 GB, ~362k files.

## Layout

```
download_data.sh     fetch + verify + unpack the two bundles into downloads/
config.sh            paths + env names; defaults work out of the box (override in env or config.local.sh)
data/                target ids, paper splits (neurips_ids.py), RNAformer 2D predictions, analysis inputs
shs/                 Synthetic Co-evolution Engine (SHS generator) + the recipes used for every input set
fold/                AF3 inference driver, input selection/re-seeding, the CRC fold submissions
eval/                scoring engine (US-align + OpenStructure), shard driver, condition registry, RhoFold+ helper
tables/              table generator, metric export, finalize script
analysis/            rebuttal / CRC analyses (paired tests, gates, ladders, construction sweep)
results/             tables/*.tex|csv, metrics_per_prediction.csv.gz, analysis outputs (DBs: download)
provenance/          versions, checks, origins, sha256 manifests, submission-era table code
env/                 conda specs: eval_analysis (synEvo), ost, shs_generator (RnaBench)
```

## Quick path: every number from the downloads

```bash
conda env create -f env/eval_analysis.yml          # python 3.9 (env name synEvo); RnaBench env optional (cd-hit-est-2d)
bash download_data.sh                              # ~47 GB download -> downloads/ (resumable, checksums verified)
bash tables/finalize.sh                            # tables, metrics export, all analyses, checksums (~1 h, 16 cores)
```

`finalize.sh` builds US-align if needed and writes `results/tables/`, `results/metrics_per_prediction.csv.gz` and
`results/analysis/` (logs in `results/analysis/logs/`). Metric DBs: your own re-scored `results/db*` if present (step 4 below),
otherwise the shipped `metric_dbs/` of the outputs bundle.
No path needs editing: everything is resolved relative to the repository and `downloads/`
(`CRC_DATA`, `CRC_INPUTS`, `CRC_OUTPUTS`, `GT_*`, ... can be overridden in the environment or in an untracked `config.local.sh`).

## Full path: predictions → metrics → tables

1. **Environments**: `env/*.yml`; build US-align with `bash eval/build_usalign.sh`; AF3 paths via `AF3_HOME` etc. (`config.sh`).
   AlphaFold 3 v3.0.0 (weights from Google DeepMind) is required only to re-fold.
2. **Inputs**: use the inputs bundle, or regenerate SHS with `shs/generate_shs.sh` (RnaBench env; recipes in `shs/recipes/`).
   Seed variants: `fold/prepare_inputs.py --seeds 1`.
3. **Fold**: `fold/fold_array.sh <input_root> <output_dir> <list>` (`--run_data_pipeline=False`);
   the runs added for the CRC are in `fold/runs/submit_crc_folds.sh`.
4. **Score**: `eval/make_shards.py` writes one shard per (condition, seed, sample) from `eval/conditions.tsv`
   (pointing into the unzipped outputs bundle; the Slurm scripts use cluster partition names, adapt the `#SBATCH` lines);
   `sbatch --array=0-<n-1> --export=CRC_ROOT=$PWD eval/run_eval_array.sh <shards.tsv>` → `results/db/<tag>_<idset>_s<seed>_m<sample>.db`
   (re-scoring `rnaformer_single_s0_m0` reproduced the shipped DB exactly; see `provenance/PROVENANCE.md`).
   RhoFold+ outputs: `eval/rhofold/rhofold_to_eval_layout.py`, then the same driver.
5. **Tables + analyses**: `bash tables/finalize.sh`.

## Results map

| Paper / rebuttal result | Script |
|---|---|
| All AF3 tables (single chains, RNA-protein complexes, splits, 1v1/multichain, controls, ordering, denoise, SHS+MSA, depth scaling) | `tables/make_crc_tables.py` |
| SHS depth scaling N=20 … 20000 (single chains, complexes) | `tables/make_crc_tables.py` → `rna_single_chain_scaling`, `rna_protein_complexes_scaling` |
| Seed robustness (submitted seed 1 vs 4×5) | `tables/make_crc_tables.py` → `results/tables/seed_robustness.csv` |
| RhoFold+ tables (incl. denoise ladder and native-MSA comparison) | `tables/make_crc_tables.py` from `metric_dbs/rhofold` (re-score: `eval/rhofold/score_rhofold.sh`) |
| Orphan parity (paired Wilcoxon + bootstrap CI) | `analysis/rebuttal/multiseed_orphan_parity.py` |
| Near-floor solvable orphans | `analysis/rebuttal/uTY4_q4_near_floor.py` |
| 2D-accuracy gates, self-consistency, complex bins, ribozymes, family leakage | `analysis/rebuttal/multiseed_nomsa_final.py`, `multiseed_nomsa_cis.py` |
| Denoising mechanism (AF3) | `analysis/rebuttal/denoise_ms.py` |
| AF3 confidence (pTM) ladder | `analysis/rebuttal/confidence_ladder.py`, `confidence_ladder_confirm.py` |
| Long RNAs (>200 nt) ladder; homolog-search runtime | `analysis/rebuttal/general2_long_run_ladder.py`, `general2_runtime.py` |
| Sequence/structure novelty of post-cutoff targets | `analysis/rebuttal/faN3_q1_similarity.py` (complex structure scan: `finalize.sh --full`) |
| 4P8Z/6WJS interventions, 91-RNA construction sweep (N=1000 and N=100) | `analysis/construction/eval_intervention_ms.py`, `construct91_eval.py` → `construct91_within.py`, `construct91_plot.py` (`CONSTRUCT_N=100`) |

## Notes

* Target sets: 91 single-chain RNAs, 377 RNA–protein complexes, 25 long RNAs (`data/ids/`); splits from `data/splits/neurips_ids.py`.
* SHS from the reference topology (`dssr`) covers 89/91 single chains (7EFG, 7VFT not folded).
* AlphaFold 3 outputs in the outputs bundle are subject to the AlphaFold 3 Output Terms of Use.
* Details, checks and origins: `provenance/PROVENANCE.md`.
