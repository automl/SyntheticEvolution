# Provenance

## What was run, with what

| Component | Version / identity |
|---|---|
| AlphaFold 3 | source `google-deepmind/alphafold3` @ `2e562355` (v3.0.0); container `alphafold3.sif` sha256 `a32b4f5d…`; weights `af3.bin` sha256 `f3c2f9c1…`; local `run_alphafold.py` has the `max_template_date` lines commented out (irrelevant here: all inputs carry MSAs/templates, inference ran with `--run_data_pipeline=False`) |
| AF3 sampling | defaults (5 diffusion samples); seeds from the input JSON `modelSeeds`: seed 1 (submission) + 0, 42, 137 (rebuttal/CRC); A100 80GB |
| Eval engine | `eval/engine/` = files of AF3InterfaceEval @ `a23a04f` (`UPSTREAM_COMMIT`), unmodified except one commented-out example path in `database/exportDB2CSV.py` (personal path replaced by `<path>`; no code change); basic tier only |
| US-align | Version 20241108 (`USalign.cpp` vendored); Linux build `c++ -O3 -ffast-math`, binary sha256 `d0b3d601…` |
| OpenStructure | 2.11.1 (`env/ost.yml`) |
| SHS generator | `shs/rna_msa_generator_base_pair.py` (submission version + additive `construct` branch); config80, seed 42 |
| RNAformer 2D | `data/rnaformer/rnaformer_all_preds.pkl` (predictor `rnaformer`) |

## Checks performed when assembling this package (2026-10-04)

* `shs/generate_shs.sh` parameters + the vendored generator reproduce an existing N=5000 complex input JSON exactly (only `modelSeeds` differs).
* The vendored engine re-scores `rnaformer_single_s0_m0` with all 91 rows / all columns identical to the rebuttal DB.
* `tables/make_crc_tables.py --submitted` (submitted seed-1 evaluation) reproduces all 507 table cell means of the submitted paper that it generates (|Δ| ≤ 0.0015).
* The copied rebuttal scripts reproduce the posted rebuttal numbers on seeds 0/42/137 (orphan parity +0.019 [−0.024, +0.061] / +0.004 [−0.011, +0.018]; denoise +0.006 p=0.23, reference +0.091).

## Where things came from (cluster paths, originals untouched)

| Here | Original |
|---|---|
| `results/db/*_s{0,42,137}_*.db` | `reproduce_rebuttal/results/benchmark_multiseed/` (rebuttal, July 2026) |
| `results/db_long/` | `reproduce_rebuttal/results/long_run/` |
| `data/published_seed1/` | iris_eval CSVs of the submission (`reproduce_rebuttal/results/published_seed1/`) |
| `data/splits/neurips_ids.py`, `provenance/submission/` | `~/dssr_parsing/` (submission table pipeline) |
| `data/gates/` | `tmp/rebuttal_bundle/reviewer_analysis/` (2D gates, family labels; `metadata_468.csv` RNA families were merged with an internal DB) |
| `analysis/rebuttal/` | `reproduce_rebuttal/scripts/analysis/` + session scratch (denoise_ms, orphan_*, struct_vs_delta) |
| `analysis/construction/` | session scratch (construct91, 4P8Z/6WJS multiseed) |
| `shs/recipes/` | `~/RnaBench/create_*.sh` (how every SHS input set was generated) |
| `eval/shards/rebuttal_original/` | shard lists of the rebuttal eval |

## Known gaps

* `dssr` (SHS from the reference topology) single chains: 89/91 — 7EFG and 7VFT were never folded.
* RhoFold+ (one deterministic run per input; no model seeds): all six runs of the paper are in the outputs bundle (`rhofold/`:
  SHS from RNAfold / SPOT-RNA / RNAformer, RNAfold-denoised, reference topology, native RhoFold+ MSA) with their inputs
  (`rhofold_inputs/`, A3M = exactly the AF3 config80 SHS). They are re-scored with the same engine
  (`eval/rhofold/score_rhofold.sh` -> `metric_dbs/rhofold`) and the tables use these DBs. Versus the submitted evaluation:
  TM and RMSD identical for every target; all-atom LDDT +0.009 to +0.014 on average (the same engine-version effect seen for AF3,
  +0.003 to +0.005). A later SHS-generator-seed sweep of the RhoFold+ runs (seeds 0/42/137) used generator defaults, not
  config80, and is deliberately not included.
* SHS depth scaling is complete for N = 20, 100, 500, 1000, 5000, 10000, 15000, 20000 (RNAformer SHS, config80, 4 seeds x 5
  samples, 91 single chains + 377 complexes); the CRC additions are in `fold/runs/submit_crc_depth_folds.sh` (inputs re-seeded
  from the folded config80 JSONs; every input JSON checked: RNA MSA depth = N, all other content identical to the N=100 inputs).
* 91-RNA construction sweep: N=1000 (original) and N=100 (CRC; `CONSTRUCT_N=100`), same pairs/config80/seeds; the N=100 build
  reproduced the N=1000 inputs byte-for-byte when run at N=1000.
* AlphaFold 3 weights, container, and databases are not redistributable.
