# SyntheticEvolution: NeurIPS 2026 reproduction

This repository contains the code for the following papers:

* J. K. H. Franke\*, I. K. Tennie\*, M. Uhl, R. Köksal, D. Matus, D. Scheuer, F. Hutter, R. Backofen, F. Runge.
  *Encoding RNA Topology into Synthetic Alignments for 3D Structure Prediction.* NeurIPS 2026.
  This branch reproduces the results of this paper.
* J. K. H. Franke\*, I. K. Tennie\*, M. Uhl, R. Köksal, D. Matus, F. Hutter, R. Backofen, F. Runge.
  *Synthetic RNA Evolution Enables Accurate Alignment-Free 3D Structure Prediction.* GEM Workshop at ICLR 2026.
  [OpenReview](https://openreview.net/forum?id=wWKuI35OW7)

\* Equal contribution.

The 3D predictions are scored with an early version of the evaluation pipeline described in
I. K. Tennie, M. Uhl, J. K. H. Franke, R. Backofen, F. Runge. *AlphaFold 3 Prediction Accuracy and Reliability across
RNA Structures and Protein-RNA Interfaces.* bioRxiv 2026, [doi:10.64898/2026.11.02.763690](https://doi.org/10.64898/2026.11.02.763690).
The full version of that pipeline is at [automl/RNA-Protein-3D-Structure-Evaluation](https://github.com/automl/RNA-Protein-3D-Structure-Evaluation).
RNAformer, the 2D predictor used to build the synthetic alignments, is at [automl/RNAformer](https://github.com/automl/RNAformer).

All AlphaFold 3 results use 4 seeds with 5 diffusion samples each, so 20 predictions per target.
A target's score is the mean over these 20 predictions.

## Data

The predictions and inputs are too large for git. They are available from
<https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/>, and
`bash download_data.sh` downloads, checks and unpacks them into `downloads/`.

| File | Size | Contents |
|---|---|---|
| [`2026_neurips_shs_inputs.zip`](https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/2026_neurips_shs_inputs.zip) | 26.8 GB | AlphaFold 3 input files with natural MSAs, SHS and control alignments; reference structures |
| [`2026_neurips_shs_outputs.zip`](https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/2026_neurips_shs_outputs.zip) | 20.1 GB | All AlphaFold 3 and RhoFold+ predictions with their confidence files, and the metric databases |
| [`2026_neurips_shs_extras.zip`](https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc/2026_neurips_shs_extras.zip) | 1.7 MB | The example structures shown in the figures (1SJ4, 4P8Z) that are not in the other two files |

SHA-256 checksums: inputs `2131fa721fa0c38dae15a957d94d08386b0eef9fd62fec381a1310ecb1463aca`,
outputs `47ca85aacba9814b532d73675dab0ddff130d5dd076e45c175835176a98e276a`,
extras `848283cc8741268598eeb69a41317a800757587d29198a7d491032738c532c5e`.
Each zip also has a `.sha256` file next to it. Unpacked, the data takes about 218 GB.

## Installation

```bash
conda env create -f env/eval_analysis.yml   # analysis and scoring (env name synEvo)
conda env create -f env/ost.yml             # OpenStructure, only needed to re-score predictions
conda env create -f env/shs_generator.yml   # only needed to generate new SHS
bash eval/build_usalign.sh                  # US-align (finalize.sh builds it too if it is missing)
```

AlphaFold 3 (v3.0.0) is only needed if you want to fold the inputs again. Set its paths in `config.sh`.

## Reproducing the tables from the downloaded data

```bash
bash download_data.sh      # about 47 GB
bash tables/finalize.sh    # about 1 hour on 16 cores
```

This writes all tables to `results/tables/`, all analyses to `results/analysis/`, and every single prediction score to
`results/metrics_per_prediction.csv.gz`. You do not need to edit any paths. If you want to point to data somewhere
else, set the variables from `config.sh` in your environment or in a `config.local.sh` file.

## Running everything from scratch

1. Take the input files from the inputs zip, or generate new SHS with `shs/generate_shs.sh` (settings in `shs/recipes/`).
2. Fold them with AlphaFold 3 using `fold/fold_array.sh`. The scripts we used are in `fold/runs/`.
3. Score the predictions. `eval/make_shards.py` splits the work into one job per condition, seed and sample, and
   `eval/run_eval_array.sh` runs these jobs as a Slurm array. The results go to `results/db/`.
   For RhoFold+, first convert the outputs with `eval/rhofold/rhofold_to_eval_layout.py`.
4. Compute the complex RMSD with `eval/crmsd/run_crmsd_array.sh` (see `eval/crmsd/README.md`).
5. Self-consistency needs DSSR, which we cannot distribute. With your own DSSR binary, run `eval/selfcons/selfcons.py`.
6. Run `bash tables/finalize.sh` as above.

Steps 4 and 5 are optional: their results are already in `results/crmsd/` and `results/selfcons/`.
The Slurm scripts use the partition names of our cluster, so you will need to change the `#SBATCH` lines.

## Where each result comes from

| Result | Script |
|---|---|
| Main AlphaFold 3 and RhoFold+ tables (all splits, controls, denoising, SHS depth, seed robustness) | `tables/make_crc_tables.py` |
| Remaining tables (orphan benchmark, long RNAs, construction steps, ribozymes, and others) | `tables/make_extra_tables.py` |
| Statistics quoted in the text (p-values, confidence intervals) | `tables/prose_numbers.py` |
| Orphan comparison and near-floor orphans | `analysis/benchmark/orphan_parity.py`, `orphan_near_floor.py` |
| Fig. 3a (3D accuracy against 2D accuracy) | `analysis/benchmark/fig3a_mcc_vs_3d.py` |
| Fig. 3b (arc diagrams for 1SJ4) | `analysis/benchmark/fig3b_arcs.py --pdb 1SJ4` (needs VARNA) |
| 2D accuracy groups, self-consistency, ribozymes, family overlap | `analysis/benchmark/twod_gates.py`, `twod_gates_ci.py` |
| Agreement of predictions with the encoded topology | `analysis/benchmark/encoded_agreement.py` |
| Denoising | `analysis/benchmark/denoise.py` |
| AlphaFold 3 confidence | `analysis/benchmark/confidence_ladder.py` |
| Long RNAs and homology search time | `analysis/benchmark/long_rna_ladder.py`, `long_rna_runtime.py` |
| Novelty of post-cutoff targets | `analysis/benchmark/postcutoff_novelty.py` (the structure scan runs with `finalize.sh --full`) |
| Differences between diffusion samples | `analysis/benchmark/across_samples.py` |
| 4P8Z and 6WJS case studies, construction sweep over 91 RNAs | `analysis/construction/` |
| 2D accuracy of the construction-sweep predictions | `eval/dssr2d/construction_2d.py` |

## Notes

* We evaluate 91 single-chain RNAs, 377 RNA-protein complexes and 25 long RNAs. The ID lists are in `data/ids/` and the
  splits in `data/splits/neurips_ids.py`.
* SHS from the reference topology exist for 89 of the 91 single chains (not for 7EFG and 7VFT).
* The complex RMSD is not the `Complex_RMSD` column of the metric databases. See `eval/crmsd/README.md`.
* The AlphaFold 3 predictions in the outputs zip fall under the AlphaFold 3 Output Terms of Use.
* Software versions and our consistency checks are listed in `provenance/PROVENANCE.md`.
