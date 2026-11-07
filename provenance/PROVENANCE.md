# Provenance

This file lists the software versions behind the results and the checks we ran to make sure the shipped data and code
give the numbers in the paper.

## Software

| Component | Version |
|---|---|
| AlphaFold 3 | v3.0.0 (`google-deepmind/alphafold3` at `2e562355`), container `alphafold3.sif` sha256 `a32b4f5d…`, weights `af3.bin` sha256 `f3c2f9c1…`. Run with `--run_data_pipeline=False`, because every input already contains its MSA and templates. |
| AlphaFold 3 sampling | Default settings, 4 seeds with 5 diffusion samples each, on A100 80GB GPUs |
| Evaluation engine | `eval/engine/` is a copy of AF3InterfaceEval at `a23a04f` (see `UPSTREAM_COMMIT`). The only change is one commented-out example path in `database/exportDB2CSV.py`. |
| US-align | Version 20241108 (`USalign.cpp` is included), built with `c++ -O3 -ffast-math`, binary sha256 `d0b3d601…` |
| OpenStructure | 2.11.1 (`env/ost.yml`) |
| SHS generator | `shs/rna_msa_generator_base_pair.py`, settings config80. The construction sweep uses its `construct` branch. |
| RNAformer 2D predictions | `data/rnaformer/rnaformer_all_preds.pkl` |

The checksums of all inputs, metric databases and reference structures are in `provenance/manifests/`
(rebuild them with `provenance/make_manifests.sh`).

## Checks

* Running `shs/generate_shs.sh` with the included generator gives exactly the shipped N=5000 complex input files
  (apart from the seed field).
* Scoring a shard again with the included engine gives exactly the shipped metric
  database, in every row and column.
* The inputs for the SHS depth sweep (N = 20 to 20000) were made from the folded config80 inputs. We checked every file:
  the RNA MSA has depth N, and everything else is identical to the N=100 inputs.
* Building the 91-RNA construction sweep at N=1000 with the N=100 code gives the original N=1000 inputs byte for byte.
* `tables/finalize.sh`, run from the downloaded data in a clean environment (`env -i`), reproduces all result files.

## Complex RMSD

The complex RMSD is computed with PyMOL `align('predictedA and polymer', 'referenceA and polymer')` over all protein and
RNA atoms. We report the RMSD before outlier removal (`align(...)[3]`). The metric databases instead store the value
after PyMOL's 5 default rejection cycles (`align(...)[0]`) in `Complex_RMSD[0]`.

`eval/crmsd/crmsd.py` computes both values for all 90,480 complex predictions (240 shards). Tables and
`metrics_per_prediction.csv.gz` use the value before outlier removal from `results/crmsd/`. As a check, the value after
outlier removal matches `Complex_RMSD[0]` in the databases for all 90,480 predictions (largest difference 5e-7).

## Self-consistency

For each prediction, `eval/selfcons/selfcons.py` computes the MCC between the RNAformer topology that was encoded in the
SHS and the base pairs that DSSR v2.2.1 finds in the AlphaFold 3 structure (full L x L matrices over the query
sequence). The topology is looked up in `rnaformer_all_preds.pkl` in the same way as in the SHS generator.
`results/selfcons/` has these values for all 20 predictions of the 91 single chains, for RNAformer SHS, no MSA, shuffled
SHS, random SHS and the natural MSA. For RNAformer SHS, this is the self-consistency score. DSSR is not distributed.

## 2D accuracy of the encoded topologies

`data/gates/rnaformer_2d_vs_gt.csv` and the `rf2d_pred`, `pred_bin` and `has2d` columns of
`data/gates/master_combined.csv` hold the MCC between the RNAformer topology encoded in the SHS and the DSSR pairs of the
reference structure. The reference pairs come from `data/construct/construct_pairs_91.json` (single chains) and from the
`pairs2` column of `data/gates/inputs/rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv` (complexes with one RNA and
one protein chain). The script is `analysis/gates/build_rnaformer_2d_gt.py --update-master`. These are the same values as
the RNAformer x-axis of Fig. 3a. The SPOT-RNA and RNAfold topologies are in `data/topologies/`.

We checked that the paired positions of every topology are the conserved columns of the matching shipped SHS. The
reference pairs are the topology of the reference-topology SHS for single chains, and match a DSSR v2.2.1 run on the
reference structures for complexes.

## 2D accuracy in the construction sweep

`results/dssr2d/construction_2d.tsv.gz` compares the DSSR v2.2.1 base pairs of every prediction in the 4P8Z and 6WJS
construction sweeps with the DSSR pairs of the reference structures (full-matrix MCC). That is 6 topology stages x 14
SHS depths x 4 seeds x 5 samples = 3360 predictions.

## Other notes

* SHS from the reference topology exist for 89 of the 91 single chains (not for 7EFG and 7VFT).
* RhoFold+ has no seeds, so every input is run once. We ran SHS from RNAfold, SPOT-RNA and RNAformer, the RNAfold-denoised
  SHS, SHS from the reference topology, and RhoFold+ with its own MSA pipeline. The inputs are in `rhofold_inputs/` (the
  A3M files are the same config80 SHS as for AlphaFold 3), and they are scored with the same engine
  (`eval/rhofold/score_rhofold.sh`, results in `metric_dbs/rhofold`).
* SHS depth sweep: N = 20, 100, 500, 1000, 5000, 10000, 15000 and 20000 (RNAformer SHS, config80, 4 seeds x 5 samples,
  91 single chains and 377 complexes). The fold jobs are in `fold/runs/`.
* The 91-RNA construction sweep was run at N=1000 and N=100 (`CONSTRUCT_N=100`), with the same pairs, settings and seeds.
* AlphaFold 3 weights, container and databases cannot be redistributed.
