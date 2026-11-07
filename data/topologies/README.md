# SPOT-RNA and RNAfold topologies

The SHS generator runs SPOT-RNA and RNAfold itself when it builds an alignment. This folder keeps a copy of those
predictions so that analyses such as Fig. 3a (`analysis/benchmark/fig3a_mcc_vs_3d.py`) can use them without rerunning
the predictors.

There is one JSON file per predictor and target set: `{rnafold,spotrna}_{single,complex1v1}.json`. The sets are the 91
single-chain RNAs and the 112 complexes with one RNA and one protein chain. Each entry has a `pdb_id`, the `sequence`
and the predicted `pairs` (0-based).

The files were made with `analysis/topologies/predict_2d.py` (ViennaRNA 2.6.3 for RNAfold, the 5-model SPOT-RNA
ensemble on CPU via `run_spotrna_array.sh`). We checked them against the shipped SHS: for all four files, the paired
positions are exactly the conserved columns of the alignments.
