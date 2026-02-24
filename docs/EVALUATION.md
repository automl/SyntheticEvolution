# Evaluation
We provide all csvs from our 3D evaluations in the ```results/csv``` directory.
However, for reproducibility, we also provide scripts for rerunning our evaluations.

**Note: Evaluations of 3D predictions might take a bit, especially when rerunning our scripts for the full benchmark.**

You can run the evaluation of 3D predictions with 
```
./run_evaluation.sh <prediction_data_dir> <true_data_dir> <sample type; either rna_rna (includes RNA monomers) or protein_rna>
```
The ouput are two csvs saved to the ```results``` directory, ```pred_protein_rna.csv``` (or ```pred_rna_rna.csv```) with the computed metrics for the prediction and ```exp_protein_rna.csv``` (or ```exp_rna_rna.csv```) for further information on each sample.

Predictions are provided in the ```evaluation/predictions``` directory for all methods evaluated.

We also provide an examle script for the evaluation of all RNAformer seeded SHS predictions in ```evaluate_all.sh```. Simply run
```
./evaluate_all.sh
```
from the project root to redo all evaluations of RNAformer seeded SHS predictions.
