#! /bin/bash

echo "Evaluate RNAformer N100 RNA only"

./run_evaluation.sh evaluation/predictions/rnaformerN100/ evaluation/predictions/gt/ rna_rna

cd results

ls -lah

mv pred_rna_rna.csv rnaformerN100_rna_rna.csv
mv exp_rna_rna.csv pdb_rna_rna_from_rnaformerN100_eval.csv

ls -lah

cd ..

echo "Evaluate RNAformer N5000 RNA only"

./run_evaluation.sh evaluation/predictions/rnaformerN5000/ evaluation/predictions/gt/ rna_rna

cd results

ls -lah

mv pred_rna_rna.csv rnaformerN5000_rna_rna.csv
mv exp_rna_rna.csv pdb_rna_rna_from_rnaformerN5000_eval.csv

ls -lah

cd ..

echo "Evaluate RNAformer N20000 RNA only"

./run_evaluation.sh evaluation/predictions/rnaformerN20000/ evaluation/predictions/gt/ rna_rna

cd results

ls -lah

mv pred_rna_rna.csv rnaformerN20000_rna_rna.csv
mv exp_rna_rna.csv pdb_rna_rna_from_rnaformerN20000_eval.csv

ls -lah

cd ..

echo "Evaluate RNAformer N100 Protein-RNA"

./run_evaluation.sh evaluation/predictions/rnaformerN100/ evaluation/predictions/gt/ protein_rna

cd results

ls -lah 

mv pred_protein_rna.csv rnaformerN100_protein_rna.csv
mv exp_protein_rna.csv pdb_protein_rna_from_rnaformerN100_eval.csv

ls -lah

cd ..

echo "Evaluate RNAformer N5000 Protein-RNA"

./run_evaluation.sh evaluation/predictions/rnaformerN5000/ evaluation/predictions/gt/ protein_rna

cd results

ls -lah

mv pred_protein_rna.csv rnaformerN5000_protein_rna.csv
mv exp_protein_rna.csv pdb_protein_rna_from_rnaformerN5000_eval.csv

ls -lah

cd ..

echo "Evaluate RNAformer N20000 Protein-RNA"

./run_evaluation.sh evaluation/predictions/rnaformerN20000/ evaluation/predictions/gt/ protein_rna

cd results

ls -lah

mv pred_protein_rna.csv rnaformerN20000_protein_rna.csv
mv exp_protein_rna.csv pdb_protein_rna_from_rnaformerN20000_eval.csv

ls -lah

cd ..

echo "Done."