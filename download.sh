#! /bin/bash

echo "download data files."

wget https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/data.tar.gz
tar -xzvf data.tar.gz
rm data.tar.gz

echo "Download predictions"

wget https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/predictions.tar.gz
tar -xzvf predictions.tar.gz
rm predictions.tar.gz
mv predictions/ evaluation/

echo "Download dssr json files"

wget https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/dssr.tar.gz
tar -xzvf dssr.tar.gz
rm dssr.tar.gz
mv dssr/ evaluation/

echo "Download results"

wget https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/results.tar.gz
tar -xzvf results.tar.gz
rm results.tar.gz

echo "Download rnaformer models"

wget https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/rnaformer_models.tar.gz
tar -xzvf rnaformer_models.tar.gz
rm rnaformer_models.tar.gz
mv models/ RNAformer

echo "Done with downloads."


