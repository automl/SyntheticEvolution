#! /bin/bash

python get_gapped_alignment.py --algorithm alphafold --data_dir data/datafiles/alphafold3 --out_dir data/gapped_alignments/alphafold3
python get_gapped_alignment.py --algorithm rnaformer --data_dir data/datafiles/rnaformer --out_dir data/gapped_alignments/rnaformer
python get_gapped_alignment.py --algorithm rnafold --data_dir data/datafiles/rnafold --out_dir data/gapped_alignments/rnafold
python get_gapped_alignment.py --algorithm spotrna --data_dir data/datafiles/spotrna --out_dir data/gapped_alignments/spotrna
python get_gapped_alignment.py --algorithm dssr --data_dir data/datafiles/dssr --out_dir data/gapped_alignments/dssr