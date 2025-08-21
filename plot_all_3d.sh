#! /bin/bash

# plot 3D structures

echo "Plotting 3D structures..."

conda activate synEvo

echo "Start with RNAformer vs AlphaFold 3"

python plot_3d_structures3.py --rnaformer-pred-dir evaluation/predictions/rnaformer --alphafold-pred-dir evaluation/predictions/alphafold --gt-dir evaluation/predictions/gt --out-dir paper_figures --buffer 25 --alg1-name AlphaFold --alg2-name RNAformer --rf-color green --af3-color blue --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --view-store-dir pymol_views --view-mode rotation

echo "Plotting RNAformer vs SpotRNA"

python plot_3d_structures3.py --rnaformer-pred-dir evaluation/predictions/rnaformer --alphafold-pred-dir evaluation/predictions/spotrna --gt-dir evaluation/predictions/gt --out-dir paper_figures --buffer 25 --alg1-name SPOT-RNA --alg2-name RNAformer --rf-color green --af3-color tv_red --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --view-store-dir pymol_views --view-mode rotation

echo "Plotting RNAformer vs RNAfold"

python plot_3d_structures3.py --rnaformer-pred-dir evaluation/predictions/rnaformer --alphafold-pred-dir evaluation/predictions/rnafold --gt-dir evaluation/predictions/gt --out-dir paper_figures --buffer 25 --alg1-name RNAfold --alg2-name RNAformer --rf-color green --af3-color deepteal --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --view-store-dir pymol_views --view-mode rotation

echo "Plotting RNAformer vs DSSR-3D"

python plot_3d_structures3.py --rnaformer-pred-dir evaluation/predictions/rnaformer --alphafold-pred-dir evaluation/predictions/dssr --gt-dir evaluation/predictions/gt --out-dir paper_figures --views-dir pymol_views --specific-id 1SJ4 --buffer 25 --alg1-name DSSR --alg2-name RNAformer --rf-color green --af3-color purpleblue --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent

python plot_3d_structures3.py --rnaformer-pred-dir evaluation/predictions/rnaformer --alphafold-pred-dir evaluation/predictions/dssr --gt-dir evaluation/predictions/gt --out-dir paper_figures --views-dir pymol_views --specific-id 4P8Z --buffer 25 --alg1-name DSSR --alg2-name RNAformer --rf-color green --af3-color purpleblue --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent

echo "done."