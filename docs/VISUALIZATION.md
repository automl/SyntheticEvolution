# Visualization and Analysis
Here we provide guidelines how to visualize different parts of our results, including analysis of the MSA-like SHS ensembles, our analysis regarding RNA families, plotting of 2D and 3D structures, and reproduction of all figures of performance. 

## Plotting and Analysis
All plots will be saved in the ```plots``` directory.

### Analyze MSA
To analyse the MSA of all algorithms, first make sure to download all data as described in [our main README.md](/README.md).
Then, run
```
./get_all_gapped_alignments.sh
```
to generate gapped alignments from the original/modified input json files for AlphaFold 3.
You can plot individual MSA comparisons afterwards using e.g.

```
python analysis/analyse_gapped_alignment.py --fasta1 data/gapped_alignments/rnaformerN100/1SJ3_gapped_rna_alignment_rnaformer.fasta --fasta2 data/gapped_alignments/rnafold/1SJ3_gapped_rna_alignment_rnafold.fasta
```
The resulting pngs are saved to the ```plots/MSA``` directory. 


### Analyze RNA families
For completeness, we also include our analysis regarding RNA families. Simply run
```
python analyse_fam.py
```
to obtain an overview of RNAfamilies in different predictions.

### 3D plotting
To visually compare the predictions of two methods, we include a plotting scipt using [pymol](https://www.pymol.org/).
Currently, we are using the free version of pymol without the requirement to obtain a license first.
To plot the comparison of the 3D predictions of vanilla AlphaFold3 and SHS seeded by RNAformer for PDB ID 2PLY, you can run
```
python plot_3d_structures.py --rnaformer-pred-dir evaluation/predictions/rnaformerN100 --alphafold-pred-dir evaluation/predictions/alphafold3 --gt-dir evaluation/predictions/gt --out-dir paper_figures --buffer 25 --alg1-name AlphaFold --alg2-name RNAformer --rf-color green --af3-color blue --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --views-dir pymol_views --specific-id 2PLY
```
The resuting structure visuals are saved to the ```paper_figures``` directory.

For other comparisons, simply run the command above with the correct paths to other predictions, and provide a valid id as well as colors of choice.
For example, to reproduce the predictions of the different SHS seeded variants for PDB id 1SJ4, the following commands were used:

#### SPOT-RNA
```
python plot_3d_structures.py --rnaformer-pred-dir evaluation/predictions/rnaformerN100 --alphafold-pred-dir evaluation/predictions/spotrna --gt-dir evaluation/predictions/gt --out-dir paper_figures --buffer 25 --alg1-name SPOT-RNAN100 --alg2-name RNAformerN100 --rf-color green --af3-color tv_red --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --views-dir pymol_views --specific-id 1SJ4
```
#### RNAfold
```
python plot_3d_structures.py --rnaformer-pred-dir evaluation/predictions/rnaformerN100 --alphafold-pred-dir evaluation/predictions/rnafold --gt-dir evaluation/predictions/gt --out-dir paper_figures --buffer 25 --alg1-name RNAfoldN100 --alg2-name RNAformerN100 --rf-color green --af3-color deepteal --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --views-dir pymol_views --specific-id 1SJ4
```
#### Ground Truth (DSSR)
```
python plot_3d_structures.py --rnaformer-pred-dir evaluation/predictions/rnaformerN100 --alphafold-pred-dir evaluation/dssr --gt-dir evaluation/predictions/gt --out-dir paper_figures --views-dir pymol_views --buffer 25 --alg1-name DSSRN100 --alg2-name RNAformerN100 --rf-color green --af3-color purpleblue --gt-color gray --gt-rna-color brightorange --pred-rna-color magenta --plot-rna-only --show-plot --transparent --views-dir pymol_views --specific-id 1SJ4
```

### Secondary structure plotting
**Note: Our secondary structure prediction script depends on RnaBench visualization with additional requirements. Before running the following commands, please additionally install the following dependencies using pip after activating the synEvo environment:**

```
pip install plotly matplotlib_venn upsetplot
```

You can also get predictions for any sequence with SPOT-RNA, RNAfold or the RNAformer using
```
python get_secondary_structure_prediction.py --predictor <choose your predictor> --sequence GAUGGCCGGCAUGGUCCCAGCCUCCUCGCUGGCGCCGGCUGGGCAACACCAUUGCACUCCGGUGGUGAAUGGGACU
```
choose the predictor from rnaformer, rnafold, spotrna.
The nucleotide interaction maps will be saved to the ```secondary_structure_predictions``` directory.

We also provide the script for plotting secondary structures used in the paper
```
python plot_secondary_structures.py
```
**However, we note that this script depends on analysis obtained by running DSSR. Please see our [install instructions](/docs/INSTALL.md) for more information on DSSR.**

### Reporduce all plots from paper figures
To reproduce the plots from the figures of our paper, run
```
python plot_results.py
```
You can find all plots in the ```plots/performance``` directory afterwards.
