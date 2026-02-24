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

### Analyze RNA families
To analyze the RNA families of the data accompanying our manuscript, run
```
python analyse_fam.py
```

### 3D plotting
To visually compare the predictions of two methods, we include a plotting scipt using [pymol](https://www.pymol.org/).
Currently, we are using the free version of pymol without the requirement to obtain a license first.
To plot all predictions, run
```
./plot_all_3d.sh
```
**Note that the script will take a while to plot all 3D structures using PYMOL.**

### Secondary structure plotting
We also provide all data for running the script for plotting secondary structures of 1SJ4 and 4P8Z by running
```
python plot_secondary_structures.py
```

You can also get predictions for any sequence with SPOT-RNA, RNAfold or the RNAformer using
```
python get_secondary_structure_prediction.py --predictor <choose your predictor> --sequence GAUGGCCGGCAUGGUCCCAGCCUCCUCGCUGGCGCCGGCUGGGCAACACCAUUGCACUCCGGUGGUGAAUGGGACU
```
choose the predictor from rnaformer, rnafold, spotrna.
The nucleotide interaction maps will be saved to the ```secondary_structure_predictions``` directory.

### Reporduce all paper figures
To reproduce the plots from the figures of our paper, run
```
python plot_results.py
```
