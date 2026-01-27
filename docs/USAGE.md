# Usage

If you haven't already done so, please start with reading the [Install Guidelines](/docs/INSTALL.md).

All following instructions expect that all files were already downloaded as described in [our main README.md](/README.md).
All calls should be run from the project root.

## Generation of Synthetic Homologous Sequences

You can run our synthetic pipeline on any monomeric RNA-protein complex by providing an RNA sequence, a protein sequence, and a RNA secondary structure in dot-bracket format.
For example run
```
python rna_msa_generator_base_pair.py --structure "..((((....))))...." --rna-seq AGCGCGUAACGAUAGCUA --protein-seq MKTIIALSYIFCLVFAGQDEIRTLVSRVELTKLSDKIAARHGLQEVNRAALGRGIVRVA --seed 42 --show_plot
```
Per default, the generation pipeline will use the parameters that were also used in our Nature Biotechnology submission.
The ```--show_plot``` option enables plotting of some MSA features.
The output is a json file that can be used as input for AlphaFold 3.

Alternatively, you can run the msa generation to modify an existing json data file of a previous AlphaFold 3 prediction.

For example, simply run
```
python rna_msa_generator_base_pair.py --structure_predictor rnafold --input_json_path data/datafiles/alphafold3/1SJ4_data.json
```
This will generate synthetic homologues using the RNAfold.
Alternatively, you can use rnaformer or spotrna as predictors.
The output is stored in the ```custom_msa_json_output``` directory. 

You can always provide a pdb_id via ```--pdb_id``` to name your output files file.

For example, you can use
```
python rna_msa_generator_base_pair.py --structure_predictor rnaformer --input_json_path data/datafiles/alphafold3/1SJ4_data.json --pdb_id myverylongcustomid
```

**Note that we only provide a runtime patch for RNAformer that allows to run it on CPU. The runtime of RNAformer is substantially increased compared to running on GPU and we observed some differences in the predictions when using the CPU variant provided here.**

## Evaluation
We also provide scripts for rerunning our evaluation. 

You can run the eval with 
```
./run_evaluation.sh <prediction_data_dir> <true_data_dir> <sample type; either rna_rna (includes RNA monomers) or protein_rna>
```
The ouput are two csvs saved to the ```results``` directory, ```pred_protein_rna.csv``` (or ```pred_rna_rna.csv```) with the computed metrics for the prediction and ```exp_protein_rna.csv``` (or ```exp_rna_rna.csv```) for further information on each sample.

## Plotting and Analysis
All plots will be saved in the ```plots``` directory.

### Analyze MSA
To analyse the MSA of different algorithms, first make sure to download all data as described above.
Then, run
```
./get_all_gapped_alignments.sh
```
to generate gapped alignments from the original/modified input json files for AlphaFold 3.
You can plot individual MSA comparisons afterwards using e.g.

```
python analysis/analyse_gapped_alignment.py --fasta1 data/gapped_alignments/rnaformer/1SJ3_gapped_rna_alignment_rnaformer.fasta --fasta2 data/gapped_alignments/rnafold/1SJ3_gapped_rna_alignment_rnafold.fasta
```

### Analyze RNA families
To analyze the RNA families of the data accompanying our Nature Biotechnology submission, run
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
All figures will be saved to the ```paper_figures``` directory.
