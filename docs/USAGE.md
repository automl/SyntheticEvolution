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
Per default, the generation pipeline will use the parameters that were also used in our manuscript.
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

For structure prediction, we currently natively support ```rnafold```, ```rnaformer```, and ```spotrna```. However, more predictors are generally available via RnaBench and can easily be included into our SHS generation script e.g. to serve benchmarking purposes.
Please see RnaBench at ```https://github.com/automl/RnaBench``` for more information and guidelines on the installation and usage of several other folding engines.


All figures will be saved to the ```paper_figures``` directory.
