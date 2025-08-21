# SyntheticEvolution
Secondary Structure Based generation of synthetic homologues that improves AlphaFold 3 Predictions

Cellular RNA function is dictated by its three-dimensional fold, yet experimental structure determination is slow and expensive.
Although protein structures can now be predicted with near-experimental fidelity, RNA 3D prediction still lags behind.
One reason is that state-of-the-art tools such as AlphaFold 3 require deep multiple sequence alignments (MSAs) that are far harder to obtain for RNA than for proteins. 
Although predicted secondary structures and large sequence foundation models can, in principle, supply covariation information, state-of-the-art predictors continue to rely almost exclusively on deep MSAs.
We remove this bottleneck by synthesising MSAs directly from secondary structure information. First, we introduce RNAformer, a transformer-based predictor that generates highly accurate secondary structures. Guided only by these predictions and evolution-inspired mutation rules, we inflate a single RNA into an arbitrarily large synthetic MSA in a fraction of a second, embedding realistic intramolecular interaction biases without relying on natural homologs.
Feeding these artificial alignments to AlphaFold 3 boosts RNA 3D accuracy across a diverse set of RNA structures, including state-of-the-art models for monomeric RNA-protein complexes. 
Secondary structure driven synthetic evolution therefore unlocks deep-alignment benefits for RNA 3D structure prediction.

**Note that this repo currently mainly serves for reproducing the results of our Nature Methods submission.**

## Install
We use [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)/[Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) to create a virtual environment for all our dependencies. If you are not familiar with using ```conda``` or ```mamba``` please find more information at <https://www.anaconda.com/docs/getting-started/miniconda/main> and <https://mamba.readthedocs.io/en/latest/index.html>. You can download miniconda for Linux, Windows, and MAC [here](https://www.anaconda.com/download) or read about Micromamba installation [here](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

Please find detailed guidelines for installation of our pipeline for the generation of synthetic homologues in our [Install Guidelines](/docs/INSTALL.md). 
For quick install, you can use our installation script by running 
```
./install.sh
```
This will install all requirements for the synthetic MSA genration, including downloading RnaBench and SPOT-RNA for structure predictions with SPOT-RNA and RNAfold.
The install script was tested on ```Ubuntu 22.04.5 LTS (Jammy Jellyfish)```.

**Note: Our full prediction pipeline requires AlphaFold 3.**
You can find more guidelines for installation in our [Install Guidelines](/docs/INSTALL.md).

## Data
You can download all required files by running
```
./download.sh
```
This will download all data files, predictions, and evaluation ```.csv``` files to reproduce the results of our Nature Methods submission.
**NOTE: The data files archive download is roughly 9GB and extracts to roughly 50GB.**

Alternatively, directly download all predictions, datafiles, and evaluations manually from
```
https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/data.tar.gz
https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/predictions.tar.gz
https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/results.tar.gz
https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/dssr.tar.gz
https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/rnaformer_models.tar.gz
```
and extract them with
```
tar -xzvf <name>.tar.gz
```
Afterwards move them to the correct directories as follows
1. ```predictions``` to ```evaluation```
2. ```dssr``` to ```evaluation```
3. ```models``` (from ```rnaformer_models.tar.gz```) to ```RNAformer```

After this, you should be ready to go.

## Usage
Please find detailed usage instructions in our [Usage Guidelines](docs/USAGE.md).
