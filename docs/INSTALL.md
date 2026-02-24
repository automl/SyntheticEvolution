# Install Guidelines
Please find the guidelines for installing required dependencies below.
Our install requires ```conda``` or ```mamba```.
For installation of these, please follow our guidelines in the repository root ```README.md```.

## Synthetic Homologous Sequences (SHS)
If you have ```conda``` or ```mamba```, we provide a ready to use install script via
```
./install.sh
```
This will install all dependencies. The script was tested on ```Ubuntu 22.04.5 LTS ((Jammy Jellyfish)```.


### Manual Install
We provide guidelines to install our environment on ```linux64``` and ```osx64```.

**Run all following commands from the project root.**

#### Linux64 Environment

For manual istall, you can run
```
conda env create -f env/environment.yml
```
Or install via our conda-lock file (requires ```conda-lock```)
```
conda-lock install --name synEvo --file env/conda-lock.yml --platform linux-64
```

Continue with the guidelines for all environments below.

#### OSX64 Environment

We provide a ```conda-lock file``` for osx. You can install the environment with
```
conda-lock install --name synEvo --file env/conda-lock.yml --platform osx-64
```
activate the environment with
```
conda activate synEvo
```

Continue with the guidelines for all environments below.

#### All Environments

Activate the envrionment with
```
conda activate synEvo
```

After that, run 
```
pip install -e .
```
to install pip dependencies.

After that, you can clone RnaBench via
```
git clone https://github.com/automl/RnaBench.git
```
and copy the lib into the correct location via

```
cd RnaBench
```
and 
```
mv RnaBench/lib .
```
For SPOT-RNA, first create the ```external_algorithms``` directory
```
mkdir -p external_algorithms
```
and move there
```
cd external_algorithms
```

You can now clone SPOT-RNA via
```
git clone https://github.com/jaswindersingh2/SPOT-RNA.git
```
For downloading SPOT-RNA models, run
```
cd SPOT-RNA && wget https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz
```
and extract the models via
```
tar -xzvf SPOT-RNA-models.tar.gz
```
Then go back to the project root.


## RNAformer
### CPU
We provide a CPU runtime patch for RNAformer so one can also run it here without any further dependencies.
However, this patch mainly serves demonstration purposes and we strongly recommend running RNAformer on GPU as desccribed below for any performance critique use cases. 

### GPU
You can use the ```infer_d``` branch of the original github repository of the RNAformer to produce secondary structure predictions using GPUs.

To do so, please follow the install instructions as desribed in the RNAformer repository: <https://github.com/automl/RNAformer>.
The models used in our manuscript can be downloaded following the instructions in the [README](/README.md)

## AlphaFold 3
We use the inference code of AlphaFold 3 to run our experiments. 

**Note that the AlphaFold 3 webserver at <https://alphafoldserver.com/> currently does not support custom MSAs.**

To install AlphaFold 3, we follow the instructions at <https://github.com/google-deepmind/alphafold3>.

**NOTE that using AlphaFold 3 requires to agree to the terms of usage and requires to request the open weights.**

## DSSR
For completeness, we also provide instructions for the usage of DSSR.
We use X3DNA-DSSR to obtain secondary structure information from 3D structures.
While not necessary for our pipeline, this step is currently required during evaluation of the 3D predictions.

However, X3DNA-DSSR is currently available under license and requires application.
You can find more infomation here: <https://x3dna.org/>.

** We provide all DSSR json files for download for reproducing our results. There is no need to install DSSR in general to generate our synthetic homologous sequences! We only use it for secondary structure evaluation!**

## US-align
We use US-align in our evaluation pipeline.
We provide links and some calls for the installs of US-align.
Please find detailed instructions on US-align at the [US-align Website](https://zhanggroup.org/US-align/).

### Source Code
To compile US-align from scratch, you can directly download the [USalign.cpp](https://zhanggroup.org/US-align/bin/module/USalign.cpp) an compile with
```
g++ -static -O3 -ffast-math -o USalign USalign.cpp
```

Alternatively, you can also directly download a zip of the binaries for your OS.
- [Linux64 Binaries](https://zhanggroup.org/US-align/bin/module/USalignLinux64.zip)
- [Win64 Binaries](https://zhanggroup.org/US-align/bin/module/USalignWin64.zip)
- [Mac Binaries](https://zhanggroup.org/US-align/bin/module/USalignMac.tar.gz)


### Help:
[USAlign Help Page](https://zhanggroup.org/US-align/help)












