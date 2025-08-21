<!---
# Source Code for the RNAformer ICLR 2025 Submission

### Model Checkpoints

The model checkpoints are stored in the 'models' directory and can be unpacked using the following command:

```bash
find models -name "*.tar.xz" -exec tar -xvf {} -C models --use-compress-program=xz \;
```

### Datasets

The datasets are stored in the 'datasets' directory and can be unpacked using the following command:

```bash
find datasets -name "*.tar.xz.aa" -print0 | xargs -0 -I {} bash -c 'basename=$(basename "$1" .tar.xz.aa); cat "${1%aa}"* | tar -xvJf - ; echo "Decompressed: $basename"' _ {}
```

### Setup Environment

The environment can be set up using the following commands (python version 3.10+):
using conda fro python 3.10, then everything works (for python 3.12 problems with torch version)

```bash
python3 -m venv venv3
source venv3/bin/activate
pip install --no-cache-dir -r requirements.txt
pip install --no-cache-dir flash-attn==2.3.4
pip install -e .
```

### Evaluating trained models

Run on a machine with one GPU (12GB+ VRAM)

#### Homology-Aware Finetuning

```bash
python evaluate_RNAformer.py -c 1 -n state_dict.pth -m models/homology_aware_finetune
```
```
sbatch --job-name "eval_ft_non_homolog" run_eval.sh ft_non_homolog
```
#### AF3-Like Finetuning

```bash
python evaluate_RNAformer.py -c 1 -n state_dict.pth -m models/af3_like_finetune
``` 
```
sbatch --job-name "eval_ft_homolog" run_eval.sh ft_homolog
```
#### Base Model

```bash
python evaluate_RNAformer.py -b -c 6 -n state_dict.pth -m models/base_model_pretrain
```
```
sbatch --job-name "eval_pretrain_base_model" run_eval.sh pretrain
```
#### Synthetic Data Trained Model

```bash
python evaluate_RNAformer.py -y -c 6 -n state_dict.pth -m models/synthetic_pretrain
``` 

```
sbatch --job-name "eval_biosynthetic_from_checkpoint" run_eval.sh bio
```

### Train biophysical Model

Run on a node with 4x A100 GPUs (40GB VRAM)

```bash
python pretrain_RNAformer.py -c train_biophysical_model
```

```
sbatch --job-name "biophysical_single_a100_4" --time "72:00:00" run_train.sh bio
```

### Pretrain model

Run on a node with 4x A100 GPUs (40GB VRAM)

```bash
python pretrain_RNAformer.py -c pretraining_base_model
``` 

```
sbatch --job-name "pretrain_single_a100_4" --time "48:00:00" run_train.sh pretrain
```

### Finetuning models

Run on a node with 4x A100 GPUs (40GB VRAM)

#### Homology-Aware Finetuning

```bash
python finetune_RNAformer.py -c finetune_homology_aware
```

```
sbatch --job-name "ft_non_homolog_single_a100_4_seed3" run_train.sh ft_non_homolog
```

#### AF3-Like Finetuning

```bash
python finetune_RNAformer.py -c finetune_af3_like
```

```
sbatch --job-name "ft_homolog_single_a100_4_seed3" run_train.sh ft_homolog
```
-->
## Inference

Create **datasets** folder with .pkl files.

Create **models** folder with subfolders containing config.yml and state_dict.pth files.

Create conda environment from file:
```
conda env create -f env.yml
```

and activate it:
```
conda activate infer
```

Run 
```
python infer_RNAformer.py -m "models/<model-name>" -p "datasets/<dataset-name>.pkl" -c 1
```

or 

```
sbatch --job-name "test_infer" run_infer.sh test
```
