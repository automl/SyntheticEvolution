#!/bin/bash 
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=3:00:00
#SBATCH --mem=50gb
#SBATCH --output=LOGS/%x.%N.%A.%a.out
#SBATCH --error=LOGS/%x.%N.%A.%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=runget@cs.uni-freiburg.de
#SBATCH --export=NONE

# source venv3/bin/activate
# source /home/fr/fr_fr/fr_dm339/miniconda3/bin/activate iclr

module load devel/cuda/12.1

module load devel/miniconda/3
source $MINICONDA_HOME/etc/profile.d/conda.sh
conda activate infer

# case "$1" in
#     "test")
#         python infer_RNAformer.py -c 1 -m models/af3_like_finetune -p a2021_rnaonly.pkl
#         ;;
#     *)
#         echo "Unknown script: $1"
#         exit 1
#         ;;
# esac

python infer_RNAformer.py -c 1 -m models/af3_like_finetune -p datasets/b2021_rnaonly_l500.pkl
