#!/bin/bash 
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=10:00
#SBATCH --mem=50gb
#SBATCH --output=LOGS//%x.%N.%A.%a.out
#SBATCH --error=LOGS//%x.%N.%A.%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matus.dominika@gmail.com

# source venv3/bin/activate
# source /home/fr/fr_fr/fr_dm339/miniconda3/bin/activate iclr
source venv3/bin/activate

case "$1" in
    "2m")
        python evaluate_RNAformer.py -y -c 6 -n state_dict.pth -m new_models/biosynthetic_2m_cyc_6
        ;;
    "8m")
        python evaluate_RNAformer.py -y -c 6 -n state_dict.pth -m new_models/biosynthetic_8m_cyc_6
        ;;
    "bio")
        python evaluate_RNAformer.py -y -c 6 -n state_dict.pth -m models/synthetic_pretrain
        ;;
    "pretrain")
        python evaluate_RNAformer.py -b -c 8 -n state_dict.pth -m models/base_model_pretrain
        ;;
    "ft_non_homolog")
        python evaluate_RNAformer.py -c 1 -n state_dict.pth -m models/homology_aware_finetune
        ;;
    "ft_homolog")
        python evaluate_RNAformer.py -c 1 -n state_dict.pth -m models/af3_like_finetune
        ;;
    *)
        echo "Unknown script: $1"
        exit 1
        ;;
esac