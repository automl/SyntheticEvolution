#!/bin/bash 
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --gres=gpu:A100:4
#SBATCH --time=38:00:00
#SBATCH --mem=250gb
#SBATCH --output=LOGS//%x.%N.%A.%a.out
#SBATCH --error=LOGS//%x.%N.%A.%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matus.dominika@gmail.com

source venv3/bin/activate
# source /home/fr/fr_fr/fr_dm339/miniconda3/bin/activate iclr

case "$1" in
    "bio")
        python pretrain_RNAformer.py -c train_biophysical_model
        ;;
    "bio_32")
        python pretrain_RNAformer.py -c train_bio_bs32
        ;;
    "bio_grad_acc")
        python pretrain_RNAformer.py -c train_bio_grad_acc
        ;;
    "bio_grad_acc_2m")
        python pretrain_RNAformer.py -c train_bio_grad_acc_2m
        ;;
    "bio_grad_acc_8m")
        python pretrain_RNAformer.py -c train_bio_grad_acc_8m
        ;;
    "pretrain")
        python pretrain_RNAformer.py -c pretraining_base_model
        ;;
    "ft_non_homolog")
        python finetune_RNAformer.py -c finetune_homology_aware
        ;;
    "ft_homolog")
        python finetune_RNAformer.py -c finetune_af3_like
        ;;
    *)
        ech/work/ws/fr_dm339-rnaformer/RNAformer_ICLR25/statso "Unknown script: $1"
        exit 1
        ;;
esac