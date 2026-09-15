#!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --time=00:15:00
#SBATCH --mem=20G
set -euo pipefail
# Python paths and manifest are supplied as arguments; array size is set by caller.
module load "${3:-bio/alphafold/3.0.1}"
python "$1" "$2"
