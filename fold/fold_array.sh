#!/bin/bash
# AF3 inference over a list of prepared input JSON directories (MSAs/templates already in the JSON).
#
# Usage: sbatch --array=0-<ntasks-1> --export=CRC_ROOT=<crc> fold/fold_array.sh <input_root> <output_dir> <list_file> [per_task]
#   list_file : one input subdirectory (relative to input_root) per line
#   per_task  : how many list entries one array task folds sequentially (default 1)
#SBATCH --job-name=crc_fold
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=300GB
#SBATCH --time=06:00:00
#SBATCH --output=%x-%A-%a.out

set -u
source "$CRC_ROOT/config.sh"
IN_ROOT=$1; OUT_DIR=$2; LIST=$3; K=${4:-1}
mkdir -p "$OUT_DIR"
type module > /dev/null 2>&1 && module load system/singularity   # HPC environment modules (optional)
cd "$AF3_HOME"

first=$((SLURM_ARRAY_TASK_ID * K))
for i in $(seq "$first" $((first + K - 1))); do
  sub=$(sed -n "$((i + 1))p" "$LIST"); [ -z "$sub" ] && continue
  echo "[$(date)] task $SLURM_ARRAY_TASK_ID item $i: $sub"
  start=$(date +%s)
  singularity exec --nv \
    --bind "$IN_ROOT/$sub":/root/af_input \
    --bind "$OUT_DIR":/root/af_output \
    --bind "$AF3_MODELS":/root/models \
    --bind "$AF3_DB":/root/public_databases \
    "$AF3_SIF" python alphafold3/run_alphafold.py \
      --input_dir=/root/af_input --model_dir=/root/models --db_dir=/root/public_databases \
      --output_dir=/root/af_output --run_inference=True --run_data_pipeline=False
  echo "FOLD_DONE $sub rc=$? runtime=$(( $(date +%s) - start ))s"
done
