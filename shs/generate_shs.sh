#!/bin/bash
# Generate SHS-conditioned AF3 inputs (SCE, paper configuration "config80") for a list of targets.
# Replaces the RNA unpairedMsa of each base AF3 JSON with N synthetic homologs from the RNAformer topology.
#
# Usage: sbatch --array=0-<ntasks-1> --export=CRC_ROOT=<crc> shs/generate_shs.sh <base_json_dir> <ids_file> <out_root> <N> [per_task] [predictor]
#   base_json_dir : <PDBID>_data.json AF3 inputs with natural protein MSAs/templates
#   out_root      : writes <out_root>/<PDBID>/<PDBID>_custom_rnamsa_N<N>_seed42_..._<predictor>.json
#SBATCH --job-name=crc_shs
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=%x-%A-%a.out

set -u
source "$CRC_ROOT/config.sh"
BASE=$1; IDS=$2; OUT=$3; N=$4; K=${5:-1}; PRED=${6:-rnaformer}
type module > /dev/null 2>&1 && module load devel/miniforge   # HPC environment modules (optional)
conda activate "$RNABENCH_ENV" 2> /dev/null
export PATH="$HOME/.conda/envs/$RNABENCH_ENV/bin:$PATH"   # activate can fail silently in batch shells
export PYTHONPATH="$RNABENCH_HOME${PYTHONPATH:+:$PYTHONPATH}"
python -c "import pandas, RnaBench" || { echo "RnaBench env not usable"; exit 1; }
cd "$CRC_ROOT/data/rnaformer"   # predictor 'rnaformer' reads ./rnaformer_all_preds.pkl

first=$((SLURM_ARRAY_TASK_ID * K))
for i in $(seq "$first" $((first + K - 1))); do
  pid=$(sed -n "$((i + 1))p" "$IDS" | tr 'a-z' 'A-Z'); [ -z "$pid" ] && continue
  mkdir -p "$OUT/$pid"
  python "$CRC_ROOT/shs/rna_msa_generator_base_pair.py" \
    --structure_predictor "$PRED" --pdb_id "$pid" --seed 42 -N "$N" --max_chains 1000 \
    --insertion-prob-loop 0.2 --deletion-prob-loop 0.8 \
    --insertion-prob-stem 0.005 --deletion-prob-stem 0.005 --wobble-prob 0.1 \
    --long-insertion-prob 0.1 --long-deletion-prob 0.05 \
    --max-insertion-fraction 0.1 --max-deletion-fraction 0.1 --stem-keep-prob 0.99 \
    --input_json_path "$BASE/${pid}_data.json" --output_json_dir "$OUT/$pid"
  echo "SHS_DONE $pid rc=$?"
done
