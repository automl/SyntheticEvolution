#!/bin/bash
# Score one (condition, seed, sample) shard per array task with run_eval (basic tier:
# US-align TM/RMSD + OpenStructure all-atom LDDT / DockQ / QS / iLDDT; no web API).
#
# Usage: sbatch --array=0-<nshards-1> --export=CRC_ROOT=<crc> eval/run_eval_array.sh <shards.tsv>
# shards.tsv columns (tab): tag  af_output_dir  seed  sample  idset(single|complex|long)  class(rna_rna|protein_rna)
# Result: $EVAL_DB_DIR/<tag>_<idset>_s<seed>_m<sample>.db  (default results/db; never into the shipped bundle DBs)
#SBATCH --job-name=crc_eval
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=06:00:00
#SBATCH --output=%x-%A-%a.out

set -u
source "$CRC_ROOT/config.sh"
SHARDS=$1
EVAL_DB_DIR=${EVAL_DB_DIR:-$CRC_ROOT/results/db}
IFS=$'\t' read -r tag exp seed sample idset cls < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SHARDS")
case "$idset" in
  single)  ids=$CRC_ROOT/data/ids/single_91.txt;  gt=$GT_SINGLE ;;
  complex) ids=$CRC_ROOT/data/ids/complex_377.txt; gt=$GT_COMPLEX ;;
  long)    ids=$CRC_ROOT/data/ids/long_25.txt;    gt=$GT_LONG ;;
esac
[[ "$exp" = /* ]] || exp="$AF_OUTPUT/$exp"
name="${tag}_${idset}_s${seed}_m${sample}"
SH="$WORK/eval/$name"
echo "shard $SLURM_ARRAY_TASK_ID: $name <- $exp"

type module > /dev/null 2>&1 && module load devel/miniforge   # HPC environment modules (optional)
conda activate "$EVAL_ENV" 2> /dev/null
export PATH="$HOME/.conda/envs/$EVAL_ENV/bin:$PATH"
python -c "import pandas, pymol, Bio" || { echo "eval env not usable"; exit 1; }
rm -rf "$SH"; mkdir -p "$SH/pred" "$SH/conf" "$EVAL_DB_DIR"
n=0
while read -r id; do
  [ -z "$id" ] && continue
  cif=$(ls "$exp/${id}_"*"/seed-${seed}_sample-${sample}/model.cif" "$exp/${id^^}_"*"/seed-${seed}_sample-${sample}/model.cif" 2>/dev/null | head -1)
  [ -z "$cif" ] && { echo "MISSING $id"; continue; }
  ln -s "$cif" "$SH/pred/${id}.cif"; n=$((n + 1))
  for c in summary_confidences confidences; do
    s="$(dirname "$cif")/${c}.json"; [ -f "$s" ] && ln -s "$s" "$SH/conf/${id}_${c}.json"
  done
done < "$ids"
echo "staged $n predictions"

export EVAL_CONFIG_NAME="config_shard_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
ENGINE=$CRC_ROOT/eval/engine
cd "$ENGINE"; export PYTHONPATH="$ENGINE"
python run_eval.py --pred "$SH/pred" --gt "$gt" --confidence "$SH/conf" --out "$SH/run" \
  --classification "$cls" --metrics basic --ost-env "$OST_ENV" --non-interactive
rc=$?
if [ -f "$SH/run/run.db" ]; then
  cp -f "$SH/run/run.db" "$EVAL_DB_DIR/$name.db"
  # keep only the small provenance files; drop staged structure copies
  find "$SH/run" -mindepth 1 -maxdepth 1 ! -name 'run.db' ! -name '*.csv' -exec rm -rf {} +
  rm -rf "$SH/pred" "$SH/conf"
fi
rm -f "$ENGINE/database/$EVAL_CONFIG_NAME"
echo "SHARD_DONE $name rc=$rc staged=$n"
