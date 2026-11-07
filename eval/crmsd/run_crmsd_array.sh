#!/bin/bash
# Complex RMSD (paper definition, see eval/crmsd/crmsd.py) for the complex shards of an eval shard file.
#   sbatch --array=0-<n-1> --export=CRC_ROOT=$PWD eval/crmsd/run_crmsd_array.sh <shards.tsv>
# Same shard format as eval/run_eval_array.sh (make_shards.py); non-complex lines are skipped.
# Writes results/crmsd/<tag>_complex_s<seed>_m<sample>.tsv.gz (override with CRMSD_DIR).
#SBATCH --job-name=crc_crmsd
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=04:00:00
#SBATCH --output=%x-%A-%a.out
set -u
source "$CRC_ROOT/config.sh"
SHARDS=$1
CRMSD_DIR=${CRMSD_DIR:-$CRC_ROOT/results/crmsd}
IFS=$'\t' read -r tag exp seed sample idset cls < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SHARDS")
[ "$idset" = complex ] || { echo "skip $tag $idset"; exit 0; }
[[ "$exp" = /* ]] || exp="$AF_OUTPUT/$exp"
type module > /dev/null 2>&1 && module load devel/miniforge   # HPC environment modules (optional)
conda activate "$EVAL_ENV" 2> /dev/null
export PATH="$HOME/.conda/envs/$EVAL_ENV/bin:$PATH"
mkdir -p "$CRMSD_DIR" "$WORK"
export TMPDIR=$(mktemp -d "$WORK/crmsd_tmp_XXXX")
python "$CRC_ROOT/eval/crmsd/crmsd.py" "$tag" "$exp" "$seed" "$sample" "$CRC_ROOT/data/ids/complex_377.txt" \
  "$GT_COMPLEX" "$CRMSD_DIR/${tag}_complex_s${seed}_m${sample}.tsv.gz"
rc=$?
rm -rf "$TMPDIR"
exit $rc
