#!/bin/bash
# CRC evaluation shards (DBs land in $DB_DIR). Seeds 0/42/137 of the original conditions are already scored
# (results/db, copied from the rebuttal run); this adds seed 1 everywhere plus the new CRC folds.
#   Stage A: seed 1 for conditions whose seed-1 folds exist           (submit now)
#   Stage B: conditions folded for the CRC (after the fold arrays)    (needs fold job ids)
# Usage: bash eval/runs/submit_crc_evals.sh <natmsa_fold> <n10k_single_fold> <n5k_cplx_fold> <n20k_cplx_fold> <n10k_cplx_fold>
set -euo pipefail
CRC_ROOT=$(cd "$(dirname "$0")/../.." && pwd); export CRC_ROOT
source "$CRC_ROOT/config.sh"
R=$CRC_ROOT/eval/runs; L=$WORK/logs; mkdir -p "$L"
PY=$HOME/.conda/envs/$EVAL_ENV/bin/python
MS="$PY $CRC_ROOT/eval/make_shards.py"
EV="sbatch --parsable --export=CRC_ROOT=$CRC_ROOT $CRC_ROOT/eval/run_eval_array.sh"
sub() {  # sub <name> <shards.tsv> [dependency]
  local n; n=$(wc -l < "$2"); [ "$n" -eq 0 ] && { echo "$1: nothing to do"; return; }
  local dep=""; [ -n "${3:-}" ] && dep="--dependency=afterany:$3"
  local j; j=$(sbatch --parsable --export=CRC_ROOT="$CRC_ROOT" --job-name="crc_eval_$1" --array=0-$((n - 1))%40 \
    --output="$L/%x-%A-%a.out" $dep "$CRC_ROOT/eval/run_eval_array.sh" "$2")
  echo "$1: $n shards -> job $j" | tee -a "$R/submitted_jobs.txt"
}

# Stage A
$MS --out "$R/stageA_seed1.tsv" --seeds 1 --missing --tags empty rnafold spotrna rnaformer dssr denoise_rnafold \
  randomSHS shuffleSHS nat_and_random_N100 nat_and_shuffle_N100 synMSA_N20 synMSA_N50 synMSA_N100 synMSA_N200 \
  synMSA_N500 synMSA_N1000 rnaformerN5000 rnaformerN20000
sub seed1_existing "$R/stageA_seed1.tsv"

# Stage B
[ $# -eq 5 ] || { echo "stage B needs the 5 fold job ids"; exit 0; }
$MS --out "$R/stageB_af3_seed1.tsv" --seeds 1 --missing --tags af3;                                   sub af3_seed1 "$R/stageB_af3_seed1.tsv" "$1"
$MS --out "$R/stageB_N10000_single_seed1.tsv" --seeds 1 --missing --tags rnaformerN10000 --idsets single; sub N10000_single_seed1 "$R/stageB_N10000_single_seed1.tsv" "$2"
$MS --out "$R/stageB_N5000_complex.tsv" --seeds 0 42 137 --missing --tags rnaformerN5000 --idsets complex; sub N5000_complex "$R/stageB_N5000_complex.tsv" "$3"
$MS --out "$R/stageB_N20000_complex.tsv" --seeds 0 42 137 --missing --tags rnaformerN20000 --idsets complex; sub N20000_complex "$R/stageB_N20000_complex.tsv" "$4"
$MS --out "$R/stageB_N10000_complex.tsv" --seeds 1 0 42 137 --missing --tags rnaformerN10000 --idsets complex; sub N10000_complex "$R/stageB_N10000_complex.tsv" "$5"
