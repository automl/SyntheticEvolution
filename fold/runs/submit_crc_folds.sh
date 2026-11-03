#!/bin/bash
# CRC runs added after the rebuttal (4 AF3 seeds x 5 samples for every main-table condition).
#   1. natural-MSA AF3, seed 1, 468 targets (same input JSONs/MSAs/templates as seeds 0/42/137)
#   2. SHS RNAformer N=10000 single chains, seed 1 (91)
#   3. depth scaling on RNA-protein complexes: N=5000 / N=20000 seeds 0/42/137 (seed 1 exists),
#      N=10000 all four seeds (SHS inputs generated here first)
# Usage: bash fold/runs/submit_crc_folds.sh   (prints job ids; writes fold/runs/*.list)
set -euo pipefail
CRC_ROOT=$(cd "$(dirname "$0")/../.." && pwd); export CRC_ROOT
source "$CRC_ROOT/config.sh"
R=$CRC_ROOT/fold/runs; L=$WORK/logs; mkdir -p "$L"
IDS=$CRC_ROOT/data/ids
SB="sbatch --parsable --export=CRC_ROOT=$CRC_ROOT"

# ---- input preparation (CPU) -------------------------------------------------------------
prep=$($SB --job-name=crc_prep_inputs --partition=cpu-single --cpus-per-task=2 --mem=32G --time=06:00:00 \
  --output=$L/%x-%j.out --wrap "
$HOME/.conda/envs/$EVAL_ENV/bin/python $CRC_ROOT/fold/prepare_inputs.py --src_root $AF_INPUT/nature_methods_full_evaluation_data_files_af3_rerun_multiseeds \
  --ids $IDS/all_468.txt --seeds 1 --dst_root $AF_INPUT/crc2026_natmsa_seed1 --list $R/natmsa_seed1.list
$HOME/.conda/envs/$EVAL_ENV/bin/python $CRC_ROOT/fold/prepare_inputs.py --src_root $AF_INPUT/rnaonly_data_nature_methods_synMSA_config80_lsHPO_N10000_rnaformer_multiseed \
  --ids $IDS/single_91.txt --seeds 1 --dst_root $AF_INPUT/crc2026_rnaonly_N10000_rnaformer_seed1 --list $R/rnaonly_N10000_seed1.list")
for N in 5000 20000; do
  python3 "$CRC_ROOT/fold/prepare_inputs.py" --src_root "$AF_INPUT/all_data_nature_methods_synMSA_config80_lsHPO_N${N}_rnaformer_multiseeds" \
    --ids "$IDS/complex_377.txt" --list "$R/complex_N${N}_seeds_0_42_137.list"
done
# N=10000 complex SHS inputs: generate (38 tasks x 10 targets), then set seeds 1/0/42/137
gen=$($SB --job-name=crc_shs_N10000_complex --array=0-37 --output=$L/%x-%A-%a.out \
  "$CRC_ROOT/shs/generate_shs.sh" "$AF_INPUT/nature_methods_full_evaluation_data_files" "$IDS/complex_377.txt" \
  "$AF_INPUT/crc2026_complex_N10000_rnaformer" 10000 10)
seed10k=$($SB --job-name=crc_seed_N10000_complex --dependency=afterok:$gen --partition=cpu-single --cpus-per-task=2 \
  --mem=32G --time=04:00:00 --output=$L/%x-%j.out --wrap "
$HOME/.conda/envs/$EVAL_ENV/bin/python $CRC_ROOT/fold/prepare_inputs.py --src_root $AF_INPUT/crc2026_complex_N10000_rnaformer \
  --ids $IDS/complex_377.txt --seeds 1 0 42 137 --dst_root $AF_INPUT/crc2026_complex_N10000_rnaformer_multiseed \
  --list $R/complex_N10000_seeds_1_0_42_137.list")

# ---- AF3 folds (GPU) -------------------------------------------------------------------------
F="$CRC_ROOT/fold/fold_array.sh"
f1=$($SB --job-name=crc_fold_natmsa_s1 --dependency=afterok:$prep --array=0-467 --time=04:00:00 \
  --output=$L/%x-%A-%a.out "$F" "$AF_INPUT/crc2026_natmsa_seed1" "$AF_OUTPUT/crc2026_natmsa_seed1" "$R/natmsa_seed1.list" 1)
f2=$($SB --job-name=crc_fold_N10000_single_s1 --dependency=afterok:$prep --array=0-90 --time=01:00:00 \
  --output=$L/%x-%A-%a.out "$F" "$AF_INPUT/crc2026_rnaonly_N10000_rnaformer_seed1" "$AF_OUTPUT/crc2026_rnaonly_N10000_rnaformer_seed1" "$R/rnaonly_N10000_seed1.list" 1)
f3=$($SB --job-name=crc_fold_N5000_complex --array=0-188 --time=12:00:00 --output=$L/%x-%A-%a.out "$F" \
  "$AF_INPUT/all_data_nature_methods_synMSA_config80_lsHPO_N5000_rnaformer_multiseeds" "$AF_OUTPUT/crc2026_complex_N5000_rnaformer_multiseed" "$R/complex_N5000_seeds_0_42_137.list" 2)
f4=$($SB --job-name=crc_fold_N20000_complex --array=0-188 --time=12:00:00 --output=$L/%x-%A-%a.out "$F" \
  "$AF_INPUT/all_data_nature_methods_synMSA_config80_lsHPO_N20000_rnaformer_multiseeds" "$AF_OUTPUT/crc2026_complex_N20000_rnaformer_multiseed" "$R/complex_N20000_seeds_0_42_137.list" 2)
f5=$($SB --job-name=crc_fold_N10000_complex --dependency=afterok:$seed10k --array=0-188 --time=16:00:00 \
  --output=$L/%x-%A-%a.out "$F" "$AF_INPUT/crc2026_complex_N10000_rnaformer_multiseed" "$AF_OUTPUT/crc2026_complex_N10000_rnaformer_multiseed" "$R/complex_N10000_seeds_1_0_42_137.list" 2)

printf 'prep=%s gen=%s seed10k=%s natmsa=%s n10k_single=%s n5k_cplx=%s n20k_cplx=%s n10k_cplx=%s\n' \
  "$prep" "$gen" "$seed10k" "$f1" "$f2" "$f3" "$f4" "$f5" | tee -a "$R/submitted_jobs.txt"
