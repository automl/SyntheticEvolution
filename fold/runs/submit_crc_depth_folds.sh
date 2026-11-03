#!/bin/bash
# CRC depth-scaling completion (SHS RNAformer only, config80, 4 AF3 seeds x 5 samples on 91 single + 377 complexes):
#   N=20, N=500   seeds 0/42/137 exist (all_data ..._rnaformerN{20,500}_multiseed)  -> add seed 1 (468)
#   N=15000       single seeds 0/42/137 exist (rnaonly ..._N15000_rnaformer_multiseed) -> add seed 1 (91);
#                 complexes never folded -> all four seeds (377)
#   N=1000        never folded -> all four seeds (468)
# Inputs are the existing config80 JSONs that were (or would have been) folded; prepare_inputs.py only replaces
# `modelSeeds` (checked: config80 seed-1 base == _multiseeds inputs in every chain incl. MSAs/templates, and the
# folded N20/N500/N15000 outputs' _data.json == these inputs).
# Usage: bash fold/runs/submit_crc_depth_folds.sh   (prints job ids; writes fold/runs/*.list)
set -euo pipefail
CRC_ROOT=$(cd "$(dirname "$0")/../.." && pwd); export CRC_ROOT
source "$CRC_ROOT/config.sh"
R=$CRC_ROOT/fold/runs; L=$WORK/logs; mkdir -p "$L"
IDS=$CRC_ROOT/data/ids
SB="sbatch --parsable --export=CRC_ROOT=$CRC_ROOT"
P="$EVAL_PY $CRC_ROOT/fold/prepare_inputs.py"
ALL=$AF_INPUT/all_data_nature_methods_synMSA_config80_lsHPO

# ---- input preparation (CPU): re-seed only ----------------------------------------------------
prep=$($SB --job-name=crc_prep_depth --partition=cpu-single --cpus-per-task=2 --mem=32G --time=04:00:00 \
  --output=$L/%x-%j.out --wrap "set -e
$P --src_root ${ALL}_N20_rnaformer_multiseeds  --ids $IDS/all_468.txt --seeds 1 --dst_root $AF_INPUT/crc2026_N20_rnaformer_seed1  --list $R/N20_seed1.list
$P --src_root ${ALL}_N500_rnaformer_multiseeds --ids $IDS/all_468.txt --seeds 1 --dst_root $AF_INPUT/crc2026_N500_rnaformer_seed1 --list $R/N500_seed1.list
$P --src_root $AF_INPUT/rnaonly_data_nature_methods_synMSA_config80_lsHPO_N15000_rnaformer_multiseed \
   --ids $IDS/single_91.txt --seeds 1 --dst_root $AF_INPUT/crc2026_rnaonly_N15000_rnaformer_seed1 --list $R/rnaonly_N15000_seed1.list
$P --src_root ${ALL}_N15000_rnaformer_multiseeds --ids $IDS/complex_377.txt --seeds 1 0 42 137 \
   --dst_root $AF_INPUT/crc2026_complex_N15000_rnaformer_multiseed --list $R/complex_N15000_seeds_1_0_42_137.list
$P --src_root ${ALL}_N1000_rnaformer_multiseeds --ids $IDS/all_468.txt --seeds 1 0 42 137 \
   --dst_root $AF_INPUT/crc2026_N1000_rnaformer_multiseed --list $R/N1000_seeds_1_0_42_137.list")

# ---- AF3 folds (GPU) ----------------------------------------------------------------------------
F="$CRC_ROOT/fold/fold_array.sh"
fold() {  # fold <name> <ntargets> <per_task> <time> <in> <out> <list>
  local n=$(( ($2 + $3 - 1) / $3 - 1 ))
  $SB --job-name=$1 --dependency=afterok:$prep --array=0-$n --time=$4 --output=$L/%x-%A-%a.out "$F" "$5" "$6" "$7" "$3"
}
f1=$(fold crc_fold_N20_s1 468 4 03:00:00 "$AF_INPUT/crc2026_N20_rnaformer_seed1" "$AF_OUTPUT/crc2026_N20_rnaformer_seed1" "$R/N20_seed1.list")
f2=$(fold crc_fold_N500_s1 468 4 03:00:00 "$AF_INPUT/crc2026_N500_rnaformer_seed1" "$AF_OUTPUT/crc2026_N500_rnaformer_seed1" "$R/N500_seed1.list")
f3=$(fold crc_fold_N15000_single_s1 91 4 03:00:00 "$AF_INPUT/crc2026_rnaonly_N15000_rnaformer_seed1" \
  "$AF_OUTPUT/crc2026_rnaonly_N15000_rnaformer_seed1" "$R/rnaonly_N15000_seed1.list")
f4=$(fold crc_fold_N15000_complex 377 2 12:00:00 "$AF_INPUT/crc2026_complex_N15000_rnaformer_multiseed" \
  "$AF_OUTPUT/crc2026_complex_N15000_rnaformer_multiseed" "$R/complex_N15000_seeds_1_0_42_137.list")
f5=$(fold crc_fold_N1000 468 2 12:00:00 "$AF_INPUT/crc2026_N1000_rnaformer_multiseed" \
  "$AF_OUTPUT/crc2026_N1000_rnaformer_multiseed" "$R/N1000_seeds_1_0_42_137.list")

printf 'depth: prep=%s N20_s1=%s N500_s1=%s N15000_single_s1=%s N15000_cplx=%s N1000=%s\n' \
  "$prep" "$f1" "$f2" "$f3" "$f4" "$f5" | tee -a "$R/submitted_jobs.txt"
