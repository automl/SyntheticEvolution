#!/bin/bash
# Quick single-seed check of the 4P8Z construction sweep: mean RNA TM-score / RMSD per step (all depths x samples).
# Usage: bash analysis/construction/construction_plateau.sh   (needs US-align: eval/build_usalign.sh)
CRC_ROOT=${CRC_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}; source "$CRC_ROOT/config.sh"
UA=${USALIGN:-$CRC_ROOT/eval/engine/evaluation/USalign/USalign}
GT=$CRC_ROOT/data/gt_intervention/4p8z_gt.cif
D=${CRC_OUTPUTS}/af3_case_studies/4P8Z_intervention_full_sweep
[ -d "$D" ] || D=$AF_OUTPUT/4P8Z_intervention_full_sweep
echo "construction sweep: mean RNA TM-score vs 4P8Z GT (single-seed, all depths x samples)"
printf "%-14s %5s %8s %8s\n" "step" "n" "mean_TM" "mean_RMSD"
for step in local_wc init_wc step1_nc step2_pk step3_multi final_all; do
  tms=(); rms=()
  for p in $(find "$D"/ -path "*intervention_${step}/*" -name model.cif 2>/dev/null); do
    out=$("$UA" "$p" "$GT" -mol RNA -ter 0 2>/dev/null)
    tm=$(echo "$out" | grep -m1 "normalized by length of Structure_2" | grep -oE "TM-score= [0-9.]+" | grep -oE "[0-9.]+")
    rm=$(echo "$out" | grep -m1 "Aligned length" | grep -oE "RMSD= *[0-9.]+" | grep -oE "[0-9.]+")
    [ -n "$tm" ] && tms+=("$tm") && rms+=("$rm")
  done
  n=${#tms[@]}
  mtm=$(printf "%s\n" "${tms[@]}" | awk '{s+=$1}END{if(NR)printf "%.3f",s/NR}')
  mrm=$(printf "%s\n" "${rms[@]}" | awk '{s+=$1}END{if(NR)printf "%.2f",s/NR}')
  printf "%-14s %5s %8s %8s\n" "$step" "$n" "$mtm" "$mrm"
done
echo "DONE"
