#!/bin/bash
# Re-score the six RhoFold+ runs (relaxed models) with the same engine as AF3 -> $EVAL_DB_DIR/<tag>_single_s1_m0.db
# (default results/db_rhofold). RhoFold+ is deterministic: one run per input, no seeds.
# Usage: bash eval/rhofold/score_rhofold.sh   (submits one Slurm array of 6 shards)
set -euo pipefail
CRC_ROOT=${CRC_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}; export CRC_ROOT
source "$CRC_ROOT/config.sh"
SRC=${RHOFOLD_DIR:-$CRC_OUTPUTS/rhofold}                 # predictions_rhofold_* folders (outputs bundle)
LAY=$WORK/rhofold_layout; SH=$WORK/rhofold_shards.tsv; mkdir -p "$LAY" "$WORK/logs"; : > "$SH"
while read -r tag dir; do
  "$EVAL_PY" "$CRC_ROOT/eval/rhofold/rhofold_to_eval_layout.py" --src "$SRC/$dir" --dest "$LAY/$tag" --tag "$tag" | grep -v '^#' >> "$SH"
done <<'EOF'
rnafoldN100Rhofold predictions_rhofold_on_shs_from_rnafold_folder
spotrnaN100Rhofold predictions_rhofold_on_shs_from_spotrna_folder
rnaformerN100Rhofold predictions_rhofold_on_shs_from_rnaformerN100_folder
rnafoldN100DenoiseRhofold predictions_rhofold_on_shs_from_rna_monomers_rnafold_denoise_folder
dssrN100Rhofold predictions_rhofold_on_shs_from_rna_monomers_dssrN100_folder
stdmsaRhofold predictions_rhofold_stdmsa_12636910
EOF
sbatch --array=0-5 --job-name=crc_eval_rhofold --export=CRC_ROOT="$CRC_ROOT",EVAL_DB_DIR="${EVAL_DB_DIR:-$CRC_ROOT/results/db_rhofold}" \
  --output="$WORK/logs/%x-%A-%a.out" "$CRC_ROOT/eval/run_eval_array.sh" "$SH"
