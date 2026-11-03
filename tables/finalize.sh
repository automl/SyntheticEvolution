#!/bin/bash
# Regenerate every table, the flat metrics export, all rebuttal/CRC analyses and checksums from the metric DBs
# and the unzipped download bundles (see README: "Reproduce the numbers").
# Usage: bash tables/finalize.sh [--full]     (--full adds the structural-novelty scan, ~32k US-align pairs)
#        sbatch tables/finalize.sh [--full]
#SBATCH --job-name=crc_finalize
#SBATCH --partition=cpu-single
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
set -uo pipefail
CRC_ROOT=${CRC_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
[ -f "$CRC_ROOT/config.sh" ] || CRC_ROOT=${SLURM_SUBMIT_DIR:-$PWD}   # sbatch runs a copy of this script
export CRC_ROOT; source "$CRC_ROOT/config.sh"
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
FULL=0; [ "${1:-}" = "--full" ] && FULL=1
PY=$EVAL_PY
for d in "$HOME/.conda/envs/$RNABENCH_ENV/bin" "${CONDA_PREFIX:-/nonexistent}/../$RNABENCH_ENV/bin"; do
  [ -x "$d/cd-hit-est-2d" ] && export PATH="$PATH:$d"   # cd-hit-est-2d for the sequence-novelty analysis
done
cd "$CRC_ROOT"
[ -x "${USALIGN:-}" ] || export USALIGN=$CRC_ROOT/eval/engine/evaluation/USalign/USalign
[ -x "$USALIGN" ] || bash eval/build_usalign.sh
LOG=results/analysis/logs; mkdir -p "$LOG"
run() {  # run <dir> <script> [args...]
  local d=$1 s=$2; shift 2; local tag=$s$( (($#)) && printf "_%s" "$@" | tr -d "-")${LOGSUF:-}
  (cd "$d" && $PY "$s.py" "$@") > "$LOG/$tag.log" 2>&1 && echo "ok   $tag" || echo "FAIL $tag (see $LOG/$tag.log)"
}
$PY eval/check_completeness.py
$PY tables/make_crc_tables.py --min-preds 20
$PY tables/make_crc_tables.py --submitted
$PY tables/export_metrics.py
for s in multiseed_orphan_parity uTY4_q4_near_floor multiseed_nomsa_final multiseed_nomsa_cis denoise_ms \
         general2_long_run_ladder general2_runtime robustness_table across_samples orphan_gated orphan_verify \
         struct_vs_delta confidence_ladder confidence_ladder_confirm; do
  run analysis/rebuttal $s
done
if command -v cd-hit-est-2d > /dev/null; then run analysis/rebuttal faN3_q1_similarity --mode sequence
else echo "skip faN3_q1_similarity sequence (cd-hit-est-2d not found; RnaBench env)"; fi
run analysis/rebuttal faN3_q1_similarity --mode structure --subtype single_chain
[ $FULL = 1 ] && run analysis/rebuttal faN3_q1_similarity --mode structure --subtype complex
run analysis/construction eval_intervention_ms
for n in 1000 100; do   # 91-RNA construction sweep at SHS depth N=1000 (original) and N=100 (paper default)
  for s in construct91_eval construct91_within construct91_plot; do CONSTRUCT_N=$n LOGSUF=_N$n run analysis/construction $s; done
done
bash provenance/make_manifests.sh
echo FINALIZE_DONE
