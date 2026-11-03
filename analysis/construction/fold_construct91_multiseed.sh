#!/bin/bash
# Fold the 91-RNA construction sweep (91 task dirs x 4 steps, seeds 1/0/42/137 x 5 samples) with the generic
# fold/fold_array.sh. Run from the repository root after build_inputs_91.sh.
# Usage: bash analysis/construction/fold_construct91_multiseed.sh [N]   (default 1000 = original run)
CRC_ROOT=${CRC_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}; source "$CRC_ROOT/config.sh"
N=${1:-1000}
[ "$N" = 1000 ] && SUB=construct91_multiseed || SUB=construct91_N${N}_multiseed
IN=$AF_INPUT/$SUB
OUT=$AF_OUTPUT/$SUB
LIST=$CRC_ROOT/results/analysis/construction/${SUB}_fold.list
seq 0 90 > "$LIST"
sbatch --array=0-90 --time=04:00:00 --mem=100GB --cpus-per-task=4 --job-name=$SUB ${SBATCH_EXTRA:-} \
  --export=CRC_ROOT="$CRC_ROOT" "$CRC_ROOT/fold/fold_array.sh" "$IN" "$OUT" "$LIST" 1
