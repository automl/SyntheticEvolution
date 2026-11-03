#!/bin/bash
# Build the 91-RNA construction AF3 inputs: one numbered task dir per RNA, each with 4 step JSONs
# (wc/nc/pk/all), config80 parameters, SHS depth N (default 1000 = the original run; seeds [1,0,42,137] are taken
# from the synthesized data JSONs).
# Needs build_construct_91.py output (or the frozen copy in data/construct/) and the RnaBench env.
# Usage: bash analysis/construction/build_inputs_91.sh [N] [out_dir]
#        default out_dir: $AF_INPUT/construct91_multiseed (N=1000) or $AF_INPUT/construct91_N<N>_multiseed
set -e
CRC_ROOT=${CRC_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}; source "$CRC_ROOT/config.sh"
SCR=$CRC_ROOT/results/analysis/construction          # build_construct_91.py writes construct91_datafiles/ here
PY=${RNABENCH_PY:-$HOME/.conda/envs/$RNABENCH_ENV/bin/python}
N=${1:-1000}
[ "$N" = 1000 ] && SUB=construct91_multiseed || SUB=construct91_N${N}_multiseed
INBASE=${2:-$AF_INPUT/$SUB}
export SHS_CONSTRUCT_PAIRS=$CRC_ROOT/data/construct/construct_pairs_91.json
export PYTHONPATH="$RNABENCH_HOME${PYTHONPATH:+:$PYTHONPATH}"
IDS=$CRC_ROOT/data/ids/single_91.txt
MAP=$SCR/${SUB}_idxmap.tsv                           # reference copy (N=1000): data/construct/construct91_idxmap.tsv
mkdir -p "$SCR"; : > "$MAP"
[ -e "$INBASE" ] && { echo "$INBASE exists; refusing to overwrite (pass another out_dir)"; exit 1; }
mkdir -p "$INBASE"
idx=0
while read pid; do
  [ -z "$pid" ] && continue
  pid=$(echo "$pid" | tr 'A-Z' 'a-z')
  PIDU=$(echo "$pid" | tr 'a-z' 'A-Z')
  dj=$SCR/construct91_datafiles/${pid}_data.json
  [ -f "$dj" ] || { echo "MISSING datajson $pid"; continue; }
  dst=$INBASE/$idx; mkdir -p "$dst"
  printf "%s\t%s\n" "$idx" "$pid" >> "$MAP"
  for step in wc nc pk all; do
    export SHS_CONSTRUCT_KEY=${pid}_${step}
    $PY "$CRC_ROOT/shs/rna_msa_generator_base_pair.py" \
      --structure_predictor construct --pdb_id "${PIDU}_${step}" --seed 42 -N "$N" \
      --insertion-prob-loop 0.2 --deletion-prob-loop 0.8 --insertion-prob-stem 0.005 --deletion-prob-stem 0.005 \
      --wobble-prob 0.1 --long-insertion-prob 0.1 --long-deletion-prob 0.05 \
      --max-insertion-fraction 0.1 --max-deletion-fraction 0.1 --stem-keep-prob 0.99 \
      --input_json_path "$dj" --output_json_dir "$dst" > /dev/null 2>&1
  done
  n=$(ls "$dst"/*.json 2>/dev/null | wc -l)
  [ "$n" -ne 4 ] && echo "WARN idx=$idx $pid produced $n/4"
  idx=$((idx+1))
  [ $((idx % 20)) -eq 0 ] && echo "  built $idx RNAs..."
done < "$IDS"
echo "DONE: $idx RNA task-dirs in $INBASE"
echo "total step jsons: $(find "$INBASE" -name '*.json' | wc -l) (expect $((idx*4)))"
echo "idxmap: $MAP"
