#!/bin/bash
# sha256 manifests of everything the results depend on (written to provenance/manifests/).
set -euo pipefail
export LC_ALL=C   # locale-independent sort order
CRC_ROOT=${CRC_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}; source "$CRC_ROOT/config.sh"
M=$CRC_ROOT/provenance/manifests; mkdir -p "$M"
cd "$CRC_ROOT"
find data -type f ! -name '.*' ! -path '*/__pycache__/*' | sort | xargs sha256sum > "$M/data.sha256"
find -H "$DB_DIR" "$DB_LONG" "$DB_RHOFOLD" -name '*.db' | sort | xargs sha256sum | sed "s|  $DB_DIR/|  db/|; s|  $DB_LONG/|  db_long/|; s|  $DB_RHOFOLD/|  db_rhofold/|" > "$M/metric_dbs.sha256"
(cd "$GT_SINGLE" && sha256sum -- *.cif) > "$M/ground_truth_single_91.sha256"
(cd "$GT_COMPLEX" && sha256sum -- *.cif) > "$M/ground_truth_complex_377.sha256"
(cd "$GT_LONG" && sha256sum -- *.cif) > "$M/ground_truth_long_25.sha256"
sha256sum eval/engine/evaluation/USalign/USalign.cpp > "$M/usalign.sha256"
[ -f eval/engine/evaluation/USalign/USalign ] && sha256sum eval/engine/evaluation/USalign/USalign >> "$M/usalign.sha256"
wc -l "$M"/*.sha256
