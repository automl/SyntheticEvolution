#!/usr/bin/env bash
# Download, verify and unpack the data bundles of the NeurIPS 2026 SHS camera-ready package.
#
#   bash download_data.sh             # -> downloads/2026_neurips_shs_{inputs,outputs,extras}/  (zips are deleted after unpacking)
#   bash download_data.sh --keep-zips # keep the .zip files
#   CRC_DATA=/big/disk bash download_data.sh   # download somewhere else (config.sh reads the same variable)
#
# Needs ~47 GB for the zips + ~218 GB unpacked (~362k files); curl or wget, unzip, sha256sum or shasum.
# Interrupted downloads resume when the script is run again.
set -euo pipefail

BASE_URL=${CRC_DATA_URL:-https://ml.informatik.uni-freiburg.de/research-artifacts/SyntheticEvolution/2026_neurips_crc}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEST=${CRC_DATA:-$ROOT/downloads}
KEEP=0; [ "${1:-}" = "--keep-zips" ] && KEEP=1

# name  sha256  (the published .sha256 files must agree with these values)
BUNDLES="2026_neurips_shs_inputs.zip 2131fa721fa0c38dae15a957d94d08386b0eef9fd62fec381a1310ecb1463aca
2026_neurips_shs_outputs.zip 47ca85aacba9814b532d73675dab0ddff130d5dd076e45c175835176a98e276a
2026_neurips_shs_extras.zip 848283cc8741268598eeb69a41317a800757587d29198a7d491032738c532c5e"

fetch() {  # fetch <url> <file>  (resumes partial downloads)
  if command -v curl > /dev/null; then curl -fL --retry 5 -C - -o "$2" "$1"
  else wget -c --tries=5 -O "$2" "$1"; fi
}
sha256() {
  if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1
  else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

mkdir -p "$DEST"; cd "$DEST"
echo "$BUNDLES" | while read -r zip sum; do
  dir=${zip%.zip}
  if [ -f "$dir/.unpacked" ]; then echo "$dir/ already present"; continue; fi
  echo "== $zip"
  fetch "$BASE_URL/$zip.sha256" "$zip.sha256.tmp" > /dev/null 2>&1 || true
  if [ -s "$zip.sha256.tmp" ] && [ "$(cut -d' ' -f1 "$zip.sha256.tmp")" != "$sum" ]; then
    echo "published checksum of $zip differs from the one in this repository" >&2; exit 1
  fi
  rm -f "$zip.sha256.tmp"
  [ -f "$zip" ] && [ "$(sha256 "$zip")" = "$sum" ] || fetch "$BASE_URL/$zip" "$zip"
  got=$(sha256 "$zip")
  [ "$got" = "$sum" ] || { echo "checksum mismatch for $zip: $got (expected $sum); delete it and rerun" >&2; exit 1; }
  echo "checksum ok, unpacking"
  unzip -q -o "$zip"
  touch "$dir/.unpacked"
  [ $KEEP = 1 ] || rm -f "$zip"
done
echo "done: $DEST"
echo "next: export DB_DIR=$DEST/2026_neurips_shs_outputs/metric_dbs/benchmark DB_LONG=$DEST/2026_neurips_shs_outputs/metric_dbs/long_rna"
echo "      bash tables/finalize.sh"
