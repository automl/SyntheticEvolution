#!/bin/bash
# Build US-align (Version 20241108, source vendored) the way AF3InterfaceEval's install.sh does on Linux.
# The DBs shipped with this package were scored with such a build (sha256 recorded in provenance/).
set -e
cd "$(dirname "$0")/engine/evaluation/USalign"
"${CXX:-c++}" -O3 -ffast-math -o USalign USalign.cpp
./USalign -h | sed -n 2p
