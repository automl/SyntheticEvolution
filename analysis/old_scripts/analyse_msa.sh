#! /bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <algorithm> <data_dir> <out_dir>"
    exit 1
fi

python get_gapped_alignment.py --algorithm $1 --data_dir $2 --out_dir $3
