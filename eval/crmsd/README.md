# Complex RMSD

In the paper, the complex RMSD (cRMSD) is the RMSD over all protein and RNA atoms after PyMOL `align`, taken before
PyMOL removes outlier atoms. The evaluation engine in `eval/engine` stores the value after outlier removal in the metric
databases (`Complex_RMSD`), so we compute the paper value here.

`crmsd.py` scores one shard and writes two values per target: `rmsd_all` (the paper value) and `rmsd_refined` (the same
value as in the metric database). `run_crmsd_array.sh <shards.tsv>` runs all complex shards of an eval shard file
(`eval/make_shards.py`) as a Slurm array and writes one file per condition, seed and sample to `results/crmsd/`.

The table scripts read `rmsd_all` from `results/crmsd/` and stop with an error if a file is missing. The checks we ran
are described in `provenance/PROVENANCE.md`.
