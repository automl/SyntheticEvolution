#!/usr/bin/env python3
"""Recompute the SPOT-RNA / RNAfold topologies that shs/rna_msa_generator_base_pair.py computes on the fly.

Same calls as the generator (RnaBench wrappers): RNAfold = RnaBench RNAFold (ViennaRNA RNA.fold -> db2pairs);
SPOT-RNA = RnaBench preparing_tfr_records + infer_spotrna (5-model ensemble, CPU), with a private tfrecords file per
call so that array tasks do not share the wrapper's fixed temporary path. Pairs are written in predictor output order
(the generator keeps the last partner per nucleotide). Query sequences are read from the shipped AF3 inputs of the
matching SHS condition (<method>_<idset>), i.e. exactly the sequence the generator saw.

Run with the RnaBench env (config.sh: RNABENCH_ENV, RNABENCH_HOME):
  python analysis/topologies/predict_2d.py --method rnafold --set single              # all targets, one process
  python analysis/topologies/predict_2d.py --method spotrna --set complex1v1 --task 3 --ntasks 20   # one shard
  python analysis/topologies/predict_2d.py --method spotrna --set complex1v1 --merge  # merge shards
Output: data/topologies/<method>_<set>.json (list of {"pdb_id", "sequence", "pairs"}, 0-based).
"""
import argparse, glob, json, os, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import CRC, INPUTS, MC_CSV  # noqa: E402

sys.path.insert(0, os.environ.get("RNABENCH_HOME") or os.path.expanduser("~/RnaBench"))
OUT = CRC / "data" / "topologies"


def targets(which):
    if which == "single":
        ids = (CRC / "data" / "ids" / "single_91.txt").read_text().split()
    else:
        import csv
        ids = sorted({r["exp_db_id"].upper() for r in csv.DictReader(open(MC_CSV)) if r["group"] == "1v1"})
    return [i.upper() for i in ids]


def query_sequence(method, idset, pid):
    js = glob.glob(str(INPUTS / "af3_inputs" / f"{method}_{idset}" / pid / "*.json"))
    if not js: return None
    return [s["rna"]["sequence"] for s in json.load(open(js[0]))["sequences"] if "rna" in s][0]


def predictor(method):
    if method == "rnafold":
        from RnaBench.lib.rna_folding_algorithms.rnafold import RNAFold
        f = RNAFold()
        return lambda seq: [[int(i), int(j)] for i, j, *_ in f(seq)[0]]
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    from RnaBench.lib.rna_folding_algorithms.DL.spotrna import preparing_tfr_records, infer_spotrna
    spot = os.path.join(os.environ.get("RNABENCH_HOME") or os.path.expanduser("~/RnaBench"), "external_algorithms", "SPOT-RNA")
    def run(seq):
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
            tfr = os.path.join(d, "spotrna.tfrecords")
            preparing_tfr_records(seq, tfr)
            return [[int(p[0]), int(p[1])] for p in infer_spotrna(seq, tfr, spot, gpu=-1, cpu=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))]
    return run


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", choices=["rnafold", "spotrna"], required=True)
    ap.add_argument("--set", choices=["single", "complex1v1"], required=True)
    ap.add_argument("--task", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args(); idset = "single" if a.set == "single" else "complex"
    OUT.mkdir(parents=True, exist_ok=True); shard_dir = OUT / "shards" / f"{a.method}_{a.set}"
    final = OUT / f"{a.method}_{a.set}.json"
    if a.merge:
        recs = [r for f in sorted(shard_dir.glob("*.json")) for r in json.load(open(f))]
        want = targets(a.set); got = {r["pdb_id"] for r in recs}
        missing = [i for i in want if i not in got]
        if missing: sys.exit(f"missing {len(missing)} targets, e.g. {missing[:5]}")
        json.dump(sorted(recs, key=lambda r: r["pdb_id"]), open(final, "w"))
        print(f"merged {len(recs)} -> {final}"); return
    ids = targets(a.set)[a.task::a.ntasks]
    run = predictor(a.method); recs = []
    for pid in ids:
        seq = query_sequence(a.method, idset, pid)
        if seq is None: print(f"  {pid}: no AF3 input for {a.method}_{idset}"); continue
        recs.append({"pdb_id": pid, "sequence": seq, "pairs": run(seq)})
        print(f"  {pid}: {len(seq)} nt, {len(recs[-1]['pairs'])} pairs", flush=True)
    if a.ntasks == 1:
        json.dump(recs, open(final, "w")); print(f"wrote {len(recs)} -> {final}")
    else:
        shard_dir.mkdir(parents=True, exist_ok=True)
        json.dump(recs, open(shard_dir / f"task{a.task:03d}.json", "w")); print(f"wrote shard {a.task}: {len(recs)}")


if __name__ == "__main__":
    main()
