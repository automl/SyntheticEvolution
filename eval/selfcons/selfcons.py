#!/usr/bin/env python3
"""Self-consistency of SHS-guided AF3 predictions (single chains), for one seed/sample.

Self-consistency = MCC between the topology injected into the SHS and the base pairs DSSR extracts from the
SHS-guided AF3 prediction; no ground truth is used. The injected topology is the RNAformer prediction exactly as the SHS
generator encodes it (shs/rna_msa_generator_base_pair.py: lookup in rnaformer_all_preds.pkl by PDB id, else by sequence,
then convert_pred_pairs to a partner map). Pair matrices are symmetric L x L (L = query length), all DSSR pairs are kept,
and TP/TN/FP/FN are counted over the full matrix.
DSSR (x3dna-dssr) is not distributed; pass the path to your own licensed binary.

  python eval/selfcons/selfcons.py <af3_output_dir> <af3_input_dir> <seed> <sample> <ids.txt> <rnaformer_all_preds.pkl> <dssr_bin> <out.tsv.gz>
The query sequence (matrix size) is read from the AF3 input JSON.
"""
import glob, gzip, json, os, subprocess, sys, tempfile
import numpy as np
import pandas as pd

exp, inp, seed, sample, ids_file, pairs_pkl, dssr, out = sys.argv[1:9]
preds = pd.read_pickle(pairs_pkl)


def injected_pairs(pdb_id, seq):
    """Partner map {i: j, j: i} as built by the SHS generator (same lookup and convert_pred_pairs)."""
    d = preds[preds["pdb_id"].str.lower() == pdb_id]
    if d.empty: d = preds[preds["sequence"] == seq]
    if d.empty: return None
    partner = {}
    for p1, p2 in d["pairs"].values[0]:
        partner[int(p1)] = int(p2); partner[int(p2)] = int(p1)
    return list(partner.items())


def query_sequence(i):
    js = glob.glob(f"{inp}/{i.upper()}/*.json") + glob.glob(f"{inp}/{i}/*.json")
    if not js: return None
    for s in json.load(open(js[0]))["sequences"]:
        if "rna" in s: return s["rna"]["sequence"].upper()
    return None


def mat(pairs, n):
    m = np.zeros((n, n), dtype=int)
    for i, j in pairs:
        if i < n and j < n: m[i, j] = m[j, i] = 1
    return m


def mcc(pred, true):
    tp = np.logical_and(pred, true).sum(); tn = np.logical_and(pred == 0, true == 0).sum()
    fp = pred.sum() - tp; fn = true.sum() - tp
    return (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)


def dssr_pairs(cif, workdir):
    out_json = os.path.join(workdir, "out.json")
    subprocess.run([dssr, f"-i={cif}", "--json", f"-o={out_json}"], cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    data = json.load(open(out_json))
    if "pairs" not in data or "nts" not in data: return []
    idx = {nt["nt_id"]: int(nt["index_chain"]) for nt in data["nts"]}
    return [[idx[p["nt1"]] - 1, idx[p["nt2"]] - 1] for p in data["pairs"]]


rows = []
for i in open(ids_file).read().split():
    cifs = sorted(glob.glob(f"{exp}/{i}_*/seed-{seed}_sample-{sample}/model.cif")
                  + glob.glob(f"{exp}/{i.upper()}_*/seed-{seed}_sample-{sample}/model.cif"))
    seq = query_sequence(i)
    injected = injected_pairs(i, seq) if seq else None
    if not cifs or injected is None:
        rows.append((i, "MISSING", "")); continue
    try:
        with tempfile.TemporaryDirectory() as wd:
            pred = dssr_pairs(cifs[0], wd)
        rows.append((i, "ok", f"{mcc(mat(pred, len(seq)), mat(injected, len(seq))):.6f}"))
    except Exception as e:
        rows.append((i, f"ERROR:{type(e).__name__}:{e}".replace("\t", " ")[:200], ""))

tmp = out + ".tmp"
with gzip.open(tmp, "wt") as f:
    f.write("id\tstatus\tselfcons_mcc\n")
    for r in rows: f.write("\t".join(r) + "\n")
os.replace(tmp, out)
print(f"s{seed} m{sample}: {sum(r[1] == 'ok' for r in rows)}/{len(rows)} ok -> {out}")
