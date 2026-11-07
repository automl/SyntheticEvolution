#!/usr/bin/env python3
"""2D accuracy of the RNAformer topology encoded in the SHS, per target (the gating variable rf2d_pred).

Prediction: the topology the SHS generator encodes for predictor `rnaformer`
  (data/rnaformer/rnaformer_all_preds.pkl; generator lookup by PDB id, else by sequence), verified against the shipped
  SHS by analysis/benchmark/fig3a_mcc_vs_3d.py (paired positions = conserved alignment columns).
Reference (DSSR, all base pairs):
  single chains: data/construct/construct_pairs_91.json (*_all; the topology of the reference-topology SHS)
  1v1 complexes: data/gates/inputs/rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv (pairs2)
Both in the frame of the AF3 query sequence. Metrics over the full symmetric L x L matrix (TN included), with the same
functions as fig3a_mcc_vs_3d.py, so rf2d_mcc equals the RNAformer x-axis of Fig. 3a.

  python analysis/gates/build_rnaformer_2d_gt.py                  # -> data/gates/rnaformer_2d_vs_gt.csv
  python analysis/gates/build_rnaformer_2d_gt.py --update-master  # also refresh rf2d_pred / pred_bin / has2d in
                                                                  #    data/gates/master_combined.csv
"""
import argparse, importlib.util, json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis"))
from paths import MC_CSV  # noqa: E402
_spec = importlib.util.spec_from_file_location("fig3a", ROOT / "analysis" / "benchmark" / "fig3a_mcc_vs_3d.py")
F = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(F)
OUT = ROOT / "data" / "gates" / "rnaformer_2d_vs_gt.csv"
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}


def references():
    ref = {r["pdb_id"].lower(): ("rna-only", F.parse_pairs(r["pairs"]))
           for r in json.load(open(ROOT / "data" / "construct" / "construct_pairs_91.json")) if r["id"].endswith("_all")}
    ids_1v1 = set(pd.read_csv(MC_CSV).query("group == '1v1'").exp_db_id.str.lower())
    d = pd.read_csv(ROOT / "data" / "gates" / "inputs" / "rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv")
    for r in d[d.method == "rnaformer"].itertuples():
        if r.pdbid.lower() in ids_1v1: ref[r.pdbid.lower()] = ("1v1", F.parse_pairs(r.pairs2))
    return ref


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update-master", action="store_true")
    a = ap.parse_args()
    preds = F.table_from(ROOT / "data" / "rnaformer" / "rnaformer_all_preds.pkl")
    rows = []
    for pid, (group, R) in sorted(references().items()):
        idset = "single" if group == "rna-only" else "complex"
        seq, _ = F.query_and_msa("rnaformer", idset, pid)
        if seq is None: continue
        P = F.lookup(preds, pid, seq); L = len(seq)
        p, t = F.matrix(P, L), F.matrix(R, L)
        tp = int((p & t).sum() // 2); fp = int((p & ~t).sum() // 2); fn = int((~p & t).sum() // 2)
        prec = tp / (tp + fp) if tp + fp else 0.0; rec = tp / (tp + fn) if tp + fn else 0.0
        rows.append(dict(exp_db_id=pid.upper(), group=group, rf2d_precision=prec, rf2d_recall=rec,
                         rf2d_f1=2 * prec * rec / (prec + rec) if prec + rec else 0.0, rf2d_mcc=F.mcc(p, t),
                         rf2d_degenerate=not R or not P, n_gt_pairs=int(t.sum() // 2), n_pred_pairs=int(p.sum() // 2),
                         ref_canonical_frac=np.mean([seq[i] + seq[j] in CANON for i, j in R if i < L and j < L]) if R else np.nan))
    out = pd.DataFrame(rows); out.to_csv(OUT, index=False)
    for g, d in out.groupby("group"):
        print(f"  {g:9s} n={len(d):3d} | f1 {d.rf2d_f1.mean():.3f}  mcc {d.rf2d_mcc.mean():.3f}  "
              f"precision {d.rf2d_precision.mean():.3f}  recall {d.rf2d_recall.mean():.3f}  "
              f"reference pairs complementary in query frame {d.ref_canonical_frac.mean():.3f}")
    print(f"-> {OUT}")
    if a.update_master:
        # rewrite only rf2d_pred / has2d / pred_bin; every other column is kept byte-identical
        mc = pd.read_csv(MC_CSV, dtype=str, keep_default_na=False)
        new = out.set_index(out.exp_db_id.str.upper()).rf2d_mcc
        v = mc.exp_db_id.str.upper().map(new)
        old = pd.to_numeric(mc.rf2d_pred, errors="coerce")
        changed = int((~np.isclose(old, v, equal_nan=True)).sum())
        mc["rf2d_pred"] = [repr(float(x)) if pd.notna(x) else "" for x in v]
        mc["has2d"] = ["True" if pd.notna(x) else "False" for x in v]
        mc["pred_bin"] = pd.cut(v, [-1.01, .5, .8, .95, 1.01], labels=["<0.5", "0.5-0.8", "0.8-0.95", ">0.95"]).astype(object).where(lambda b: b.notna(), "")
        mc.to_csv(MC_CSV, index=False)
        print(f"updated rf2d_pred in {MC_CSV}: {changed} rows changed, {int(v.notna().sum())} non-null")


if __name__ == "__main__":
    main()
