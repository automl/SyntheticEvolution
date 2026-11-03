#!/usr/bin/env python3
"""Compute the TRUE RNAformer 2D-prediction accuracy vs ground truth (gt_dssr), per eval id.

Correct gating variable for "SHS quality depends on how well RNAformer predicts the 2D".
NOT self_consistency (seed-vs-realised-3D) and NOT _dssr (realised-2D-vs-GT): this is the
RNAformer SEED prediction scored against gt_dssr.

Prediction pairs (the generator's own lookup, rna_msa_generator_base_pair.py:223-225):
    df = read_pickle('...rnaformer2025_pdb_id2pairs_mapping.pkl'); df[df.sequence == seq].pairs
  keyed by the EXACT eval sequence (seq = sequence.upper()). So every eval id has a prediction.
Ground-truth pairs: pairs2 (reference_method == gt_dssr) from the SS archive (monomers) and the
  rp1v1 with-pairs csv (1v1). The GT (dssr) sequence can be +-1 nt vs the modelled sequence
  (a known PDB artefact), so the GT is ALIGNED into the eval-sequence frame by substring before
  scoring (indices shifted; pairs falling outside the eval sequence are dropped).
Output: reviewer_analysis/rnaformer_2d_vs_gt.csv  [exp_db_id, group, rf2d_precision/recall/f1/mcc].
"""
import pickle, ast, difflib
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
NE = HERE.parent
PKL = NE / "2026_neurips_code" / ("synthetic_msa_selection_max_rna_len_200_min_rna_len_15_"
                                  "max_protein_len_200_rnaformer2025_pdb_id2pairs_mapping.pkl")
SS = NE / "ss_results" / "secondary_structure_neurips_2026_rna_rna" / "all_secondary_structure_data_rna_rna.csv"
MONO = NE / "RNA_Monomer_2D_self_consistency.csv"
RP = NE / "rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv"


def up(s):
    return str(s).upper()


def norm(pairs, n):
    out = set()
    for p in pairs:
        try:
            i, j = int(p[0]), int(p[1])
        except Exception:
            continue
        if i == j or i < 0 or j < 0 or i >= n or j >= n:
            continue
        out.add((min(i, j), max(i, j)))
    return out


def align_gt(pred_seq, gt_seq, gt_pairs):
    """Return GT pairs expressed in the pred-sequence index frame.

    The GT (dssr) sequence can differ from the modelled sequence by a +-1 nt PDB artefact or a
    single substitution/indel. We map GT positions -> pred positions via the longest-matching
    blocks (difflib), so any small edit is tolerated; GT pairs whose endpoints fall outside a
    matched block are dropped (they cannot be scored in the eval frame)."""
    if not gt_pairs:
        return []
    gi2pi = {}
    for blk in difflib.SequenceMatcher(None, gt_seq, pred_seq, autojunk=False).get_matching_blocks():
        for k in range(blk.size):
            gi2pi[blk.a + k] = blk.b + k
    out = []
    for p in gt_pairs:
        try:
            i, j = int(p[0]), int(p[1])
        except Exception:
            continue
        if i in gi2pi and j in gi2pi:
            out.append((gi2pi[i], gi2pi[j]))
    return out


def metrics(pred, gt, n):
    tp = len(pred & gt); fp = len(pred - gt); fn = len(gt - pred)
    tn = n * (n - 1) // 2 - tp - fp - fn
    # degenerate case: RNAformer predicted NO pairs (tp+fp==0). That is a total 2D failure vs a GT
    # that HAS pairs, so score it as 0 (worst) rather than leaving precision/mcc undefined -- this
    # keeps such ids (e.g. 7MKT) in the low-2D-quality bin instead of dropping them.
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else np.nan          # GT has no pairs -> recall undefined
    f1 = 2 * prec * rec / (prec + rec) if (prec and rec and not np.isnan(rec)) else 0.0
    den = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / den if den else 0.0        # degenerate 0/0 -> 0 (paper convention)
    return prec, rec, f1, mcc


def pred_in_frame(es, exact, items):
    """RNAformer predicted pairs mapped into the eval-sequence (es) index frame.
    Exact sequence hit first (monomers, by construction); else align the closest pkl sequence
    (rp1v1 can differ by the same +-1 nt PDB artefact) and remap pair indices. None if no match."""
    if es in exact:
        return norm(exact[es], len(es))
    best, br = None, 0.0
    for s, pairs in items:
        if abs(len(s) - len(es)) > 3:
            continue
        r = difflib.SequenceMatcher(None, s, es, autojunk=False).ratio()
        if r > br:
            br, best = r, (s, pairs)
    if best is None or br < 0.9:
        return None
    s, pairs = best
    m = {}
    for blk in difflib.SequenceMatcher(None, s, es, autojunk=False).get_matching_blocks():
        for k in range(blk.size):
            m[blk.a + k] = blk.b + k
    out = set()
    for p in pairs:
        i, j = int(p[0]), int(p[1])
        if i in m and j in m:
            out.add((min(m[i], m[j]), max(m[i], m[j])))
    return out


def main():
    pk = pickle.load(open(PKL, "rb"))
    pk["us"] = pk.sequence.map(up)
    exact = {s: g.pairs.iloc[0] for s, g in pk.groupby("us")}
    items = list(exact.items())

    # ---- ground-truth records: (exp_db_id, group, gt_seq, gt_pairs) ----
    recs = []
    ss = pd.read_csv(SS).dropna(subset=["pairs2", "sequence"])
    ss["base"] = ss.pdbid.astype(str).str.slice(0, 4).str.upper()
    gt_ss = {b: (up(g.sequence.iloc[0]), ast.literal_eval(g.pairs2.iloc[0])) for b, g in ss.groupby("base")}
    # eval sequence per monomer (matches the pkl exactly, by construction)
    mono = pd.read_csv(MONO)
    mono["base"] = mono.pdbid.str.upper().str.slice(0, 4)
    eval_seq = {r.base: up(r.sequence) for _, r in mono.drop_duplicates("base").iterrows()}
    for b, es in eval_seq.items():
        if b in gt_ss:
            recs.append((b, "rna-only", es, gt_ss[b]))
    if RP.exists():
        rp = pd.read_csv(RP).dropna(subset=["pairs2", "sequence"])
        for _, r in rp.drop_duplicates("pdbid").iterrows():
            s = up(r.sequence)
            recs.append((str(r.pdbid).upper(), "1v1", s, (s, ast.literal_eval(r.pairs2))))

    # PAPER CONVENTION: keep degenerate cases (GT has 0 base pairs, or RNAformer predicted 0 pairs)
    # scored as mcc=0 -> they land in the <0.5 (untrustworthy) bin instead of being dropped.
    rows, no_pred, degen = [], [], 0
    for eid, grp, es, gt in recs:
        pred = pred_in_frame(es, exact, items)
        if pred is None:
            no_pred.append(eid); continue
        n = len(es)
        gt_seq, gt_pairs = gt
        gt_set = norm(align_gt(es, gt_seq, gt_pairs), n)
        is_degen = (not gt_set) or (not pred)   # GT has 0 pairs, or RNAformer predicted 0 pairs
        degen += int(is_degen)
        prec, rec, f1, mcc = metrics(pred, gt_set, n)
        rows.append({"exp_db_id": eid, "group": grp, "n_res": n,
                     "rf2d_precision": prec, "rf2d_recall": rec, "rf2d_f1": f1, "rf2d_mcc": mcc,
                     "rf2d_degenerate": is_degen, "n_gt_pairs": len(gt_set), "n_pred_pairs": len(pred)})
    out = pd.DataFrame(rows).drop_duplicates("exp_db_id")
    out.to_csv(HERE / "rnaformer_2d_vs_gt.csv", index=False)
    print(f"matched {len(out)} ids | no-pred {len(no_pred)} {no_pred} | degenerate(mcc=0) {degen}")
    for g, d in out.groupby("group"):
        print(f"  {g:9s} n={len(d):3d} | f1 {d.rf2d_f1.mean():.3f}  mcc {d.rf2d_mcc.mean():.3f}  "
              f"prec {d.rf2d_precision.mean():.3f}  rec {d.rf2d_recall.mean():.3f}")


if __name__ == "__main__":
    main()
