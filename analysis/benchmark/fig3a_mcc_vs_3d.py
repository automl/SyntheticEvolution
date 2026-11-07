#!/usr/bin/env python3
"""Fig. 3a: 3D accuracy of SHS-guided AF3 predictions binned by the 2D accuracy of the topology encoded in the SHS.

One panel per topology source (RNAformer, SPOT-RNA, RNAfold). x = MCC of the encoded topology vs the reference base pairs
(bins <0.5, 0.5-0.7, 0.7-0.9, >=0.9); y = per-target mean over 4 AF3 seeds x 5 samples (same loader as the paper tables,
tables/make_crc_tables.py:load_condition).

  --set complex1v1 (default): one-RNA-one-protein complexes, y = DockQ; reference pairs from
                              data/gates/inputs/rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv (pairs2, DSSR, all pairs)
  --set single:               single-chain RNAs, y = TM-score; reference pairs from data/construct/construct_pairs_91.json (*_all)

Topologies (exactly what shs/rna_msa_generator_base_pair.py encoded):
  RNAformer: data/rnaformer/rnaformer_all_preds.pkl, looked up as in the generator (PDB id, else sequence).
  SPOT-RNA / RNAfold: computed on the fly by the generator; supply them with --spotrna / --rnafold
                      (default data/topologies/{spotrna,rnafold}_<set>.json). Accepted formats:
                      .json  list of {"pdb_id"|"id"|"pdbid": .., "pairs": [[i, j], ..]} or {id: pairs}
                      .pkl / .csv  table with columns pdb_id (or id/pdbid) and pairs [, sequence]
                      0- or 1-based indices (detected from the base-pair complementarity of the target sequence).
Each supplied topology is checked against the shipped SHS of that condition: paired positions must be the conserved
columns of the alignment (stem-keep 0.99), reported as the mean Jaccard; --strict aborts below 0.9.

Outputs: results/analysis/fig3a/fig3a_<set>.{pdf,png,svg} (3 panels), fig3a_<set>_<source>.{pdf,png,svg} (single
panels, 600 dpi, tight crop shared by the panels of a set: same size and axes rectangle), fig3a_<set>_per_target.csv, fig3a_<set>_bins.csv
"""
import argparse, ast, glob, importlib.util, json, os, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import CRC, INPUTS, MC_CSV, OUT_ROOT  # noqa: E402
import numpy as np, pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_spec = importlib.util.spec_from_file_location("crc_tables", CRC / "tables" / "make_crc_tables.py")
_tables = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_tables)

SOURCES = [("rnaformer", "SHS (RNAformer)"), ("spotrna", "SHS (SPOT-RNA)"), ("rnafold", "SHS (RNAfold)")]
BINS = [(-np.inf, 0.5, "< 0.5"), (0.5, 0.7, "0.5 - 0.7"), (0.7, 0.9, "0.7 - 0.9"), (0.9, np.inf, ">= 0.9")]
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}
SETS = {"complex1v1": dict(idset="complex", metric="Dockq", ylabel="DockQ"),
        "single": dict(idset="single", metric="RNA_TM", ylabel="TM-score")}


def parse_pairs(p):
    if isinstance(p, str): p = ast.literal_eval(p)
    return [(int(x[0]), int(x[1])) for x in p if len(x) >= 2 and int(x[0]) != int(x[1])]


def matrix(pairs, L):
    m = np.zeros((L, L), bool)
    for i, j in pairs:
        if 0 <= i < L and 0 <= j < L: m[i, j] = m[j, i] = True
    return m


def mcc(pred, true):
    """MCC over the full symmetric L x L matrix (TN included), as in eval/selfcons."""
    tp = (pred & true).sum(); tn = (~pred & ~true).sum(); fp = (pred & ~true).sum(); fn = (~pred & true).sum()
    den = np.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    return float((tp * tn - fp * fn) / den) if den else 0.0


def query_and_msa(tag, idset, pid):
    """(query sequence, SHS rows without insertions) from the shipped AF3 input of condition <tag>."""
    js = glob.glob(str(INPUTS / "af3_inputs" / f"{tag}_{idset}" / pid.upper() / "*.json"))
    if not js: return None, []
    rna = [s["rna"] for s in json.load(open(js[0]))["sequences"] if "rna" in s][0]
    rows = [re.sub("[a-z]", "", l) for l in (rna.get("unpairedMsa") or "").splitlines() if l and not l.startswith(">")]
    return rna["sequence"].upper(), rows


def table_from(path):
    path = Path(path)
    if path.suffix == ".json":
        d = json.load(open(path))
        if isinstance(d, dict): d = [{"pdb_id": k, "pairs": v} for k, v in d.items()]
        df = pd.DataFrame(d)
    elif path.suffix == ".pkl": df = pd.read_pickle(path)
    else: df = pd.read_csv(path, sep=None, engine="python")
    idc = next(c for c in ("pdb_id", "id", "pdbid", "exp_db_id") if c in df.columns)
    df = df.rename(columns={idc: "pdb_id"})
    if "sequence" not in df.columns: df["sequence"] = ""
    return df[["pdb_id", "pairs", "sequence"]]


def lookup(df, pid, seq):
    """Generator lookup: PDB id (first four characters, case-insensitive), else exact sequence."""
    q = df[df.pdb_id.astype(str).str.lower() == pid.lower()]
    if q.empty: q = df[df.sequence == seq]
    return parse_pairs(q.pairs.values[0]) if len(q) else None


def index_offset(topos, seqs):
    """0 or 1: the indexing under which the supplied pairs are more often complementary in the target sequence."""
    score = {}
    for off in (0, 1):
        ok = n = 0
        for pid, P in topos.items():
            s = seqs[pid]
            for i, j in P:
                i, j = i - off, j - off
                if 0 <= i < len(s) and 0 <= j < len(s): ok += (s[i] + s[j]) in CANON; n += 1
        score[off] = ok / max(n, 1)
    off = max(score, key=score.get)
    print(f"  index base: {off}-based (complementary pairs: 0-based {score[0]:.2f}, 1-based {score[1]:.2f})")
    return off


def encoded_jaccard(P, seq, rows):
    """Jaccard between the paired positions of P and the conserved (>90% identical among non-gap rows) columns of the shipped SHS."""
    rows = [r for r in rows[1:] if len(r) == len(seq)]
    if not rows: return np.nan
    A = np.array([list(r) for r in rows]); ungapped = A != "-"            # deletions are not substitutions
    ident = ((A == np.array(list(seq))) & ungapped).sum(0) / np.maximum(ungapped.sum(0), 1)
    cons = set(np.flatnonzero(ident > 0.9)); paired = {x for p in P for x in p}
    return len(cons & paired) / len(cons | paired) if cons | paired else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=SETS, default="complex1v1")
    ap.add_argument("--spotrna"); ap.add_argument("--rnafold")
    ap.add_argument("--min-preds", type=int, default=20)
    ap.add_argument("--strict", action="store_true", help="abort if a supplied topology does not match the shipped SHS")
    a = ap.parse_args(); cfg = SETS[a.set]; idset = cfg["idset"]
    out = Path(OUT_ROOT) / "fig3a"; out.mkdir(parents=True, exist_ok=True)

    # targets + reference pairs
    if a.set == "complex1v1":
        ids = set(pd.read_csv(MC_CSV).query("group == '1v1'").exp_db_id.str.lower())
        d = pd.read_csv(CRC / "data" / "gates" / "inputs" / "rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv")
        ref = {r.pdbid.lower(): parse_pairs(r.pairs2) for r in d[d.method == "rnaformer"].itertuples() if r.pdbid.lower() in ids}
    else:
        ref = {r["pdb_id"].lower(): parse_pairs(r["pairs"]) for r in json.load(open(CRC / "data" / "construct" / "construct_pairs_91.json"))
               if r["id"].endswith("_all")}
    print(f"[{a.set}] targets with reference pairs: {len(ref)}")

    sources = {"rnaformer": CRC / "data" / "rnaformer" / "rnaformer_all_preds.pkl",
               "spotrna": Path(a.spotrna) if a.spotrna else CRC / "data" / "topologies" / f"spotrna_{a.set}.json",
               "rnafold": Path(a.rnafold) if a.rnafold else CRC / "data" / "topologies" / f"rnafold_{a.set}.json"}
    rows, panels = [], []
    for tag, label in SOURCES:
        if not sources[tag].is_file():
            print(f"  {label}: topology file {sources[tag]} not found -> panel left empty"); panels.append((tag, label, None)); continue
        print(f"  {label}: {sources[tag]}")
        df = table_from(sources[tag]); y3d, _ = _tables.load_condition(tag, idset, a.min_preds)
        seqs, msas, topos = {}, {}, {}
        for pid in sorted(ref):
            seq, msa = query_and_msa(tag, idset, pid)
            if seq is None: continue
            P = lookup(df, pid, seq)
            if P is None: continue
            seqs[pid], msas[pid], topos[pid] = seq, msa, P
        off = index_offset(topos, seqs) if tag != "rnaformer" else 0
        jac = []
        for pid, P in topos.items():
            P = [(i - off, j - off) for i, j in P]; L = len(seqs[pid])
            j = encoded_jaccard(P, seqs[pid], msas[pid]); jac.append(j)
            v = y3d.get(pid, {}).get(cfg["metric"])
            rows.append(dict(source=tag, id=pid, length=L, n_pairs=len(P), mcc=mcc(matrix(P, L), matrix(ref[pid], L)),
                             encoded_jaccard=j, y=v))
        jm = np.nanmean(jac); low = int(np.sum(np.array(jac) < 0.9))
        print(f"    matched {len(topos)}/{len(ref)} targets; SHS-encoding check: mean Jaccard {jm:.3f}, {low} targets < 0.9")
        if a.strict and jm < 0.9: sys.exit(f"{label}: supplied topologies do not match the shipped SHS (Jaccard {jm:.3f})")
        panels.append((tag, label, True))

    t = pd.DataFrame(rows).dropna(subset=["y"])
    t.to_csv(out / f"fig3a_{a.set}_per_target.csv", index=False)
    summ = []
    for tag, label, ok in panels:
        g = t[t.source == tag]
        if not ok or g.empty: continue
        rho = spearmanr(g.mcc, g.y)
        print(f"  {label}: n={len(g)} mean MCC {g.mcc.mean():.3f}; Spearman(MCC, {cfg['ylabel']}) = {rho.correlation:.2f} (p={rho.pvalue:.1e})")
        for lo, hi, lab in BINS:
            b = g[(g.mcc >= lo) & (g.mcc < hi)]
            summ.append(dict(source=tag, bin=lab, n=len(b), median=b.y.median(), mean=b.y.mean(), rho_all=rho.correlation))
    pd.DataFrame(summ).to_csv(out / f"fig3a_{a.set}_bins.csv", index=False)

    # plot (style of the submitted Fig. 3a: light-blue boxes, red median, serif font, n per bin; no top/right spines)
    plt.rcParams.update({"font.family": "serif", "font.size": 11})

    def draw(ax, tag, label, ok):
        ax.set_title(label, fontsize=14, loc="left"); ax.set_xlabel("MCC", fontsize=13); ax.set_ylabel(cfg["ylabel"], fontsize=13)
        ax.set_ylim(-0.05, 1.03); ax.set_yticks(np.arange(0, 1.01, 0.2)); ax.yaxis.set_major_formatter("{x:.1f}")
        ax.grid(True, color="#e6e6e6", linewidth=0.8); ax.set_axisbelow(True)
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        g = t[t.source == tag]
        data = [g[(g.mcc >= lo) & (g.mcc < hi)].y.values for lo, hi, _ in BINS]
        ax.boxplot([x if len(x) else [np.nan] for x in data], widths=0.75, patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor="#a7d8ee", edgecolor="#7f7f7f"), medianprops=dict(color="#d62728", linewidth=1.2),
                   whiskerprops=dict(color="#7f7f7f"), capprops=dict(color="#7f7f7f"))
        ax.set_xlim(0.4, 4.6); ax.set_xticks(range(1, 5))
        ax.set_xticklabels([f"{lab}\n(n = {len(x)})" for (_, _, lab), x in zip(BINS, data)], fontsize=10)
        if not ok: ax.text(2.5, 0.45, "topologies\nnot supplied", ha="center", va="center", color="#999999")

    # combined overview (3 in a row)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 2.9))
    for ax, (tag, label, ok) in zip(axes, panels): draw(ax, tag, label, ok)
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"): fig.savefig(out / f"fig3a_{a.set}.{ext}", bbox_inches="tight", dpi=600)
    plt.close(fig)
    # individual panels: identical canvas and axes rectangle, cropped to the union of the tight bounding boxes of all
    # panels of the set, so every panel is tight and has the same size and axes position (they line up side by side)
    figs = []
    for tag, label, ok in panels:
        fig = plt.figure(figsize=(3.6, 2.9))
        ax = fig.add_axes([0.20, 0.30, 0.76, 0.60])          # [left, bottom, width, height] in figure fractions
        draw(ax, tag, label, ok); figs.append((tag, fig))
    from matplotlib.transforms import Bbox
    boxes = []
    for _, fig in figs:
        r = fig.canvas.get_renderer(); boxes.append(fig.get_tightbbox(r))
    crop = Bbox.union(boxes).padded(0.02)                     # inches
    for tag, fig in figs:
        for ext in ("pdf", "png", "svg"): fig.savefig(out / f"fig3a_{a.set}_{tag}.{ext}", dpi=600, bbox_inches=crop)
        plt.close(fig)
    print(f"  -> {out}/fig3a_{a.set}.pdf and fig3a_{a.set}_<source>.pdf")

if __name__ == "__main__":
    main()
