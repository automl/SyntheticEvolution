#!/usr/bin/env python3
"""Fig. 3b, top row: VARNA line (arc) diagrams of the reference topology and the topologies encoded in the SHS.

Same drawing as the submitted figure (RnaBench RnaVisualizer.visualize_rna, line layout): pairs are written to a CT file
(one partner per nucleotide, last pair wins) and drawn by VARNA with `-algorithm line -bpStyle lw -resolution 8.0`;
multiplets, non-canonical and lone pairs are passed as auxiliary pairs and drawn green (#00FF00). The helper functions
below are the RnaBench versions used for the figure (SPOT-RNA2-derived type_pairs / ct_file_output).

Sources (all in this repo): reference = DSSR pairs (1v1 complexes: pairs2 of data/gates/inputs/rp1v1_...csv;
single chains: data/construct/construct_pairs_91.json *_all); RNAformer = data/rnaformer/rnaformer_all_preds.pkl
(generator lookup); SPOT-RNA / RNAfold = data/topologies/<method>_<set>.json. Needs the `varna` command (bioconda
package varna 3.93, e.g. in the RnaBench env; set VARNA=/path/to/varna otherwise).

  python analysis/benchmark/fig3b_arcs.py                     # 1SJ4 (1v1 complex)
  python analysis/benchmark/fig3b_arcs.py --pdb 6WJS --set single
Output: results/analysis/fig3b/<PDB>_<source>_line.png (+ .ct), and the MCC of each topology vs the reference.
"""
import argparse, importlib.util, json, os, shutil, subprocess, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("fig3a", ROOT / "analysis" / "benchmark" / "fig3a_mcc_vs_3d.py")
F = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(F)
AUX_COLOR = "#00FF00"


# --- RnaBench helpers (as used for the submitted figure) ---------------------------------------------------------------
def lone_pair(pairs):
    lone = []
    pairs.sort()
    for p in pairs:
        if [p[0] - 1, p[1] + 1] not in pairs and [p[0] + 1, p[1] - 1] not in pairs:
            lone.append(p)
    return lone


def type_pairs(pairs, sequence):
    s = [c.upper() for c in sequence]; wc, gu, other = [], [], []
    for p in pairs:
        b = [s[p[0]], s[p[1]]]
        if b in (["A", "U"], ["U", "A"], ["G", "C"], ["C", "G"]): wc.append(p)
        elif b in (["G", "U"], ["U", "G"]): gu.append(p)
        else: other.append(p)
    return wc, gu, other


def get_multiplets(pairs):
    d = defaultdict(list)
    for p in pairs: d[p[0]].append(p[1]); d[p[1]].append(p[0])
    return list({tuple(sorted([k, q])) for k, v in d.items() if len(v) > 1 for q in v})


def ct_file_output(pairs, seq, pid, path, algorithm):
    n = len(seq); partner = np.zeros(n, dtype=int)
    for i, j in pairs: partner[i] = j + 1; partner[j] = i + 1
    rows = [f"{k + 1}\t\t{seq[k]}\t\t{k}\t\t{(k + 2) if k + 1 < n else 0}\t\t{partner[k]}\t\t{k + 1}" for k in range(n)]
    Path(path).write_text(f"{n}\t\t{pid}\t\t{algorithm} output\n\n" + "\n".join(rows) + "\n")


def draw(pairs, seq, pid, tag, outdir, varna):
    pairs = sorted([list(x) for x in {tuple(sorted(p)) for p in pairs}], key=lambda x: x[0])   # lists, as in RnaBench
    aux = get_multiplets(pairs) + [tuple(p) for p in type_pairs(pairs, seq)[2]] + [tuple(p) for p in lone_pair(pairs)]
    aux = sorted({tuple(x) for x in aux})
    aux_arg = ";".join(f"({i + 1},{j + 1}):color=\"\"{AUX_COLOR}\"\"" for i, j in aux)
    ct = outdir / f"{pid}_{tag}.ct"; ct_file_output(pairs, seq, pid, ct, tag)
    png = outdir / f"{pid}_{tag}_line.png"
    subprocess.run([varna, "-i", str(ct.resolve()), "-o", str(png.resolve()), "-algorithm", "line", "-resolution", "8.0",
                    "-bpStyle", "lw", "-auxBPs", aux_arg], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    return png, len(pairs), len(aux)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", default="1SJ4"); ap.add_argument("--set", choices=["complex1v1", "single"], default="complex1v1")
    a = ap.parse_args(); pid = a.pdb.lower(); idset = "complex" if a.set == "complex1v1" else "single"
    varna = os.environ.get("VARNA") or shutil.which("varna") or os.path.expanduser(
        f"~/.conda/envs/{os.environ.get('RNABENCH_ENV', 'RnaBench')}/bin/varna")
    out = ROOT / "results" / "analysis" / "fig3b"; out.mkdir(parents=True, exist_ok=True)
    if a.set == "complex1v1":
        d = pd.read_csv(ROOT / "data" / "gates" / "inputs" / "rp1v1_2D_full_with_pairs_af3_rnaformerN100_dssr.csv")
        ref = F.parse_pairs(d[(d.pdbid.str.lower() == pid) & (d.method == "rnaformer")].pairs2.values[0])
    else:
        ref = next(F.parse_pairs(r["pairs"]) for r in json.load(open(ROOT / "data" / "construct" / "construct_pairs_91.json"))
                   if r["id"] == f"{pid}_all")
    seq, _ = F.query_and_msa("rnaformer", idset, pid); L = len(seq)
    topo = {"reference": ref,
            "rnaformer": F.lookup(F.table_from(ROOT / "data" / "rnaformer" / "rnaformer_all_preds.pkl"), pid, seq)}
    for m in ("spotrna", "rnafold"):
        topo[m] = F.lookup(F.table_from(ROOT / "data" / "topologies" / f"{m}_{a.set}.json"), pid, seq)
    for tag, P in topo.items():
        png, n, naux = draw(P, seq, pid.upper(), tag, out, varna)
        print(f"  {tag:9} MCC vs reference {F.mcc(F.matrix(P, L), F.matrix(ref, L)):5.2f} | {n} pairs, {naux} drawn green -> {png}")


if __name__ == "__main__":
    main()
