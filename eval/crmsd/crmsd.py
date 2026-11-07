#!/usr/bin/env python3
"""Complex RMSD (cRMSD) as reported in the paper, for one eval shard.

cRMSD = PyMOL `align` of all polymer atoms (protein + RNA) of the prediction onto the reference, RMSD over all aligned
atoms BEFORE PyMOL's outlier rejection (`align(...)[3]`).
The metric DBs (eval/engine, AF3InterfaceEval a23a04f) store the RMSD AFTER the 5 default rejection cycles
(`align(...)[0]`) in Complex_RMSD[0]; both are written here, and rmsd_refined reproduces the DB value exactly.
Object preparation follows eval/engine/evaluation/extractEvalMetrics.py (Biopython mmCIF -> PDB, load as
referenceA/predictedA, reset=1).

  python eval/crmsd/crmsd.py <tag> <af3_output_dir> <seed> <sample> <ids.txt> <gt_dir> <out.tsv.gz>
"""
import glob, gzip, os, sys, tempfile
from Bio.PDB import MMCIFParser, PDBIO
from pymol import cmd

tag, exp, seed, sample, ids_file, gt_dir, out = sys.argv[1:8]


def to_pdb(path):
    f = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False); f.close()
    io = PDBIO(); io.set_structure(MMCIFParser(QUIET=True).get_structure("x", path)); io.save(f.name)
    return f.name


rows = []
for i in open(ids_file).read().split():
    cifs = sorted(glob.glob(f"{exp}/{i}_*/seed-{seed}_sample-{sample}/model.cif")
                  + glob.glob(f"{exp}/{i.upper()}_*/seed-{seed}_sample-{sample}/model.cif")
                  + glob.glob(f"{exp}/{i}/seed-{seed}_sample-{sample}/model.cif"))
    if not cifs:
        rows.append((i, "MISSING", "", "", "", "")); continue
    ref = tgt = None
    try:
        ref, tgt = to_pdb(f"{gt_dir}/{i}.cif"), to_pdb(cifs[0])
        cmd.reinitialize(); cmd.load(ref, "referenceA"); cmd.load(tgt, "predictedA")
        for r in ("MG", "K", "NA", "CA", "ZN", "CL", "SO4"):
            cmd.remove(f"resn {r}")
        a = cmd.align("predictedA and polymer", "referenceA and polymer", quiet=1, reset=1)
        rows.append((i, "ok", f"{a[0]:.6f}", f"{a[3]:.6f}", a[1], a[4]))
    except Exception as e:
        rows.append((i, f"ERROR:{type(e).__name__}:{e}".replace("\t", " ")[:200], "", "", "", ""))
    finally:
        for p in (ref, tgt):
            if p and os.path.exists(p): os.remove(p)

tmp = out + ".tmp"
with gzip.open(tmp, "wt") as f:
    f.write("id\tstatus\trmsd_refined\trmsd_all\tn_atoms_refined\tn_atoms_all\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
os.replace(tmp, out)
print(f"{tag} s{seed} m{sample}: {sum(r[1] == 'ok' for r in rows)}/{len(rows)} ok -> {out}")
