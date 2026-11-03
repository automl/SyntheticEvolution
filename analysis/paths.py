"""Shared locations + helpers for all analysis/table scripts.

Everything is resolved relative to the repository or from environment variables set by config.sh:
  CRC_DATA      folder holding the unzipped download bundles (default: <repo>/downloads)
  CRC_INPUTS    2026_neurips_shs_inputs/   (AF3 input JSONs, reference mmCIFs)
  CRC_OUTPUTS   2026_neurips_shs_outputs/  (AF3 structures + summary confidences, metric DBs)
  DB_DIR        per-shard metric DBs        (default: <repo>/results/db, else the bundle's metric_dbs/benchmark)
Optional, only for the authors' original cluster layout: AF_OUTPUT / AF_INPUT (+ eval/conditions.tsv dir names).
"""
import csv, glob, importlib.util, os, re, sqlite3, statistics as st
from pathlib import Path

CRC = Path(os.environ.get("CRC_ROOT", Path(__file__).resolve().parents[1]))
DATA_ROOT = Path(os.environ.get("CRC_DATA", CRC / "downloads"))
INPUTS = Path(os.environ.get("CRC_INPUTS", DATA_ROOT / "2026_neurips_shs_inputs"))
OUTPUTS = Path(os.environ.get("CRC_OUTPUTS", DATA_ROOT / "2026_neurips_shs_outputs"))
AF_OUTPUT = os.environ.get("AF_OUTPUT", "")   # original layout (optional)
AF_INPUT = os.environ.get("AF_INPUT", "")

GT_SINGLE = os.environ.get("GT_SINGLE", str(INPUTS / "ground_truth" / "single_91"))
GT_COMPLEX = os.environ.get("GT_COMPLEX", str(INPUTS / "ground_truth" / "complex_377"))
GT_LONG = os.environ.get("GT_LONG", str(INPUTS / "ground_truth" / "long_25"))
GT_INTERVENTION = str(CRC / "data" / "gt_intervention")
RNABENCH_HOME = os.environ.get("RNABENCH_HOME", "")

def _db(own, shipped):  # own re-scored DBs if present, else the ones shipped in the outputs bundle (same rule as config.sh)
    o, s = CRC / "results" / own, OUTPUTS / "metric_dbs" / shipped
    return str(s if s.is_dir() and not o.is_dir() else o)
DB_DIR = os.environ.get("DB_DIR", _db("db", "benchmark"))
DB_LONG = os.environ.get("DB_LONG", _db("db_long", "long_rna"))
DB_RHOFOLD = os.environ.get("DB_RHOFOLD", _db("db_rhofold", "rhofold"))   # RhoFold+ runs: <tag>_single_s1_m0.db
OUT_ROOT = str(CRC / "results" / "analysis")
CONSTRUCT_DIR = str(CRC / "results" / "analysis" / "construction")   # construction-sweep outputs
CONSTRUCT_DATA = str(CRC / "data" / "construct")                      # construction-sweep inputs
USALIGN = os.environ.get("USALIGN", str(CRC / "eval" / "engine" / "evaluation" / "USalign" / "USalign"))
IDS_DIR = CRC / "data" / "ids"
MC_CSV = str(CRC / "data" / "gates" / "master_combined.csv")
MD_CSV = str(CRC / "data" / "gates" / "metadata_468.csv")
FAM_CSV = str(CRC / "data" / "gates" / "family_leakage.csv")
TIMING_CSV = str(CRC / "data" / "long_rna" / "pipeline_timing_comparison.csv")
os.makedirs(CONSTRUCT_DIR, exist_ok=True)

SEEDS = (1, 0, 42, 137)   # 4 AF3 seeds x 5 diffusion samples = 20 predictions per target

_spec = importlib.util.spec_from_file_location("neurips_ids", CRC / "data" / "splits" / "neurips_ids.py")
NEURIPS_IDS = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(NEURIPS_IDS)
CONDITIONS = {(r["tag"], r["idset"]): r for r in csv.DictReader(open(CRC / "eval" / "conditions.tsv"), delimiter="\t")}


def bundle_layout():
    """True if the unzipped outputs bundle is available (preferred over the original cluster layout)."""
    return (OUTPUTS / "af3").is_dir()


def prediction_root(tag, idset, seed):
    """Directory whose sub-directories are the AF3 target dirs holding seed-<seed>_sample-*/ for this condition."""
    if bundle_layout():
        base = OUTPUTS / "af3" / f"{tag}_{idset}"
        for sub in (["seed_1"] if seed == 1 else []) + ["seeds_0_42_137", "seeds_1_0_42_137"]:
            if (base / sub).is_dir():
                return str(base / sub)
        raise FileNotFoundError(base)
    if not AF_OUTPUT:
        raise FileNotFoundError("unzip 2026_neurips_shs_outputs.zip into CRC_DATA (or set AF_OUTPUT for the original layout)")
    row = CONDITIONS[(tag, idset)]
    return f"{AF_OUTPUT}/{row['seed1_dir'] if seed == 1 else row['multiseed_dir']}"


def prediction_files(tag, idset, fname="model.cif", seeds=SEEDS):
    """All files <target>/seed-<s>_sample-<k>/<fname> of a condition, each seed taken from exactly one source."""
    out = []
    for s in seeds:
        out += glob.glob(f"{prediction_root(tag, idset, s)}/*/seed-{s}_sample-*/{fname}")
    return out


def case_study_dir(name):
    """AF3 output dir of a case study (4P8Z/6WJS sweeps, construct91_multiseed); runs not in the bundle: AF_OUTPUT."""
    if (OUTPUTS / "af3_case_studies" / name).is_dir():
        return str(OUTPUTS / "af3_case_studies" / name)
    return f"{AF_OUTPUT}/{name}"


CONSTRUCT_N = int(os.environ.get("CONSTRUCT_N", "1000"))   # SHS depth of the 91-RNA construction run to analyse
CONSTRUCT_RUN = "construct91_multiseed" if CONSTRUCT_N == 1000 else f"construct91_N{CONSTRUCT_N}_multiseed"
CONSTRUCT_PREFIX = "construct91" if CONSTRUCT_N == 1000 else f"construct91_N{CONSTRUCT_N}"   # output file prefix


def input_jsons(idset):
    """{PDBID: path} of the natural-MSA AF3 input JSONs (sequences) for single|complex targets."""
    if (INPUTS / "af3_inputs").is_dir():
        files = glob.glob(str(INPUTS / "af3_inputs" / f"af3_{idset}" / "*" / "*.json"))
    else:
        files = glob.glob(f"{AF_INPUT}/nature_methods_full_evaluation_data_files/*_data.json")
    return {os.path.basename(p)[:4].upper(): p for p in files}


def split_ids(name):
    """Target ids (lower case) of a paper split. Accepts neurips_ids variable names
    ('rna_monomers_a2021_orphan_ids') or the old split-CSV names ('RNA_Monomers_a2021_alphafold_orphan')."""
    if not hasattr(NEURIPS_IDS, name):
        m = re.match(r"(All|RNA-Protein|RNA_Monomers)(?:_(a2021|b2021))?_[A-Za-z0-9]+?(?:_(orphan|non_orphan|1v1|multichain))?$", name)
        if not m:
            raise KeyError(name)
        base = {"All": "all_data", "RNA-Protein": "protein_rna", "RNA_Monomers": "rna_monomers"}[m.group(1)]
        name = "_".join(x for x in (base, m.group(2), m.group(3)) if x) + "_ids"
    return {i.lower() for i in getattr(NEURIPS_IDS, name)}


def idof(s):
    m = re.search(r"([0-9][A-Za-z0-9]{3})", s or "")
    return m.group(1).lower() if m else None


def db_values(tag, idset, col, seeds=SEEDS, db_dir=None):
    """{id: [values over seeds x samples]} from <tag>_<idset>_s<seed>_m<sample>.db (RNA_TM == 0 treated as missing)."""
    tbl = "pred_rna_rna" if idset in ("single", "long") else "pred_protein_rna"
    acc = {}
    for seed in seeds:
        for db in sorted(glob.glob(f"{db_dir or DB_DIR}/{tag}_{idset}_s{seed}_m*.db")):
            for fn, v in sqlite3.connect(f"file:{db}?mode=ro", uri=True).execute(f"SELECT FileName,{col} FROM {tbl}"):
                i = idof(fn)
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    continue
                if i and not (col == "RNA_TM" and v == 0):
                    acc.setdefault(i, []).append(v)
    return acc


def ms_mean(tag, idset, col, seeds=SEEDS):
    """Per-target mean over all seeds x samples."""
    return {i: st.mean(v) for i, v in db_values(tag, idset, col, seeds).items() if v}
