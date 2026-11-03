#!/usr/bin/env python3
"""run_eval.py -- drop-in evaluation wrapper for AF3InterfaceEval.

Clone -> install -> point this at a folder of AlphaFold3 predictions (and, optionally, a
folder of deposited ground-truth mmCIFs) and it produces a self-contained SQLite DB that
feeds the regression/plotting pipeline. It orchestrates the existing sub-scripts the exact
way ``evaluation/mainExecution.py`` does (config.json-driven, PYTHONPATH=<repo>, OST steps
via ``conda run -n ost``), but adds:

  * arbitrary input filenames via a ``--pair-map`` CSV or a filename-stem convention,
  * per-run isolation (a fresh DB + a temporarily-swapped config.json that is restored),
  * tool/GT capability detection with skip-and-log (``--strict`` turns skips into errors),
  * quick / basic / full metric presets (or an explicit ``--metrics`` list),
  * era tagging (``--cutoff-date``) from the GT deposition date.

Nothing here re-implements a metric: it stages inputs into the layout the sub-scripts expect
(``<workdir>/data/<classification>/{af,pdb,dssr}``) and drives them. Anything whose input or
external tool is absent is skipped and logged rather than treated as an error.

Run it from within the eval conda env, e.g.::

    conda run -n rnabench python run_eval.py --pred ./preds --gt ./gt --out ./eval_run \
        --metrics basic

See the README "Inputs" section for the pair-map / sidecar file formats.
"""

import argparse
import csv
import glob
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DATABASE_DIR = REPO_ROOT / "database"
EVAL_DIR = REPO_ROOT / "evaluation"
UTILS_DIR = REPO_ROOT / "utils"
# env-overridable so concurrent shards isolate their config (see startConfig.py)
CONFIG_PATH = DATABASE_DIR / os.environ.get("EVAL_CONFIG_NAME", "config.json")

CLASSIFICATIONS = ("protein_rna", "rna_rna")

# Metric preset -> ordering rank. A step belongs to a run if its introduced-at rank is <=
# the selected preset rank (quick < basic < full). Metadata steps are opt-in flags on top.
PRESET_RANK = {"quick": 1, "basic": 2, "full": 3}


# --------------------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------------------
def log(msg):
    print(f"[run_eval] {msg}", flush=True)


def warn(msg):
    print(f"[run_eval] WARNING: {msg}", flush=True)


def skip(step, reason):
    print(f"[run_eval] SKIP {step}: needs {reason}", flush=True)


def die(msg, code=1):
    print(f"[run_eval] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


# --------------------------------------------------------------------------------------
# id / classification / date helpers
# --------------------------------------------------------------------------------------
def sanitize_id(raw):
    """Reduce an arbitrary id to a single underscore-free alphanumeric token.

    ``parseAFcif2DB``/``af3Metrics_new`` recover the id as ``filename.split('_')[1]``, so the
    staged id must contain no ``_``. We strip every non-alphanumeric char; the DB stores it
    upper-cased (both parsers ``.upper()`` it), so the join key is the upper form.
    """
    token = re.sub(r"[^A-Za-z0-9]", "", str(raw))
    return token


# AF3 output naming: (fold_)<id>(_s<seed>)_{model,summary_confidences,full_data,data}(_<k>).
# af3_core strips that scaffolding to the shared job core, so a prediction
# (fold_<id>_s1_model_0.cif) and its confidence/data JSONs collapse to the same key -- this is
# how sidecars are grouped to their prediction (exact, never fuzzy).
_AF3_AFFIX_RE = re.compile(
    r"^(?:fold_)?(.+?)(?:_s\d+)?"
    r"_(?:model|summary_confidences|full_data|confidences|data)(?:_\d+)?$")

_SEPS = set("_-. ")

# A classic 4-char PDB id (digit-first, not embedded in a longer alnum run).
_PDB_RE = re.compile(r"(?<![A-Za-z0-9])([0-9][A-Za-z0-9]{3})(?![A-Za-z0-9])")


def af3_core(stem):
    """Job core of a filename stem: strip the AF3 fold_/_s<seed>_<kind>_<k> scaffolding."""
    m = _AF3_AFFIX_RE.match(stem)
    if m:
        return m.group(1)
    return stem[len("fold_"):] if stem.startswith("fold_") else stem


def extract_pdb(stem):
    """Return the (upper-cased) classic PDB id token in a stem, or None."""
    m = _PDB_RE.search(stem or "")
    return m.group(1).upper() if m else None


def _sep_trim_prefix(pre):
    """Keep a common prefix only up to and including its last separator (token boundary)."""
    idx = max((i for i, ch in enumerate(pre) if ch in _SEPS), default=-1)
    return pre[:idx + 1] if idx >= 0 else ""


def _sep_trim_suffix(suf):
    """Keep a common suffix only from its first separator on (token boundary)."""
    idx = next((i for i, ch in enumerate(suf) if ch in _SEPS), len(suf))
    return suf[idx:] if idx < len(suf) else ""


def strip_common_affix(cores):
    """Learn one folder's naming template: strip the longest prefix/suffix shared by *all*
    stems, but only at separator boundaries so ids that merely share leading characters
    (1abc / 1abd) are never over-stripped. Returns (new_cores, (prefix, suffix))."""
    valid = [c for c in cores if c]
    if len(valid) < 2:
        return list(cores), ("", "")
    pre = _sep_trim_prefix(os.path.commonprefix(valid))
    suf = _sep_trim_suffix(os.path.commonprefix([c[::-1] for c in valid])[::-1])
    out = []
    for c in cores:
        if not c:
            out.append(c)
            continue
        s = c
        if pre and s.startswith(pre):
            s = s[len(pre):]
        if suf and s.endswith(suf):
            s = s[:len(s) - len(suf)]
        out.append(s or c)  # never collapse to empty
    return out, (pre, suf)


def stem_cores(paths, id_regex=None):
    """Map each path to its match core (upper-cased): an --id-regex capture, else af3_core
    followed by learned common-affix stripping across the set."""
    stems = {p: Path(p).stem for p in paths}
    if id_regex:
        out = {}
        for p, s in stems.items():
            m = re.search(id_regex, s)
            out[p] = ((m.group(1) if m.groups() else m.group(0)).upper() if m else None)
        return out
    raw = [af3_core(stems[p]) for p in paths]
    stripped, _ = strip_common_affix(raw)
    return {p: (c.upper() if c else None) for p, c in zip(paths, stripped)}


def _unique_index(kv):
    """Build key -> path, dropping keys that map to more than one distinct path (ambiguous)."""
    seen = {}
    ambiguous = set()
    for key, path in kv:
        if key is None:
            continue
        if key in seen and str(seen[key]) != str(path):
            ambiguous.add(key)
        else:
            seen[key] = path
    return {k: v for k, v in seen.items() if k not in ambiguous}


def detect_classification(cif_path):
    """protein_rna vs rna_rna from a predicted CIF, dependency-free.

    Reads ``_entity_poly.type`` tokens straight from the text so the wrapper needs no Bio
    import. protein_rna if a polypeptide entity co-occurs with RNA; rna_rna if only RNA.
    """
    try:
        text = Path(cif_path).read_text(errors="ignore")
    except OSError as exc:
        warn(f"cannot read {cif_path} for classification: {exc}")
        return None
    has_protein = "polypeptide(L)" in text
    has_rna = "polyribonucleotide" in text
    if has_protein and has_rna:
        return "protein_rna"
    if has_rna and not has_protein:
        return "rna_rna"
    if has_protein and not has_rna:
        warn(f"{Path(cif_path).name}: protein without RNA -- not an RNA-interface case.")
        return None
    return None


_DATE_RE = re.compile(r"(\d{4}-\d{2}-\d{2})")


def read_deposition_date(cif_path):
    """Best-effort deposition date (YYYY-MM-DD) for era tagging.

    Prefers ``_pdbx_database_status.recvd_initial_deposition_date``; falls back to the first
    ``_pdbx_audit_revision_history.revision_date`` (what parsePDBCIF2DB records as the
    deposited date). Returns ``None`` when neither is parseable.
    """
    try:
        lines = Path(cif_path).read_text(errors="ignore").splitlines()
    except OSError:
        return None
    for key in ("_pdbx_database_status.recvd_initial_deposition_date",
                "_pdbx_audit_revision_history.revision_date"):
        for i, line in enumerate(lines):
            if key in line:
                # key/value on the same line (single datablock form)
                m = _DATE_RE.search(line)
                if m:
                    return m.group(1)
                # otherwise the value is in the loop rows just below
                for follow in lines[i + 1:i + 8]:
                    m = _DATE_RE.search(follow)
                    if m:
                        return m.group(1)
    return None


def rewrite_entry_id(cif_path, new_id):
    """Force ``_entry.id`` to ``new_id`` so exp.PDBId matches the AF-side exp_db_id join key."""
    path = Path(cif_path)
    try:
        lines = path.read_text(errors="ignore").splitlines(keepends=True)
    except OSError as exc:
        warn(f"cannot rewrite _entry.id in {cif_path}: {exc}")
        return
    out = []
    replaced = False
    data_idx = None
    for i, line in enumerate(lines):
        if line.lstrip().startswith("_entry.id") and not replaced:
            eol = "\n" if line.endswith("\n") else ""
            out.append(f"_entry.id   {new_id}{eol}")
            replaced = True
        else:
            out.append(line)
            if data_idx is None and line.startswith("data_"):
                data_idx = len(out)
    if not replaced:
        insert_at = data_idx if data_idx is not None else 0
        out.insert(insert_at, f"_entry.id   {new_id}\n")
    path.write_text("".join(out))


# --------------------------------------------------------------------------------------
# tool / capability detection
# --------------------------------------------------------------------------------------
def conda_has_env(env_name):
    conda = shutil.which("conda")
    if not conda:
        return False
    try:
        res = subprocess.run([conda, "env", "list"], capture_output=True, text=True, check=False)
    except OSError:
        return False
    for line in res.stdout.splitlines():
        first = line.split()
        if first and first[0] == env_name:
            return True
    return False


def module_importable(name):
    try:
        __import__(name)
        return True
    except Exception:
        return False


def detect_capabilities(args, has_gt, dssr_dir):
    caps = {}
    usalign = EVAL_DIR / "USalign" / ("USalign.exe" if os.name == "nt" else "USalign")
    caps["usalign"] = usalign.exists()
    caps["ost"] = conda_has_env(args.ost_env)
    # The eval scripts do `from pymol import cmd`, so they need the importable Python module in
    # the ambient env -- a pymol launcher on PATH is not enough (don't fall back to which()).
    caps["pymol"] = module_importable("pymol")
    caps["freesasa"] = module_importable("freesasa")
    caps["gt"] = has_gt
    caps["dssr"] = dssr_dir is not None and dssr_dir.is_dir() and any(dssr_dir.iterdir())
    token = (args.token or os.environ.get("SWISSMODEL_OST_TOKEN", "")).strip()
    caps["token"] = bool(token)
    caps["_token_value"] = token
    return caps


# --------------------------------------------------------------------------------------
# pairing
# --------------------------------------------------------------------------------------
class Pair:
    __slots__ = ("pred", "raw_core", "match_key", "id_token", "join_id",
                 "gt", "gt_forced", "data_json", "confidence_files",
                 "motif", "domain", "era", "staged_confidence", "staged_data_json")

    def __init__(self, pred):
        self.pred = pred
        self.raw_core = af3_core(Path(pred).stem)  # for exact sidecar grouping
        self.match_key = None                      # id used for the GT/DB join
        self.id_token = None                       # sanitised (underscore-free) form
        self.join_id = None                        # id_token upper-cased
        self.gt = None
        self.gt_forced = False                     # set by pair-map / interactive overlay
        self.data_json = None                      # resolved *_data.json path
        self.confidence_files = []                 # resolved confidence JSON paths
        self.motif = None
        self.domain = None
        self.era = None
        self.staged_confidence = 0
        self.staged_data_json = False

    def set_key(self, key):
        self.match_key = key
        self.id_token = sanitize_id(key)
        self.join_id = self.id_token.upper()


def _resolve(base_dir, value):
    """Resolve a pair-map cell to an absolute path (relative to the CSV's dir if needed)."""
    if not value:
        return None
    p = Path(value)
    if not p.is_absolute():
        p = (base_dir / p)
    return p if p.exists() else None


# ---- sidecar resolution (exact job-key first, id-containment only as a last resort) --------
# AF3 confidence sidecars come in two namings; both stage to the fold_<ID>_s1_<kind>_<k>.json
# form af3Metrics_new.py globs:
#   old  : *_summary_confidences_[0-4].json  +  *_full_data_[0-4].json
#   new  : *_summary_confidences.json        +  *_confidences.json   (plain = the full-data file)
_CONF_OLD_RE = re.compile(r"_(summary_confidences|full_data)_([0-4])\.json$")
_CONF_NEW_SUMMARY_RE = re.compile(r"summary_confidences\.json$")


def classify_confidence(name):
    """Map a confidence filename to (kind, index) where kind is summary_confidences|full_data,
    or None if it is not a confidence file. New-style plain `*confidences.json` -> full_data 0."""
    m = _CONF_OLD_RE.search(name)
    if m:
        return m.group(1), int(m.group(2))
    if _CONF_NEW_SUMMARY_RE.search(name):
        return "summary_confidences", 0
    if name.endswith("confidences.json"):     # plain confidences.json (not summary_)
        return "full_data", 0
    return None


def find_confidence_files(source, raw_core, id_token):
    """Confidence JSONs for one prediction: prefer an exact af3_core (job-key) match, fall back
    to id containment. ``source`` may be a folder or a single file. New-style files that live in
    per-seed subfolders are ignored in favour of the top-level (best-ranked) ones."""
    if source is None:
        return []
    source = Path(source)
    if source.is_file():
        return [source] if classify_confidence(source.name) else []
    if not source.is_dir():
        return []
    cands = [c for c in sorted(source.glob("*.json")) if classify_confidence(c.name)]
    exact = [c for c in cands if af3_core(c.stem).upper() == (raw_core or "").upper()]
    if exact:
        return exact
    return [c for c in cands if id_token and id_token in sanitize_id(c.stem).upper()]


def find_data_json(source, raw_core, id_token):
    """The AF3 *_data.json (MSA/templates) for one prediction, exact job-key first."""
    if source is None:
        return None
    source = Path(source)
    if source.is_file():
        return source if source.name.endswith("_data.json") else None
    if not source.is_dir():
        return None
    cands = [c for c in sorted(source.glob("*_data.json")) if not classify_confidence(c.name)]
    for c in cands:
        if af3_core(c.stem).upper() == (raw_core or "").upper():
            return c
    for c in cands:
        if id_token and id_token in sanitize_id(c.stem).upper():
            return c
    return None


def resolve_sidecars(pr, conf_source, data_source):
    pr.confidence_files = find_confidence_files(conf_source, pr.raw_core, pr.join_id)
    pr.data_json = find_data_json(data_source, pr.raw_core, pr.join_id)


# ---- prediction <-> ground-truth resolution ladder ----------------------------------------
def resolve_gt(preds, gt_paths, id_regex):
    """Assign each prediction its GT match and its DB id, via: PDB-id token -> learned core.
    Ambiguous keys (one key, several GTs) resolve to no match. pair-map/interactive-forced
    predictions are left untouched. Predictions with no confident match keep gt=None."""
    pred_core = stem_cores([p.pred for p in preds], id_regex)
    gt_core = stem_cores(gt_paths, id_regex)
    gt_by_core = _unique_index((gt_core[g], g) for g in gt_paths)
    # PDB-token rung is a heuristic; an explicit --id-regex overrides it entirely.
    use_pdb = id_regex is None
    gt_by_pdb = (_unique_index((extract_pdb(Path(g).stem), g) for g in gt_paths)
                 if use_pdb else {})

    for pr in preds:
        if pr.gt_forced:
            if pr.match_key is None:
                pr.set_key(pred_core.get(pr.pred) or af3_core(Path(pr.pred).stem).upper())
            continue
        ppdb = (extract_pdb(Path(pr.pred).stem) or extract_pdb(pr.raw_core)) if use_pdb else None
        pcore = pred_core.get(pr.pred)
        gt, key = None, None
        if ppdb and ppdb in gt_by_pdb:
            gt, key = gt_by_pdb[ppdb], ppdb
        elif pcore and pcore in gt_by_core:
            gt, key = gt_by_core[pcore], pcore
        pr.gt = gt
        pr.set_key(key or pcore or ppdb or af3_core(Path(pr.pred).stem).upper())


def _preview_unmatched(unmatched, preds, gt_paths, id_regex, limit=12):
    used = {str(p.gt) for p in preds if p.gt}
    gt_free = [g for g in gt_paths if str(g) not in used]
    pc = stem_cores([p.pred for p in unmatched], id_regex)
    gc = stem_cores(gt_free, id_regex)
    print(f"[run_eval] {len(unmatched)} prediction(s) with no ground-truth match:")
    for p in unmatched[:limit]:
        print(f"    pred  core={pc.get(p.pred) or '?':<16} {Path(p.pred).name}")
    if len(unmatched) > limit:
        print(f"    ... (+{len(unmatched) - limit} more)")
    print(f"[run_eval] {len(gt_free)} unmatched ground-truth file(s):")
    for g in gt_free[:limit]:
        print(f"    gt    core={gc.get(g) or '?':<16} {Path(g).name}")
    if len(gt_free) > limit:
        print(f"    ... (+{len(gt_free) - limit} more)")


def _apply_pairmap_overlay(preds, csv_path):
    """Fill GT (and optionally id) for predictions named in a pair-map, matched by pred path."""
    path = Path(csv_path).expanduser()
    if not path.exists():
        warn(f"pair-map {csv_path!r} not found.")
        return
    by_path = {str(Path(p.pred).resolve()): p for p in preds}
    by_name = {Path(p.pred).name: p for p in preds}
    base = path.resolve().parent
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        if "pred" not in cols:
            warn(f"pair-map {csv_path}: no 'pred' column; ignored.")
            return
        applied = 0
        for row in reader:
            predval = row[cols["pred"]].strip()
            rp = _resolve(base, predval)
            pr = (by_path.get(str(rp.resolve())) if rp else None) or by_name.get(Path(predval).name)
            if pr is None:
                continue
            if "gt" in cols and row[cols["gt"]].strip():
                pr.gt = _resolve(base, row[cols["gt"]].strip())
                pr.gt_forced = True
            if "id" in cols and row[cols["id"]].strip():
                pr.set_key(row[cols["id"]].strip())
            applied += 1
        log(f"pair-map overlay applied to {applied} prediction(s).")


def interactive_resolve(preds, gt_dir, gt_paths, args):
    """Human-in-the-loop cap on the ladder: when GT was provided but some predictions did not
    match, let the user continue GT-less, supply a pair-map / regex, re-scan, or abort. Falls
    back to continue-GT-less (or abort under --strict) when stdin is not interactive."""
    while True:
        unmatched = [p for p in preds if p.gt is None]
        if not unmatched:
            return
        interactive = sys.stdin.isatty() and not args.non_interactive
        if not interactive:
            if args.strict:
                die(f"{len(unmatched)} prediction(s) have no ground-truth match (--strict).")
            warn(f"{len(unmatched)} prediction(s) have no GT match; continuing "
                 f"confidence-only for them.")
            return
        _preview_unmatched(unmatched, preds, gt_paths, args.id_regex)
        print("[run_eval] how to proceed?\n"
              "    1) continue without GT for these (confidence-only)\n"
              "    2) provide/update a pair-map CSV\n"
              "    3) enter a new --id-regex\n"
              "    4) re-scan the GT folder (after renaming files)\n"
              "    5) abort")
        try:
            choice = input("[run_eval] choose 1-5: ").strip()
        except EOFError:
            choice = "1"
        if choice == "1":
            warn(f"continuing confidence-only for {len(unmatched)} prediction(s).")
            return
        elif choice == "2":
            _apply_pairmap_overlay(preds, input("    pair-map CSV path: ").strip())
            resolve_gt(preds, gt_paths, args.id_regex)
        elif choice == "3":
            args.id_regex = input("    new --id-regex (capture group 1 = id): ").strip() or None
            resolve_gt(preds, gt_paths, args.id_regex)
        elif choice == "4":
            input("    rename files in the GT folder, then press Enter to re-scan...")
            gt_paths = sorted(Path(gt_dir).glob("*.cif"))
            resolve_gt(preds, gt_paths, args.id_regex)
        elif choice == "5":
            die("aborted at ground-truth matching.")
        else:
            print("    invalid choice.")


def build_pairs_from_map(csv_path):
    csv_path = Path(csv_path).resolve()
    base = csv_path.parent
    pairs = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        if "pred" not in cols:
            die(f"pair-map {csv_path} must have a 'pred' column (got {reader.fieldnames}).")
        for rownum, row in enumerate(reader, start=2):
            pred = _resolve(base, row[cols["pred"]].strip())
            if pred is None:
                warn(f"pair-map row {rownum}: pred not found ({row[cols['pred']]!r}); skipping.")
                continue
            pr = Pair(pred)
            key = (row[cols["id"]].strip() if "id" in cols and row[cols["id"]].strip()
                   else pr.raw_core)
            pr.set_key(key)
            if "gt" in cols and row[cols["gt"]].strip():
                pr.gt = _resolve(base, row[cols["gt"]].strip())
                pr.gt_forced = True
            conf_source = (_resolve(base, row[cols["confidence"]].strip())
                           if "confidence" in cols and row[cols["confidence"]].strip()
                           else Path(pred).parent)
            data_source = (_resolve(base, row[cols["data_json"]].strip())
                           if "data_json" in cols and row[cols["data_json"]].strip()
                           else Path(pred).parent)
            resolve_sidecars(pr, conf_source, data_source)
            for attr in ("motif", "domain"):
                if attr in cols and row[cols[attr]].strip():
                    setattr(pr, attr, _resolve(base, row[cols[attr]].strip()))
            pairs.append(pr)
    return pairs


def build_pairs_from_dirs(args):
    pred_dir = Path(args.pred).resolve()
    if not pred_dir.is_dir():
        die(f"--pred {pred_dir} is not a directory.")
    gt_dir = Path(args.gt).resolve() if args.gt else None
    data_dir = Path(args.data_json).resolve() if args.data_json else None
    conf_dir = Path(args.confidence).resolve() if args.confidence else None

    preds = [Pair(p) for p in sorted(pred_dir.glob("*.cif"))]
    if not preds:
        return preds
    gt_paths = sorted(gt_dir.glob("*.cif")) if gt_dir and gt_dir.is_dir() else []
    if gt_dir and not gt_paths:
        warn(f"--gt {gt_dir} contains no *.cif; proceeding confidence-only.")

    resolve_gt(preds, gt_paths, args.id_regex)
    if gt_paths:
        interactive_resolve(preds, gt_dir, gt_paths, args)

    for pr in preds:
        resolve_sidecars(pr, conf_dir or Path(pr.pred).parent,
                         data_dir or Path(pr.pred).parent)
    return preds


# --------------------------------------------------------------------------------------
# staging
# --------------------------------------------------------------------------------------
def stage(pairs, workdir, classification, rewrite_ids, cutoff_date):
    """Copy inputs into <workdir>/data/<classification>/{af,pdb} under the fold_<ID>_s1_*
    naming convention. Never moves user files. Returns the list of staged Pair records and
    writes staging_manifest.csv."""
    data_root = workdir / "data" / classification
    af_dir = data_root / "af"
    pdb_dir = data_root / "pdb"
    af_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    seen = {}
    staged = []
    for pr in pairs:
        if not pr.id_token:
            warn(f"{Path(pr.pred).name}: empty id after sanitising; skipping.")
            continue
        if pr.join_id in seen:
            warn(f"id collision on {pr.join_id} ({Path(pr.pred).name} vs "
                 f"{Path(seen[pr.join_id]).name}); skipping the duplicate.")
            continue
        seen[pr.join_id] = pr.pred

        # prediction -> af/fold_<ID>_s1_model_0.cif
        pred_dst = af_dir / f"fold_{pr.id_token}_s1_model_0.cif"
        shutil.copyfile(pr.pred, pred_dst)

        # confidence sidecars -> af/fold_<ID>_s1_{summary_confidences,full_data}_<k>.json
        pr.staged_confidence = _stage_confidence(pr, af_dir)

        # AF3 input sidecar (MSA/templates) -> af/<ID>_s1_data.json  (first-4-char id rule)
        if pr.data_json and Path(pr.data_json).exists():
            shutil.copyfile(pr.data_json, af_dir / f"{pr.id_token}_s1_data.json")
            pr.staged_data_json = True

        # ground truth -> pdb/<ID>.cif  (+ force _entry.id so the DB join works)
        if pr.gt and Path(pr.gt).exists():
            gt_dst = pdb_dir / f"{pr.id_token}.cif"
            shutil.copyfile(pr.gt, gt_dst)
            if rewrite_ids:
                rewrite_entry_id(gt_dst, pr.join_id)
            date = read_deposition_date(gt_dst)
            if date:
                pr.era = "a" if date >= cutoff_date else "b"  # a = post-cutoff, b = pre
        staged.append(pr)

    _write_manifest(workdir, staged, classification)
    return staged


def _stage_confidence(pr, af_dir):
    """Copy this prediction's resolved confidence JSONs, renamed to the fold_<ID>_s1_<kind>_<k>
    form af3Metrics globs (normalising old and current AF3 namings). Returns the count staged."""
    staged = 0
    seen = set()
    for cand in pr.confidence_files:
        kv = classify_confidence(Path(cand).name)
        if kv is None or kv in seen:   # skip non-confidence and duplicate (kind,index) files
            continue
        seen.add(kv)
        kind, k = kv
        shutil.copyfile(cand, af_dir / f"fold_{pr.id_token}_s1_{kind}_{k}.json")
        staged += 1
    return staged


def _write_manifest(workdir, staged, classification):
    with open(workdir / "staging_manifest.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "join_id", "classification", "pred", "gt",
                    "confidence_files", "has_data_json", "era"])
        for pr in staged:
            w.writerow([
                pr.id_token, pr.join_id, classification,
                pr.pred, pr.gt or "",
                pr.staged_confidence,
                "yes" if pr.staged_data_json else "no",
                pr.era or "",
            ])


# --------------------------------------------------------------------------------------
# config.json management (per-run isolation)
# --------------------------------------------------------------------------------------
@contextmanager
def run_config(classification, workdir):
    """Temporarily point database/config.json at this run (fresh DB), restoring the original
    on exit so the repo's rbpDatabase.db / config.json are never clobbered."""
    backup = None
    if CONFIG_PATH.exists():
        backup = CONFIG_PATH.read_text()
    cfg = {
        "classification": classification,
        "base_data_path": str((workdir / "data").resolve()),
        "db_path": str((workdir / "run.db").resolve()),
    }
    CONFIG_PATH.write_text(json.dumps(cfg, indent=4))
    log(f"config.json -> classification={classification}, db={cfg['db_path']}")
    try:
        yield cfg
    finally:
        if backup is not None:
            CONFIG_PATH.write_text(backup)
            log("restored original database/config.json")
        else:
            CONFIG_PATH.unlink(missing_ok=True)


# --------------------------------------------------------------------------------------
# step execution
# --------------------------------------------------------------------------------------
def _base_env():
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + existing if existing else "")
    return env


def run_ambient(script_rel, script_args, cwd=None):
    """Run a repo script in the current (eval) conda env, PYTHONPATH=<repo>."""
    cmd = [sys.executable, str(REPO_ROOT / script_rel), *script_args]
    log("exec: " + " ".join(cmd))
    res = subprocess.run(cmd, cwd=str(cwd or REPO_ROOT), env=_base_env(), check=False)
    return res.returncode == 0


def run_ost(script_rel, script_args, ost_env):
    """Run an OST script inside the 'ost' conda env (mirrors mainExecution.run_ost)."""
    conda = shutil.which("conda") or "conda"
    cmd = [conda, "run", "--no-capture-output", "-n", ost_env,
           "python", str(REPO_ROOT / script_rel), *script_args]
    # redact any secret (e.g. SWISS-MODEL --token) before logging so keys never land in logs
    safe = ["***REDACTED***" if i > 0 and cmd[i - 1] == "--token" else c
            for i, c in enumerate(cmd)]
    log("exec (ost): " + " ".join(safe))
    res = subprocess.run(cmd, cwd=str(REPO_ROOT), env=_base_env(), check=False)
    return res.returncode == 0


# --------------------------------------------------------------------------------------
# pipeline construction
# --------------------------------------------------------------------------------------
class Step:
    def __init__(self, sid, phase, rank, runner, requires=(), note=""):
        self.sid = sid          # stable name (usable in --metrics)
        self.phase = phase
        self.rank = rank        # introduced at this preset rank; 99 => flag-only
        self.runner = runner    # callable() -> bool
        self.requires = tuple(requires)
        self.note = note


def build_steps(args, classification, caps, af_path, pdb_path, dssr_dir, preset_rank,
                selected_names):
    is_pr = classification == "protein_rna"
    token = caps.get("_token_value", "")
    steps = []

    def add(sid, phase, rank, runner, requires=(), note=""):
        steps.append(Step(sid, phase, rank, runner, requires, note))

    # ---- parse ----
    add("parse_gt", "parse", 1,
        lambda: run_ambient("database/parsePDBCIF2DB.py", [pdb_path]),
        requires=["gt"], note="ground-truth mmCIFs")
    add("parse_af", "parse", 1,
        lambda: run_ambient("database/parseAFcif2DB.py", [af_path]))

    # ---- AF3 confidence (pred-only, no GT needed) ----
    add("af3", "confidence", 1,
        lambda: run_ambient("evaluation/af3Metrics_new.py", [af_path, "--no-plots"]),
        note="AF3 *_summary_confidences_/_full_data_ JSONs")

    # ---- structural accuracy (US-align) ----
    # quick preset: extractEvalMetrics computes its own LDDT (no OST). Otherwise pass -lddt so
    # LDDT comes from the OST scripts (basic/full), matching mainExecution.
    quick_only = selected_names is None and preset_rank <= PRESET_RANK["quick"]
    eval_args = [] if quick_only else ["-lddt"]
    add("eval_metrics", "accuracy", 1,
        lambda: run_ambient("evaluation/extractEvalMetrics.py", eval_args, cwd=EVAL_DIR),
        requires=["usalign", "gt"], note="US-align + ground truth")

    # ---- OST scoring (basic+) ----
    add("ost_lddt", "ost", 2,
        lambda: run_ost("evaluation/ost_extractAllAtomLDDT.py",
                        [f"pred_{classification}"], args.ost_env),
        requires=["ost", "gt"], note="'ost' conda env + ground truth")
    add("ost_allchains", "ost", 2,
        lambda: run_ost("evaluation/ost_calcRNA_AllChains.py",
                        [] if is_pr else ["--lddt"], args.ost_env),
        requires=["ost", "gt"], note="'ost' conda env + ground truth")

    # ---- interface + binding + RBP (full) ----
    # extractInteractions.py imports pymol at module load, so it is gated on PyMOL too.
    add("interactions_af", "interface", 3,
        lambda: run_ambient("evaluation/extractInteractions.py", [af_path]),
        requires=["pymol"], note="PyMOL")
    if is_pr:
        add("interface_area_af", "interface", 3,
            lambda: run_ambient("evaluation/extractInterfaceArea.py", [af_path]),
            requires=["pymol", "freesasa"], note="PyMOL + freesasa")
        add("binding_af", "interface", 3,
            lambda: run_ambient("evaluation/estimateBindingAffinity.py",
                                [f"pred_{classification}"]))
    if caps.get("gt"):
        add("interactions_gt", "interface", 3,
            lambda: run_ambient("evaluation/extractInteractions.py", [pdb_path]),
            requires=["gt", "pymol"], note="ground truth + PyMOL")
        if is_pr:
            add("interface_area_gt", "interface", 3,
                lambda: run_ambient("evaluation/extractInterfaceArea.py", [pdb_path]),
                requires=["gt", "pymol", "freesasa"], note="PyMOL + freesasa")
            add("binding_gt", "interface", 3,
                lambda: run_ambient("evaluation/estimateBindingAffinity.py",
                                    [f"exp_{classification}"]), requires=["gt"])
    if is_pr and not args.no_rbp:
        add("ost_rbp", "ost", 3,
            lambda: run_ost("evaluation/ost_RBPinterfaces.py",
                            ["--token", token], args.ost_env),
            requires=["ost", "gt", "token"],
            note="'ost' env + ground truth + SWISS-MODEL token (or pass -noRBP to skip)")

    # ---- metadata add-ons (opt-in flags, rank 99 so presets never auto-include) ----
    if args.msa:
        add("msa", "metadata", 99,
            lambda: run_ambient("utils/addMSA_fromJSON.py", [af_path]),
            note="AF3 *_data.json sidecars (4-char PDB ids only)")
    if args.templates:
        add("templates", "metadata", 99,
            lambda: run_ambient("database/extractTemplatesAF3.py", [af_path]),
            note="AF3 *_data.json sidecars (4-char PDB ids only)")
    if args.db_update:
        _add_db_update_steps(add, dssr_dir)
    if args.domains:
        if not args.domain_tsv:
            warn("--domains given without --domain-tsv; domains step will be skipped.")
        add("domains", "metadata", 99,
            lambda: run_ambient("database/addDomain.py",
                                (["--tsv", args.domain_tsv] if args.domain_tsv else [])),
            requires=["_domain_tsv"], note="a domain TSV via --domain-tsv")
    if args.protein_family:
        pf_args = []
        if args.family_email:
            pf_args += ["--email", args.family_email]
        add("protein_family", "metadata", 99,
            lambda: run_ambient("database/insert_proteinFamily2DB.py", pf_args),
            requires=["gt"], note="internet + a contact email (--family-email)")
    if args.motif_csv:
        add("motif", "metadata", 99,
            lambda: run_ambient("rnaCompete/motif_from_csv.py", ["--csv", args.motif_csv]),
            requires=["_motif_csv"], note="a motif CSV via --motif-csv")
        caps["_motif_csv"] = bool(Path(args.motif_csv).exists())

    # capability flag for the domain TSV (validated as a normal requirement)
    caps["_domain_tsv"] = bool(args.domain_tsv and Path(args.domain_tsv).exists())

    # ---- selection: preset rank OR explicit --metrics names ----
    if selected_names is not None:
        chosen = [s for s in steps if s.sid in selected_names]
        missing = selected_names - {s.sid for s in steps}
        for name in sorted(missing):
            warn(f"--metrics: unknown/unavailable step '{name}' (ignored).")
    else:
        chosen = [s for s in steps if s.rank <= preset_rank]
        # metadata flag-steps (rank 99) are always kept when their flag added them
        chosen += [s for s in steps if s.rank == 99]
    return chosen


def _add_db_update_steps(add, dssr_dir):
    """The +dssr --db-update metadata block, mirroring mainExecution ordering."""
    dssr_arg = str(dssr_dir) if dssr_dir else ""
    add("dssr_conformation", "metadata", 99,
        lambda: run_ambient("evaluation/dssrConformation.py", [dssr_arg]), requires=["dssr"],
        note="a populated dssr/ folder")
    add("rna_metrics", "metadata", 99,
        lambda: run_ambient("evaluation/extractRnaMetrics.py", [dssr_arg]), requires=["dssr"],
        note="a populated dssr/ folder")
    add("dssr_pairs", "metadata", 99,
        lambda: run_ambient("evaluation/extractDSSRpairs.py", [dssr_arg]), requires=["dssr"],
        note="a populated dssr/ folder")
    add("secondary_feat", "metadata", 99,
        lambda: run_ambient("evaluation/extractSecondaryFeat.py", [dssr_arg]), requires=["dssr"],
        note="a populated dssr/ folder")
    add("add_ph", "metadata", 99,
        lambda: run_ambient("database/addPH.py", []))
    add("electrostatics", "metadata", 99,
        lambda: run_ambient("evaluation/extractElectrostaticPotential.py", []))
    add("gc_complexity", "metadata", 99,
        lambda: run_ambient("evaluation/extractGCcomplexity.py", []))
    add("ph_correlation", "metadata", 99,
        lambda: run_ambient("evaluation/phCorrelation.py", ["db_update"]))
    add("ions_ligands", "metadata", 99,
        lambda: run_ambient("database/add_ions_ligands_to_db.py", []))
    add("exp_details", "metadata", 99,
        lambda: run_ambient("database/addExpDetails.py", []))


# --------------------------------------------------------------------------------------
# execution driver
# --------------------------------------------------------------------------------------
def execute_pipeline(steps, caps, args):
    ran, skipped, failed = [], [], []
    current_phase = None
    for step in steps:
        if step.phase != current_phase:
            current_phase = step.phase
            log("=" * 60)
            log(f"phase: {current_phase}")
        missing = [r for r in step.requires if not caps.get(r)]
        if missing:
            reason = step.note or ", ".join(missing)
            if args.strict:
                die(f"{step.sid} requires {reason} (--strict).")
            skip(step.sid, reason)
            skipped.append(step.sid)
            continue
        if args.dry_run:
            log(f"DRY-RUN would execute: {step.sid}")
            ran.append(step.sid)
            continue
        ok = step.runner()
        if ok:
            ran.append(step.sid)
        else:
            warn(f"{step.sid} failed (continuing).")
            failed.append(step.sid)
    return ran, skipped, failed


def export_csvs(classification, workdir):
    log("=" * 60)
    log("exporting DB tables to CSV")
    run_ambient("database/exportDB2CSV.py", [f"pred_{classification}"])
    run_ambient("database/exportDB2CSV.py", [f"exp_{classification}"])
    results_dir = REPO_ROOT / "results"
    for tbl in (f"pred_{classification}", f"exp_{classification}"):
        src = results_dir / f"{tbl}.csv"
        if src.exists():
            shutil.copyfile(src, workdir / f"{tbl}.csv")


# --------------------------------------------------------------------------------------
# era split (optional)
# --------------------------------------------------------------------------------------
def era_split(run_db, workdir, staged, classification):
    """Post-process run.db into a<name>.db / b<name>.db by deposition era so downstream
    preprocessing (infer_post_cutoff_from_db_tag: tag.startswith('a')) can split pre/post.
    Rows with no era are dropped from both (their date was unknown)."""
    eras = {"a": set(), "b": set()}
    for pr in staged:
        if pr.era in eras:
            eras[pr.era].add(pr.join_id)
    present = [e for e in ("a", "b") if eras[e]]
    if len(present) < 2:
        log("era-split: fewer than two eras present; keeping single run.db only.")
        return []
    pred_tbl, exp_tbl = f"pred_{classification}", f"exp_{classification}"
    outputs = []
    for era in present:
        keep = eras[era]
        dst = workdir / f"{era}_run.db"
        shutil.copyfile(run_db, dst)
        conn = sqlite3.connect(dst)
        cur = conn.cursor()
        placeholders = ",".join("?" for _ in keep) or "''"
        for tbl, col in ((pred_tbl, "exp_db_id"), (exp_tbl, "PDBId")):
            try:
                cur.execute(
                    f"DELETE FROM {tbl} WHERE UPPER({col}) NOT IN ({placeholders})",
                    tuple(keep))
            except sqlite3.Error as exc:
                warn(f"era-split {tbl}: {exc}")
        conn.commit()
        conn.close()
        log(f"era-split: wrote {dst.name} ({len(keep)} complexes)")
        outputs.append(dst)
    return outputs


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def build_arg_parser():
    p = argparse.ArgumentParser(
        prog="run_eval.py",
        description="Drop-in AF3InterfaceEval wrapper: predictions (+ optional GT) -> SQLite DB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # quick (parse + AF3 confidence + US-align TM/RMSD/LDDT, no OST/DSSR)
  python run_eval.py --pred ./preds --gt ./gt --out ./run --metrics quick
  # basic (+ OST all-atom LDDT / RNA AllChains)
  python run_eval.py --pred ./preds --gt ./gt --out ./run --metrics basic
  # full (+ interfaces, binding, RBP DockQ/QS) with a SWISS-MODEL token
  python run_eval.py --pred ./preds --gt ./gt --out ./run --metrics full --token TOKEN
  # arbitrary filenames via a pair-map CSV (pred,gt[,id,data_json,confidence,motif,domain])
  python run_eval.py --pair-map pairs.csv --out ./run --metrics basic
  # confidence-only (no GT): AF3 metrics only, GT-dependent steps skipped
  python run_eval.py --pred ./preds --out ./run --metrics quick
""")
    io = p.add_argument_group("inputs / outputs")
    io.add_argument("--pred", help="Folder of predicted *.cif (required unless --pair-map).")
    io.add_argument("--gt", help="Folder of deposited ground-truth *.cif (optional; "
                                 "absent => confidence-only, GT-dependent steps skipped).")
    io.add_argument("--out", "--workdir", dest="out", required=True,
                    help="Run working directory (staging + run.db + CSVs land here).")
    io.add_argument("--pair-map", dest="pair_map",
                    help="CSV mapping arbitrary filenames to inputs; columns "
                         "pred,gt[,id,data_json,confidence,motif,domain].")
    io.add_argument("--data-json", dest="data_json",
                    help="Folder of AF3 *_data.json sidecars (MSA/templates).")
    io.add_argument("--confidence", help="Folder of AF3 confidence JSONs "
                                         "(else auto-detected beside each prediction).")
    io.add_argument("--dssr", dest="dssr_dir",
                    help="Folder of DSSR output for +dssr metadata (--db-update).")

    idg = p.add_argument_group("identity / classification")
    idg.add_argument("--classification", choices=(*CLASSIFICATIONS, "auto"), default="auto",
                     help="Force protein_rna / rna_rna, or auto-detect from predicted CIFs.")
    idg.add_argument("--id-regex", dest="id_regex",
                     help="Regex applied to each prediction stem to extract the id "
                          "(group 1 if present, else whole match).")
    idg.add_argument("--rewrite-entry-id", dest="rewrite_entry_id",
                     action="store_true", default=True,
                     help="Force GT _entry.id to the shared id so the DB join works (default on).")
    idg.add_argument("--no-rewrite-entry-id", dest="rewrite_entry_id", action="store_false",
                     help="Leave GT _entry.id untouched (only if filenames already match PDBId).")
    idg.add_argument("--cutoff-date", dest="cutoff_date", default="2021-09-30",
                     help="Era cutoff (YYYY-MM-DD); GT on/after => 'a' (post), before => 'b' (pre).")
    idg.add_argument("--era-split", dest="era_split", action="store_true",
                     help="Also emit a_run.db / b_run.db split by GT deposition era.")

    sel = p.add_argument_group("metric selection")
    sel.add_argument("--metrics", default="basic",
                     help="Preset (quick|basic|full) or a comma-separated list of step ids "
                          "(e.g. parse_af,af3,eval_metrics).")
    sel.add_argument("--token", help="SWISS-MODEL API token for ost_RBPinterfaces "
                                     "(overrides SWISSMODEL_OST_TOKEN).")
    sel.add_argument("-noRBP", "--no-rbp", dest="no_rbp", action="store_true",
                     help="Skip ost_RBPinterfaces (protein-RNA RBP DockQ/QS).")
    sel.add_argument("--ost-env", dest="ost_env", default="ost",
                     help="conda env name for OST scoring steps (default: ost).")

    meta = p.add_argument_group("metadata add-ons (opt-in, gated on inputs)")
    meta.add_argument("--msa", action="store_true", help="MSA depth from *_data.json.")
    meta.add_argument("--templates", action="store_true", help="Template counts from *_data.json.")
    meta.add_argument("--db-update", dest="db_update", action="store_true",
                      help="pH/ions/ligands/GC-MFE/secondary-structure/electrostatics/DSSR block.")
    meta.add_argument("--domains", action="store_true", help="Domain annotations (needs --domain-tsv).")
    meta.add_argument("--domain-tsv", dest="domain_tsv", help="Domain TSV for --domains.")
    meta.add_argument("--protein-family", dest="protein_family", action="store_true",
                      help="Protein family/organism via web APIs.")
    meta.add_argument("--family-email", dest="family_email",
                      help="Contact email for the InterPro/EBI family API.")
    meta.add_argument("--motif-csv", dest="motif_csv",
                      help="CSV of RNA-motif inputs (pred_id,rna_sequence,exp_motif[,pred_motif,"
                           "interface_area,protein_length]) -> RNAmotif_score / RNAMotif_inExp / "
                           "InterfaceArea_norm.")

    run = p.add_argument_group("execution")
    run.add_argument("--strict", action="store_true",
                     help="Abort (instead of skip-and-log) when a step's tool/input is missing, "
                          "or when predictions have no ground-truth match.")
    run.add_argument("--non-interactive", "--yes", dest="non_interactive", action="store_true",
                     help="Never prompt: unmatched predictions continue confidence-only "
                          "(or abort under --strict).")
    run.add_argument("--dry-run", dest="dry_run", action="store_true",
                     help="Print the resolved pipeline without executing anything.")
    return p


def parse_metrics_selection(value):
    """Return (preset_rank, selected_names_or_None)."""
    value = value.strip()
    if value in PRESET_RANK:
        return PRESET_RANK[value], None
    names = {v.strip() for v in value.split(",") if v.strip()}
    if not names:
        die(f"--metrics: empty selection {value!r}.")
    return PRESET_RANK["full"], names  # rank irrelevant when explicit names are given


def resolve_classification(args, pairs):
    if args.classification != "auto":
        return args.classification
    for pr in pairs:
        cls = detect_classification(pr.pred)
        if cls:
            log(f"auto-detected classification: {cls} (from {Path(pr.pred).name})")
            return cls
    die("could not auto-detect classification; pass --classification protein_rna|rna_rna.")


def main():
    args = build_arg_parser().parse_args()
    if not args.pair_map and not args.pred:
        die("provide --pred (a predictions folder) or --pair-map.")
    if args.no_rbp:
        args.token = None
        os.environ.pop("SWISSMODEL_OST_TOKEN", None)

    workdir = Path(args.out).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    log(f"workdir: {workdir}")

    pairs = (build_pairs_from_map(args.pair_map) if args.pair_map
             else build_pairs_from_dirs(args))
    if not pairs:
        die("no (prediction) inputs found.")
    log(f"resolved {len(pairs)} prediction(s); "
        f"{sum(1 for p in pairs if p.gt)} with matching ground truth.")

    classification = resolve_classification(args, pairs)
    has_gt = any(p.gt for p in pairs)
    dssr_dir = Path(args.dssr_dir).resolve() if args.dssr_dir else None

    staged = stage(pairs, workdir, classification, args.rewrite_entry_id, args.cutoff_date)
    if not staged:
        die("nothing staged (all inputs skipped).")
    af_path = str((workdir / "data" / classification / "af").resolve())
    pdb_path = str((workdir / "data" / classification / "pdb").resolve())

    caps = detect_capabilities(args, has_gt, dssr_dir)
    log("capabilities: " + ", ".join(
        f"{k}={'yes' if v else 'no'}" for k, v in caps.items() if not k.startswith("_")))

    preset_rank, selected_names = parse_metrics_selection(args.metrics)

    with run_config(classification, workdir):
        steps = build_steps(args, classification, caps, af_path, pdb_path, dssr_dir,
                            preset_rank, selected_names)
        log(f"pipeline ({len(steps)} steps): " + ", ".join(s.sid for s in steps))
        ran, skipped, failed = execute_pipeline(steps, caps, args)
        if not args.dry_run:
            export_csvs(classification, workdir)
            outputs = []
            if args.era_split and has_gt:
                outputs = era_split(workdir / "run.db", workdir, staged, classification)

    log("=" * 60)
    log(f"done. ran={len(ran)} skipped={len(skipped)} failed={len(failed)}")
    if skipped:
        log("skipped: " + ", ".join(skipped))
    if failed:
        log("failed: " + ", ".join(failed))
    if not args.dry_run:
        log(f"run DB: {workdir / 'run.db'}")
        log(f"manifest: {workdir / 'staging_manifest.csv'}")


if __name__ == "__main__":
    main()
