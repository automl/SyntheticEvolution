"""Run DSSR over a predicted structure and read its base pairs back out, matching the JSON layout
already consumed by plot_secondary_structures.py."""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import CommandSpec, run_command

# AF3 writes the model under a few different names depending on version and on whether the run came
# from the standalone pipeline or the server; check the specific ones before falling back to a glob.
_CIF_PATTERNS = ["*_model_0.cif", "*_model.cif", "*model*.cif", "*.cif"]


def resolve_job_dir(af3_root: Path, job_id: str) -> Optional[Path]:
    """Find the directory AF3 actually wrote for a job.

    AF3 derives the output directory from the job name and sanitises it, so the result is not
    always the lowercased name the pipeline predicts. Try the predicted name first, then a
    case-insensitive match, then a sanitised-character match.
    """
    if not af3_root.exists():
        return None
    predicted = af3_root / job_id.lower()
    if predicted.is_dir():
        return predicted
    wanted = job_id.lower()
    candidates = [d for d in af3_root.iterdir() if d.is_dir()]
    for directory in candidates:
        if directory.name.lower() == wanted:
            return directory
    squashed = "".join(c for c in wanted if c.isalnum())
    for directory in candidates:
        if "".join(c for c in directory.name.lower() if c.isalnum()) == squashed:
            return directory
    return None


def find_model_cif(job_output_dir: Optional[Path]) -> Optional[Path]:
    """Locate the ranked-first model in an AF3 output directory."""
    if job_output_dir is None or not job_output_dir.exists():
        return None
    for pattern in _CIF_PATTERNS:
        hits = sorted(job_output_dir.rglob(pattern))
        if hits:
            return hits[0]
    return None


def find_summary_confidences(job_output_dir: Path) -> Optional[Path]:
    for pattern in ["*summary_confidences_0.json", "*summary_confidences.json"]:
        hits = sorted(job_output_dir.rglob(pattern))
        if hits:
            return hits[0]
    return None


def read_confidences(path: Optional[Path]) -> Dict[str, float]:
    """Pull the headline AF3 confidence numbers, tolerating naming differences across versions."""
    if path is None or not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception as exc:
        logging.warning("Could not read %s: %s", path, exc)
        return {}
    out = {}
    for key in ("ptm", "iptm", "ranking_score", "fraction_disordered", "has_clash"):
        value = data.get(key)
        if isinstance(value, (int, float)):
            out[key] = float(value)
    return out


def run_dssr(cif: Path, out_json: Path, spec: CommandSpec, dry_run: bool = False,
             log_path: Optional[Path] = None) -> bool:
    """Annotate one CIF. Returns True when the JSON exists afterwards."""
    out_json.parent.mkdir(parents=True, exist_ok=True)
    command = spec.render(cif=str(cif), out_json=str(out_json), out_dir=str(out_json.parent))
    code = run_command(command, spec.timeout_s, dry_run=dry_run, log_path=log_path)
    if dry_run:
        return True
    return code == 0 and out_json.exists()


def parse_pairs(dssr_json: Path, chain: Optional[str] = None) -> Tuple[str, List[Tuple[int, int]]]:
    """Return (sequence, 0-based pairs) for one RNA chain of a DSSR annotation.

    Positions come from `index_chain`, which DSSR numbers from 1 within each chain, so pairs are
    only meaningful within a single chain; pairs spanning two chains are dropped.
    """
    data = json.loads(Path(dssr_json).read_text())
    chains = data.get("chains") or {}
    if not chains:
        raise ValueError(f"{dssr_json} has no 'chains' key; DSSR found no annotatable structure")

    if chain is None:
        # Longest chain, which for an RNA-protein complex prediction is the RNA of interest
        # whenever DSSR only emits nucleic acid chains (it does).
        chain = max(chains, key=lambda k: len(chains[k].get("bseq", "")))
    if chain not in chains:
        raise KeyError(f"chain {chain} not in {sorted(chains)}")

    sequence = chains[chain].get("bseq", "").replace("&", "").upper()

    nts = {nt["nt_id"]: nt for nt in data.get("nts", [])}
    pairs: List[Tuple[int, int]] = []
    for pair in data.get("pairs", []):
        nt1, nt2 = nts.get(pair.get("nt1")), nts.get(pair.get("nt2"))
        if not nt1 or not nt2:
            continue
        if nt1.get("chain_name") != chain or nt2.get("chain_name") != chain:
            continue
        i, j = int(nt1["index_chain"]) - 1, int(nt2["index_chain"]) - 1
        if 0 <= i < len(sequence) and 0 <= j < len(sequence) and i != j:
            pairs.append((min(i, j), max(i, j)))
    return sequence, sorted(set(pairs))
