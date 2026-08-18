"""Run DSSR over a predicted structure and read its base pairs back out, matching the JSON layout
already consumed by plot_secondary_structures.py."""

import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

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
             log_path: Optional[Path] = None, scratch: Optional[Path] = None,
             root: Optional[Path] = None) -> bool:
    """Annotate one CIF. Returns True when the JSON exists afterwards.

    DSSR drops fixed-name auxiliary files (dssr-pairs.pdb, dssr-stems.pdb, ...) into its working
    directory, so it is given a private one. Sharing a directory makes concurrent runs fail with
    "remove_file failed: dssr-pairs.pdb" and, worse, exit non-zero having written nothing. Because
    that moves the command out of the invocation directory, templates must not use paths relative
    to it -- use {root} for anything inside the repository.
    """
    out_json.parent.mkdir(parents=True, exist_ok=True)
    # Absolute, because the command runs in a private working directory: a path relative to where
    # the pipeline was invoked would no longer resolve.
    cif, out_json = cif.resolve(), out_json.resolve()
    command = spec.render(cif=str(cif), out_json=str(out_json), out_dir=str(out_json.parent),
                          root=str(root or Path.cwd()))
    if scratch is not None:
        scratch.mkdir(parents=True, exist_ok=True)
    code = run_command(command, spec.timeout_s, dry_run=dry_run, log_path=log_path, cwd=scratch)
    if dry_run:
        return True
    return code == 0 and out_json.exists()


# DSSR's `name` for a pair; the canonical set is what the published evaluation scores on.
_CANONICAL_NAMES = {"WC", "Wobble"}


def _bseq(chains: Dict[str, dict], key: str) -> str:
    return chains[key].get("bseq", "").replace("&", "").upper()


# DSSR keys `chains` by model and chain ("m1_chain_A") but reports the bare chain id ("A") as
# `chain_name` on every nucleotide. Comparing the two directly matches nothing, which silently
# yields an empty pair list rather than an error.
_CHAIN_KEY = re.compile(r"^m\d+_chain_(?P<name>.+)$")


def chain_name(key: str) -> str:
    match = _CHAIN_KEY.match(key)
    return match.group("name") if match else key


def select_chain(chains: Dict[str, dict], prefer: Optional[Sequence[str]] = None,
                 match_sequence: Optional[str] = None) -> str:
    """Pick the chain a job is actually about.

    `prefer` are the AF3 chain ids we wrote the SHS into, and are authoritative for a prediction.
    A ground-truth structure carries the depositor's own chain labelling, which need not agree
    with ours, so there the sequence is what identifies the chain. Falling back to the longest
    chain is only right for a single-RNA structure; in a complex holding two different RNAs it
    silently scores the other one.
    """
    by_name = {chain_name(k): k for k in chains}
    for candidate in prefer or []:
        if candidate in chains:
            return candidate
        if candidate in by_name:
            return by_name[candidate]
    if match_sequence:
        wanted = match_sequence.upper()
        for key in chains:
            if _bseq(chains, key) == wanted:
                return key
        same_length = [k for k in chains if len(_bseq(chains, k)) == len(wanted)]
        if len(same_length) == 1:
            return same_length[0]
        if same_length:
            return sorted(same_length)[0]
    return max(chains, key=lambda k: len(_bseq(chains, k)))


def parse_pairs(dssr_json: Path, chain: Optional[str] = None,
                prefer: Optional[Sequence[str]] = None,
                match_sequence: Optional[str] = None,
                canonical_only: bool = False) -> Tuple[str, List[Tuple[int, int]]]:
    """Return (sequence, 0-based pairs) for one RNA chain of a DSSR annotation.

    Positions come from `index_chain`, which DSSR numbers from 1 within each chain, so pairs are
    only meaningful within a single chain; pairs spanning two chains are dropped.
    """
    data = json.loads(Path(dssr_json).read_text())
    chains = data.get("chains") or {}
    if not chains:
        raise ValueError(f"{dssr_json} has no 'chains' key; DSSR found no annotatable structure")

    if chain is None:
        chain = select_chain(chains, prefer=prefer, match_sequence=match_sequence)
    if chain not in chains:
        by_name = {chain_name(k): k for k in chains}
        if chain not in by_name:
            raise KeyError(f"chain {chain} not in {sorted(chains)}")
        chain = by_name[chain]

    sequence = _bseq(chains, chain)
    wanted = chain_name(chain)

    nts = {nt["nt_id"]: nt for nt in data.get("nts", [])}
    pairs: List[Tuple[int, int]] = []
    for pair in data.get("pairs", []):
        if canonical_only and pair.get("name") not in _CANONICAL_NAMES:
            continue
        nt1, nt2 = nts.get(pair.get("nt1")), nts.get(pair.get("nt2"))
        if not nt1 or not nt2:
            continue
        if nt1.get("chain_name") != wanted or nt2.get("chain_name") != wanted:
            continue
        i, j = int(nt1["index_chain"]) - 1, int(nt2["index_chain"]) - 1
        if 0 <= i < len(sequence) and 0 <= j < len(sequence) and i != j:
            pairs.append((min(i, j), max(i, j)))
    return sequence, sorted(set(pairs))
