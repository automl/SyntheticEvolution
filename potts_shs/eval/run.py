"""Four-stage evaluation driver: plan the SHS inputs, fold them with AlphaFold 3, annotate the
predictions with DSSR, and score the recovered base pairs. Every stage is resumable."""

import argparse
import csv
import json
import logging
import sys
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT_DIR = Path(__file__).resolve().parents[2]
GENERATOR_DIR = ROOT_DIR / "SHS-Generator"
POTTS_DIR = ROOT_DIR / "potts_shs"
for path in (str(ROOT_DIR), str(GENERATOR_DIR), str(POTTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from .config import Config, ConfigError, run_command, submit_slurm_arrays
from .dssr import (find_model_cif, find_summary_confidences, parse_pairs, read_confidences,
                   resolve_job_dir, run_dssr)
from .metrics import canonical_fraction, filter_canonical, pair_metrics, summarise

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Returned by the fold stage when work was handed to a scheduler rather than run to completion.
SUBMITTED = 100


@dataclass
class Job:
    job_id: str
    target: str            # 4-character PDB id, used to look up ground truth
    target_key: str        # unique per pickle row; `target` plus a _cN suffix when chains collide
    approach: str
    shs_seed: int
    af3_seed: int
    sequence: str
    seeded_pairs: List[Tuple[int, int]]
    input_json: str
    # AF3 chain ids carrying this row's RNA. Scoring reads the first of these that DSSR emitted,
    # because "the longest chain" is the wrong chain whenever a complex holds two different RNAs.
    chain_ids: List[str] = field(default_factory=list)

    @property
    def af3_dir_name(self) -> str:
        # AF3 lowercases the job name when it creates the output directory.
        return self.job_id.lower()


# -- layout ----------------------------------------------------------------

class Workspace:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.inputs = root / "inputs"
        self.af3 = root / "af3"
        self.dssr = root / "dssr"
        self.logs = root / "logs"
        # DSSR litters its working directory with fixed-name aux files, so each job gets its own.
        self.scratch = root / "dssr_scratch"
        self.dssr_gt = root / "dssr_gt"
        self.manifest = root / "manifest.jsonl"
        for directory in (self.inputs, self.af3, self.dssr, self.logs, self.scratch,
                          self.dssr_gt):
            directory.mkdir(parents=True, exist_ok=True)

    def af3_dir(self, job: Job) -> Path:
        """Where AF3 is told to write. Use found_af3_dir() when reading results back."""
        return self.af3 / job.af3_dir_name

    def found_af3_dir(self, job: Job) -> Optional[Path]:
        return resolve_job_dir(self.af3, job.job_id)

    def dssr_json(self, job: Job) -> Path:
        return self.dssr / f"{job.job_id}.json"

    def read_manifest(self) -> List[Job]:
        if not self.manifest.exists():
            raise SystemExit(f"No manifest at {self.manifest}. Run the 'plan' stage first.")
        known = {f.name for f in fields(Job)}
        jobs = []
        for line in self.manifest.read_text().splitlines():
            if line.strip():
                data = {k: v for k, v in json.loads(line).items() if k in known}
                data["seeded_pairs"] = [tuple(p) for p in data["seeded_pairs"]]
                jobs.append(Job(**data))
        return jobs


# -- targets ---------------------------------------------------------------

def _find_source_json(input_json_dir: Path, target: str) -> Optional[Path]:
    for candidate in (f"{target.lower()}_s1_data.json", f"{target.upper()}_data.json",
                      f"{target.lower()}_data.json"):
        path = input_json_dir / candidate
        if path.is_file():
            return path
    return None


def load_targets(pickle_path: Optional[Path], input_json_dir: Optional[Path],
                 limit: Optional[int], keep_ids: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Targets come from the RNAformer pair mapping, optionally paired with the original AF3
    input JSON so that protein chains and modifications survive into the run.

    One pickle row is one RNA chain, so several rows share a PDB id: 919 rows collapse onto 742
    four-character ids. Each row therefore gets its own `key` (the PDB id plus a _cN suffix from
    the second row of a given id onwards) and, when an input JSON is reused, the ids of the RNA
    chains in that JSON whose sequence this row actually is. Without both, colliding rows
    overwrite each other's inputs and a two-RNA complex has the wrong chain's pairs applied.
    """
    if pickle_path is None:
        raise SystemExit("--targets-pickle is required to plan jobs.")
    import pandas as pd

    frame = pd.read_pickle(pickle_path)
    for column in ("pdb_id", "sequence", "pairs"):
        if column not in frame.columns:
            raise SystemExit(f"{pickle_path} has no '{column}' column; got {list(frame.columns)}")

    wanted: Optional[set] = None
    if keep_ids is not None:
        wanted = {line.strip().upper() for line in keep_ids.read_text().splitlines() if line.strip()}
        logging.info("Restricting to the %d target ids in %s", len(wanted), keep_ids)

    targets: List[Dict[str, Any]] = []
    seen: Dict[str, int] = {}
    unmatched = 0
    # A per-chain source table lists a homodimer's identical copies as separate rows. They would
    # fold to byte-identical jobs, and one row's SHS already goes into every chain sharing its
    # sequence, so keep only the first of each (target, sequence).
    seen_sequences: set = set()
    duplicates = 0
    for _, row in frame.iterrows():
        raw_id = str(row["pdb_id"])
        target = raw_id.split(":")[0].split("_")[0].upper()[:4]
        if wanted is not None and target not in wanted:
            continue
        sequence = str(row["sequence"]).upper()
        if (target, sequence) in seen_sequences:
            duplicates += 1
            continue
        seen_sequences.add((target, sequence))
        pairs = [(int(p[0]), int(p[1])) for p in row["pairs"]]

        seen[target] = seen.get(target, 0) + 1
        key = target if seen[target] == 1 else f"{target}_c{seen[target]}"

        source = _find_source_json(input_json_dir, target) if input_json_dir else None
        chain_ids: List[str] = []
        if source is not None:
            try:
                payload = json.loads(source.read_text())
            except Exception as exc:
                logging.warning("%s: cannot read %s (%s); falling back to an RNA-only input",
                                key, source, exc)
                source = None
            else:
                chain_ids = [c["rna"].get("id") for c in payload.get("sequences", [])
                             if "rna" in c and c["rna"].get("sequence", "").upper() == sequence]
                # Flatten AF3's list-valued ids (one entry per copy of the same chain).
                flat: List[str] = []
                for cid in chain_ids:
                    flat.extend(cid if isinstance(cid, list) else [cid])
                chain_ids = [c for c in flat if c]
                if not chain_ids:
                    # The JSON holds a different RNA than this row. Folding it would apply this
                    # row's pairs to the wrong sequence, so use an RNA-only input instead.
                    unmatched += 1
                    logging.debug("%s: no RNA chain in %s matches this row's %d-nt sequence; "
                                  "using an RNA-only input", key, source.name, len(sequence))
                    source = None

        targets.append({"target": target, "key": key, "sequence": sequence, "pairs": pairs,
                        "source_json": source, "chain_ids": chain_ids or ["A"]})
        if limit and len(targets) >= limit:
            break

    collisions = sum(1 for t in targets if "_c" in t["key"])
    logging.info("Loaded %d target rows over %d PDB ids (%d duplicate chains dropped, %d share a "
                 "PDB id with a different sequence, %d had no matching RNA chain in their input "
                 "JSON)", len(targets), len({t["target"] for t in targets}), duplicates,
                 collisions, unmatched)
    return targets


def generate_msa(target: Dict[str, Any], approach: str, shs_seed: int,
                 extra_flags: List[str]) -> str:
    """Run the real generator on this row's sequence and return the SHS alignment.

    Always driven RNA-only from the pickle sequence: the alignment depends on the sequence, the
    pairs and the seed, not on which complex the chain came from, so the result is identical to
    what the CLI produces on the full input JSON — and it is generated once instead of once per
    AF3 seed.
    """
    import shs_generator

    interactions = json.dumps([[i, j, 1.0] for i, j in target["pairs"]])
    argv = ["shs_generator.py",
            "--interactions", interactions,
            "--pair-mutation-approach", approach,
            "--seed", str(shs_seed),
            "--pdb_id", target["target"],
            "--rna-seq", target["sequence"],
            "--protein-seq", "M"] + extra_flags

    saved, sys.argv = sys.argv, argv
    try:
        args = shs_generator.parse_args()
    finally:
        sys.argv = saved

    payload = shs_generator.MsaGenerator(args).process(write=False)
    if payload is None:
        raise RuntimeError(f"generator returned nothing for {target['key']}")
    for chain in payload["sequences"]:
        if "rna" in chain:
            return chain["rna"]["unpairedMsa"]
    raise RuntimeError(f"generator produced no RNA chain for {target['key']}")


def build_shs_json(target: Dict[str, Any], msa: str, af3_seed: int,
                   job_id: str) -> Dict[str, Any]:
    """Put the generated alignment into an AF3 input for this job.

    With a source JSON the SHS goes into exactly the RNA chains this row is, leaving the rest of
    the complex (proteins, other RNAs, ligands, modifications) untouched. Without one the target
    is folded as a single RNA chain.
    """
    if target["source_json"] is None:
        payload = {
            "dialect": "alphafold3", "version": 1,
            "sequences": [{"rna": {"id": "A", "sequence": target["sequence"],
                                   "modifications": [], "unpairedMsa": msa}}],
        }
    else:
        payload = json.loads(Path(target["source_json"]).read_text())
        wanted = set(target["chain_ids"])
        patched = 0
        for chain in payload.get("sequences", []):
            if "rna" not in chain:
                continue
            cid = chain["rna"].get("id")
            ids = set(cid if isinstance(cid, list) else [cid])
            if ids & wanted:
                chain["rna"]["unpairedMsa"] = msa
                # An MSA path would otherwise win over the inline alignment we just wrote.
                chain["rna"].pop("unpairedMsaPath", None)
                patched += 1
        if not patched:
            raise RuntimeError(f"no RNA chain of {target['source_json']} matched "
                               f"{sorted(wanted)} for {job_id}")

    payload["name"] = job_id
    payload["modelSeeds"] = [af3_seed]
    return payload


# -- stages ----------------------------------------------------------------

def stage_plan(ws: Workspace, config: Config, args) -> int:
    if args.arm_label and len(args.approaches) != 1:
        raise SystemExit("--arm-label names one arm, so pass exactly one --approaches entry.")
    targets = load_targets(args.targets_pickle, args.input_json_dir, args.limit,
                           args.target_ids)
    logging.info("Planning over %d targets x %d approaches x %d SHS seeds x %d AF3 seeds",
                 len(targets), len(args.approaches), len(args.shs_seeds), len(config.af3_seeds))

    jobs: List[Job] = []
    failures = 0
    for target in targets:
        for approach in args.approaches:
            arm = args.arm_label or approach
            for shs_seed in args.shs_seeds:
                job_ids = [f"{target['key']}_{arm}_s{shs_seed}_a{af3_seed}"
                           for af3_seed in config.af3_seeds]
                paths = [ws.inputs / f"{job_id}.json" for job_id in job_ids]

                # The alignment is the same for every AF3 seed; only modelSeeds differs.
                msa = None
                if args.force or not all(p.exists() for p in paths):
                    try:
                        msa = generate_msa(target, approach, shs_seed, args.generator_flags)
                    except Exception as exc:
                        logging.error("plan failed for %s/%s/s%d: %s",
                                      target["key"], approach, shs_seed, exc)
                        failures += 1
                        continue

                for af3_seed, job_id, out_path in zip(config.af3_seeds, job_ids, paths):
                    if msa is not None and (args.force or not out_path.exists()):
                        try:
                            payload = build_shs_json(target, msa, af3_seed, job_id)
                        except Exception as exc:
                            logging.error("plan failed for %s: %s", job_id, exc)
                            failures += 1
                            continue
                        out_path.write_text(json.dumps(payload, indent=2))
                    jobs.append(Job(job_id=job_id, target=target["target"],
                                    target_key=target["key"], approach=arm,
                                    shs_seed=shs_seed, af3_seed=af3_seed,
                                    sequence=target["sequence"],
                                    seeded_pairs=target["pairs"], input_json=str(out_path),
                                    chain_ids=target["chain_ids"]))

    # Merge rather than overwrite: a multi-arm run plans one arm per pass into the same work-dir,
    # and each pass must keep the arms already there. Jobs this pass re-planned win.
    kept = []
    if ws.manifest.exists():
        fresh = {job.job_id for job in jobs}
        kept = [job for job in ws.read_manifest() if job.job_id not in fresh]
    with ws.manifest.open("w") as fh:
        for job in kept + jobs:
            fh.write(json.dumps(asdict(job)) + "\n")
    arms = sorted({job.approach for job in kept + jobs})
    logging.info("Planned %d jobs this pass (%d failed); manifest now holds %d jobs over arms %s",
                 len(jobs), failures, len(kept) + len(jobs), arms)
    return 1 if failures and not jobs else 0


def select_arms(jobs: List[Job], arms: Optional[List[str]]) -> List[Job]:
    """Narrow a manifest to the requested arms, so parallel per-arm jobs share one work-dir."""
    if not arms:
        return jobs
    wanted = set(arms)
    chosen = [job for job in jobs if job.approach in wanted]
    missing = wanted - {job.approach for job in jobs}
    if missing:
        raise SystemExit(f"No jobs for arm(s) {sorted(missing)}; manifest holds "
                         f"{sorted({job.approach for job in jobs})}")
    logging.info("Restricted to arms %s: %d of %d jobs", sorted(wanted), len(chosen), len(jobs))
    return chosen


def pending_fold(ws: Workspace, jobs: List[Job], force: bool) -> List[Job]:
    return [j for j in jobs if force or find_model_cif(ws.found_af3_dir(j)) is None]


def stage_fold(ws: Workspace, config: Config, args) -> int:
    jobs = select_arms(ws.read_manifest(), args.arms)
    todo = pending_fold(ws, jobs, args.force)
    logging.info("%d/%d jobs need folding", len(todo), len(jobs))
    if not todo:
        return 0

    if config.scheduler.kind == "slurm":
        task_file = ws.root / "af3_tasks.txt"
        lines = []
        for job in todo:
            out_dir = ws.af3_dir(job)
            out_dir.mkdir(parents=True, exist_ok=True)
            lines.append(config.af3.render(input_json=job.input_json, input_dir=str(ws.inputs),
                                           output_dir=str(ws.af3), job=job.job_id,
                                           job_output_dir=str(out_dir), seed=job.af3_seed,
                                           root=str(ROOT_DIR)))
        task_file.write_text("\n".join(lines) + "\n")
        submit_slurm_arrays(config, task_file, len(lines), ws.root,
                            dry_run=args.dry_run, cwd=ROOT_DIR)
        logging.info("Task list: %s", task_file)
        # Submission returns immediately, so downstream stages would run against results that do
        # not exist yet. SUBMITTED tells the 'all' driver to stop here. A dry run submitted
        # nothing, so there is no array to wait for.
        return 0 if args.dry_run else SUBMITTED

    failures = 0
    for index, job in enumerate(todo, 1):
        out_dir = ws.af3_dir(job)
        out_dir.mkdir(parents=True, exist_ok=True)
        command = config.af3.render(input_json=str(Path(job.input_json).resolve()),
                                    input_dir=str(ws.inputs),
                                    output_dir=str(ws.af3), job=job.job_id,
                                    job_output_dir=str(out_dir), seed=job.af3_seed,
                                    root=str(ROOT_DIR))
        logging.info("[%d/%d] folding %s", index, len(todo), job.job_id)
        code = run_command(command, config.af3.timeout_s, dry_run=args.dry_run,
                           log_path=ws.logs / f"af3_{job.job_id}.log")
        if code != 0:
            failures += 1
    succeeded = len(todo) - failures
    logging.info("Folding done: %d succeeded, %d failed", succeeded, failures)
    if failures and not succeeded:
        # Every fold failed, which is a misconfiguration rather than a few bad targets. Returning
        # 0 here would let the batch job exit successfully having produced nothing.
        logging.error("Every fold failed. The first log under %s has the reason.", ws.logs)
        return 1
    if failures:
        logging.warning("%d folds failed; rerun this stage to retry just those.", failures)
    return 0


def annotate_ground_truth(ws: Workspace, config: Config, args, jobs: List[Job]) -> None:
    """Run DSSR over the deposited structures once, so the accuracy columns have a reference.

    The repository ships ground truth as mmCIF (evaluation/predictions/gt/<PDB>_gt.cif), not as
    DSSR JSON, so without this step --gt-dssr-dir has nothing to point at.
    """
    targets = sorted({job.target for job in jobs})
    done = missing = failures = 0
    for target in targets:
        out_json = ws.dssr_gt / f"{target.upper()}_gt_dssr.json"
        if out_json.exists() and not args.force:
            done += 1
            continue
        cif = next((p for p in (args.gt_cif_dir / f"{target.upper()}_gt.cif",
                                args.gt_cif_dir / f"{target.lower()}_gt.cif",
                                args.gt_cif_dir / f"{target.upper()}.cif") if p.is_file()), None)
        if cif is None:
            missing += 1
            continue
        ok = run_dssr(cif, out_json, config.dssr, dry_run=args.dry_run,
                      log_path=ws.logs / f"dssr_gt_{target}.log",
                      scratch=ws.scratch / f"gt_{target}", root=ROOT_DIR)
        done += bool(ok)
        failures += (not ok)
    logging.info("Ground truth DSSR: %d annotated, %d failed, %d without a CIF in %s",
                 done, failures, missing, args.gt_cif_dir)


def stage_annotate(ws: Workspace, config: Config, args) -> int:
    jobs = select_arms(ws.read_manifest(), args.arms)
    if args.gt_cif_dir:
        annotate_ground_truth(ws, config, args, jobs)
    done = missing = failures = 0
    for job in jobs:
        out_json = ws.dssr_json(job)
        if out_json.exists() and not args.force:
            done += 1
            continue
        cif = find_model_cif(ws.found_af3_dir(job))
        if cif is None:
            missing += 1
            continue
        ok = run_dssr(cif, out_json, config.dssr, dry_run=args.dry_run,
                      log_path=ws.logs / f"dssr_{job.job_id}.log",
                      scratch=ws.scratch / job.job_id, root=ROOT_DIR)
        done += bool(ok)
        failures += (not ok)
    logging.info("DSSR: %d annotated, %d failed, %d still missing a model CIF",
                 done, failures, missing)
    return 0


def find_ground_truth(gt_dir: Optional[Path], target: str) -> Optional[Path]:
    if gt_dir is None:
        return None
    for candidate in (f"{target.upper()}_gt_dssr.json", f"{target.lower()}_gt_dssr.json",
                      f"{target.upper()}.json", f"{target.lower()}.json"):
        path = gt_dir / candidate
        if path.is_file():
            return path
    return None


def align_reference(reference: str, query: str,
                    pairs: List[Tuple[int, int]]) -> Optional[List[Tuple[int, int]]]:
    """Put reference pair indices into the query's coordinate frame.

    A deposited structure rarely covers exactly the modelled construct -- an extra 5' G, a
    modified residue, a disordered terminus -- so its chain is typically the query plus or minus
    a residue or two. Requiring equal lengths throws away most of the benchmark, while a single
    offset recovers it: 87 of the 91 RNA-only targets are an exact substring match one way or the
    other. Pairs falling outside the query after the shift are dropped, and None means the two
    differ internally and no offset can align them.
    """
    if reference == query:
        return list(pairs)
    if query in reference:
        shift = -reference.index(query)
    elif reference in query:
        shift = query.index(reference)
    else:
        return None
    out = []
    for i, j in pairs:
        a, b = i + shift, j + shift
        if 0 <= a < len(query) and 0 <= b < len(query):
            out.append((min(a, b), max(a, b)))
    return sorted(set(out))


def load_ground_truth(gt_dir: Optional[Path], job: Job,
                      canonical_only: bool) -> Optional[List[Tuple[int, int]]]:
    """Ground-truth pairs in the query's coordinate frame, or None if unusable."""
    path = find_ground_truth(gt_dir, job.target)
    if path is None:
        return None
    try:
        # No `prefer` here: a deposited structure carries the depositor's chain labelling, which
        # need not agree with the ids we assigned in the AF3 input. The sequence identifies it.
        sequence, pairs = parse_pairs(path, match_sequence=job.sequence,
                                      canonical_only=canonical_only)
    except Exception as exc:
        logging.debug("ground truth unreadable for %s: %s", job.target, exc)
        return None
    aligned = align_reference(sequence, job.sequence, pairs)
    if aligned is None:
        logging.debug("%s: ground truth (%d nt) does not align to the query (%d nt)",
                      job.target, len(sequence), len(job.sequence))
    return aligned


def stage_score(ws: Workspace, config: Config, args) -> int:
    jobs = select_arms(ws.read_manifest(), args.arms)
    rows: List[Dict[str, Any]] = []
    skipped = mismatched = 0

    for job in jobs:
        dssr_path = ws.dssr_json(job)
        if not dssr_path.exists():
            skipped += 1
            continue
        try:
            sequence, predicted = parse_pairs(dssr_path, chain=args.chain,
                                              prefer=job.chain_ids,
                                              match_sequence=job.sequence,
                                              canonical_only=args.canonical_only)
        except Exception as exc:
            logging.warning("could not parse %s: %s", dssr_path, exc)
            skipped += 1
            continue

        # A length mismatch means DSSR indices do not line up with the seeded pairs, which would
        # silently produce meaningless metrics.
        if len(sequence) != len(job.sequence):
            logging.warning("%s: DSSR chain length %d != query length %d; skipping",
                            job.job_id, len(sequence), len(job.sequence))
            mismatched += 1
            continue

        row: Dict[str, Any] = {
            "job_id": job.job_id, "target": job.target, "target_key": job.target_key,
            "approach": job.approach, "shs_seed": job.shs_seed, "af3_seed": job.af3_seed,
            "length": len(job.sequence),
        }
        # The reference is filtered the same way as the prediction, so a canonical-only run does
        # not simply penalise recall on non-canonical reference pairs.
        seeded_pairs = (filter_canonical(job.sequence, job.seeded_pairs)
                        if args.canonical_only else job.seeded_pairs)
        closure = pair_metrics(predicted, seeded_pairs, len(job.sequence))
        row.update({f"closure_{k}": v for k, v in closure.items()})
        row["canonical_fraction"] = canonical_fraction(sequence, predicted)
        found = ws.found_af3_dir(job)
        row.update(read_confidences(find_summary_confidences(found) if found else None))

        truth = load_ground_truth(args.gt_dssr_dir, job, args.canonical_only)
        if truth is not None:
            accuracy = pair_metrics(predicted, truth, len(job.sequence))
            row.update({f"accuracy_{k}": v for k, v in accuracy.items()})
            seeded = pair_metrics(seeded_pairs, truth, len(job.sequence))
            row["input_f1"] = seeded["f1"]
        rows.append(row)

    if not rows:
        logging.error("Nothing to score: %d skipped, %d length-mismatched.", skipped, mismatched)
        return 1

    out_csv = ws.root / "scores.csv"
    fields = sorted({k for r in rows for k in r})
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %d scored jobs to %s (%d skipped, %d length-mismatched)",
                 len(rows), out_csv, skipped, mismatched)

    keys = ["closure_f1", "closure_mcc", "closure_precision", "closure_recall",
            "canonical_fraction", "accuracy_f1", "accuracy_mcc", "input_f1", "ranking_score"]
    by_approach: Dict[str, List[Dict[str, float]]] = {}
    for row in rows:
        by_approach.setdefault(row["approach"], []).append(row)

    def label(key: str) -> str:
        short = key.replace("closure_", "S.").replace("accuracy_", "A.")
        short = short.replace("canonical_fraction", "canonical").replace("ranking_score", "rank")
        return short[:13]

    print(f"\n{'approach':<20}{'n':>5}" + "".join(f"{label(k):>14}" for k in keys))
    print("-" * (25 + 14 * len(keys)))
    for approach, group in sorted(by_approach.items()):
        stats = summarise(group, keys)
        cells = "".join(f"{stats[k]:>14.4f}" if stats[k] == stats[k] else f"{'-':>14}" for k in keys)
        print(f"{approach:<20}{len(group):>5}{cells}")
    print("\nS.* = loop closure against the seeded pairs; A.* = accuracy against ground-truth DSSR.")
    return 0


# -- CLI -------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["plan", "fold", "annotate", "score", "all"])
    parser.add_argument("--config", type=Path, required=True, help="Cluster config (YAML or JSON)")
    parser.add_argument("--work-dir", type=Path, default=Path("potts_shs/runs/current"))
    parser.add_argument("--targets-pickle", type=Path,
                        default=ROOT_DIR / "synthetic_msa_selection_max_rna_len_200_min_rna_len_15"
                                           "_max_protein_len_200_rnaformer2025_pdb_id2pairs_mapping.pkl")
    parser.add_argument("--input-json-dir", type=Path,
                        help="Original AF3 data files, so protein chains survive into the run")
    parser.add_argument("--gt-cif-dir", type=Path,
                        help="Deposited structures to annotate once during 'annotate', so the "
                             "accuracy columns have a reference. In this repository: "
                             "evaluation/predictions/gt")
    parser.add_argument("--gt-dssr-dir", type=Path,
                        help="Ready-made ground-truth DSSR annotations. Defaults to whatever "
                             "--gt-cif-dir produced under <work-dir>/dssr_gt")
    parser.add_argument("--approaches", nargs="+", default=["potts", "watson_crick_cov"])
    parser.add_argument("--arm-label",
                        help="Name this arm in job ids and the 'approach' column instead of the "
                             "generator's approach name. Lets the same approach appear twice under "
                             "different --generator-flags, e.g. potts at two couplings. Only valid "
                             "with a single --approaches entry; plan merges into any existing "
                             "manifest, so run one pass per arm.")
    parser.add_argument("--shs-seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--generator-flags", nargs=argparse.REMAINDER, default=[],
                        help="Everything after this flag is passed through to shs_generator.py")
    parser.add_argument("--arms", nargs="+",
                        help="Restrict fold/annotate/score to these arms of the manifest. Lets one "
                             "job per arm fold in parallel against a shared work-dir, without "
                             "splitting the manifest or the ground truth.")
    parser.add_argument("--chain",
                        help="Force a DSSR chain for every job; by default the chain the target "
                             "row actually is, falling back to the longest")
    parser.add_argument("--canonical-only", action="store_true",
                        help="Score only Watson-Crick and wobble pairs, on both the predicted and "
                             "the reference side. Matches the protocol behind the existing AF3 and "
                             "RhoFold MCC numbers; leave it off and those are not comparable")
    parser.add_argument("--target-ids", type=Path,
                        help="File of PDB ids, one per line; only those targets are planned. Use "
                             "it to scope a run to a curated subset, e.g. the RNA-only entries.")
    parser.add_argument("--limit", type=int, help="Only plan the first N targets")
    parser.add_argument("--force", action="store_true", help="Redo work that already has outputs")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the AF3 and DSSR commands instead of running them")
    args = parser.parse_args(argv)

    # argparse.REMAINDER swallows everything after the flag, so a misplaced --generator-flags
    # quietly eats the pipeline's own options instead of erroring.
    stray = [f for f in args.generator_flags
             if f in {"--config", "--work-dir", "--approaches", "--shs-seeds", "--limit",
                      "--input-json-dir", "--gt-dssr-dir", "--targets-pickle", "--canonical-only"}]
    if stray:
        logging.error("--generator-flags must come last; it captured %s, which are pipeline "
                      "options, not generator options.", stray)
        return 2

    try:
        config = Config.load(args.config)
    except ConfigError as exc:
        logging.error("%s", exc)
        return 2

    # Absolute, because these paths are handed to container --bind arguments, whose destinations
    # must be absolute, and to commands that may not run in this working directory.
    args.work_dir = args.work_dir.resolve()
    ws = Workspace(args.work_dir)
    if args.gt_dssr_dir is None and any(ws.dssr_gt.iterdir()):
        args.gt_dssr_dir = ws.dssr_gt
        logging.info("Using ground truth annotated into %s", ws.dssr_gt)

    stages = {"plan": stage_plan, "fold": stage_fold, "annotate": stage_annotate,
              "score": stage_score}
    order = ["plan", "fold", "annotate", "score"] if args.stage == "all" else [args.stage]
    for name in order:
        logging.info("=== stage: %s ===", name)
        code = stages[name](ws, config, args)
        if code == SUBMITTED:
            logging.info(
                "Folding was submitted to the scheduler and is still running. Wait for the array "
                "to finish, then continue with:\n"
                "  python -m potts_shs.eval.run annotate --config %s --work-dir %s\n"
                "  python -m potts_shs.eval.run score    --config %s --work-dir %s",
                args.config, args.work_dir, args.config, args.work_dir)
            return 0
        if code != 0:
            return code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
