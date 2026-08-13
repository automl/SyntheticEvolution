"""Four-stage evaluation driver: plan the SHS inputs, fold them with AlphaFold 3, annotate the
predictions with DSSR, and score the recovered base pairs. Every stage is resumable."""

import argparse
import csv
import json
import logging
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

ROOT_DIR = Path(__file__).resolve().parents[1]
GENERATOR_DIR = ROOT_DIR / "SHS-Generator"
for path in (str(ROOT_DIR), str(GENERATOR_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from .config import Config, ConfigError, run_command, submit_slurm_arrays
from .dssr import (find_model_cif, find_summary_confidences, parse_pairs, read_confidences,
                   resolve_job_dir, run_dssr)
from .metrics import canonical_fraction, pair_metrics, summarise

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Returned by the fold stage when work was handed to a scheduler rather than run to completion.
SUBMITTED = 100


@dataclass
class Job:
    job_id: str
    target: str
    approach: str
    shs_seed: int
    af3_seed: int
    sequence: str
    seeded_pairs: List[Tuple[int, int]]
    input_json: str

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
        self.manifest = root / "manifest.jsonl"
        for directory in (self.inputs, self.af3, self.dssr, self.logs):
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
        jobs = []
        for line in self.manifest.read_text().splitlines():
            if line.strip():
                data = json.loads(line)
                data["seeded_pairs"] = [tuple(p) for p in data["seeded_pairs"]]
                jobs.append(Job(**data))
        return jobs


# -- targets ---------------------------------------------------------------

def load_targets(pickle_path: Optional[Path], input_json_dir: Optional[Path],
                 limit: Optional[int]) -> List[Dict[str, Any]]:
    """Targets come from the RNAformer pair mapping, optionally paired with the original AF3
    input JSON so that protein chains and modifications survive into the run."""
    if pickle_path is None:
        raise SystemExit("--targets-pickle is required to plan jobs.")
    import pandas as pd

    frame = pd.read_pickle(pickle_path)
    for column in ("pdb_id", "sequence", "pairs"):
        if column not in frame.columns:
            raise SystemExit(f"{pickle_path} has no '{column}' column; got {list(frame.columns)}")

    targets = []
    for _, row in frame.iterrows():
        raw_id = str(row["pdb_id"])
        target = raw_id.split(":")[0].split("_")[0].upper()[:4]
        pairs = [(int(p[0]), int(p[1])) for p in row["pairs"]]
        source = None
        if input_json_dir:
            for candidate in (f"{target.lower()}_s1_data.json", f"{target.upper()}_data.json",
                              f"{target.lower()}_data.json"):
                path = input_json_dir / candidate
                if path.is_file():
                    source = path
                    break
        targets.append({"target": target, "sequence": str(row["sequence"]).upper(),
                        "pairs": pairs, "source_json": source})
        if limit and len(targets) >= limit:
            break
    return targets


def build_shs_json(target: Dict[str, Any], approach: str, shs_seed: int, af3_seed: int,
                   job_id: str, extra_flags: List[str]) -> Dict[str, Any]:
    """Drive the real generator so the evaluated MSAs are exactly what the CLI would produce."""
    import shs_generator

    interactions = json.dumps([[i, j, 1.0] for i, j in target["pairs"]])
    argv = ["shs_generator.py",
            "--interactions", interactions,
            "--pair-mutation-approach", approach,
            "--seed", str(shs_seed),
            "--pdb_id", target["target"]]
    if target["source_json"] is not None:
        argv += ["--input_json_path", str(target["source_json"])]
    else:
        argv += ["--rna-seq", target["sequence"], "--protein-seq", "M"]
    argv += extra_flags

    saved, sys.argv = sys.argv, argv
    try:
        args = shs_generator.parse_args()
    finally:
        sys.argv = saved

    generator = shs_generator.MsaGenerator(args)
    payload = generator.process(write=False)
    if payload is None:
        raise RuntimeError(f"generator returned nothing for {job_id}")
    if target["source_json"] is None:
        # Drop the placeholder protein so an RNA-only target is folded as RNA only.
        payload["sequences"] = [c for c in payload["sequences"] if "rna" in c]
    payload["name"] = job_id
    payload["modelSeeds"] = [af3_seed]
    return payload


# -- stages ----------------------------------------------------------------

def stage_plan(ws: Workspace, config: Config, args) -> int:
    targets = load_targets(args.targets_pickle, args.input_json_dir, args.limit)
    logging.info("Planning over %d targets x %d approaches x %d SHS seeds x %d AF3 seeds",
                 len(targets), len(args.approaches), len(args.shs_seeds), len(config.af3_seeds))

    jobs: List[Job] = []
    failures = 0
    for target in targets:
        for approach in args.approaches:
            for shs_seed in args.shs_seeds:
                for af3_seed in config.af3_seeds:
                    job_id = f"{target['target']}_{approach}_s{shs_seed}_a{af3_seed}"
                    out_path = ws.inputs / f"{job_id}.json"
                    if not out_path.exists() or args.force:
                        try:
                            payload = build_shs_json(target, approach, shs_seed, af3_seed,
                                                     job_id, args.generator_flags)
                        except Exception as exc:
                            logging.error("plan failed for %s: %s", job_id, exc)
                            failures += 1
                            continue
                        out_path.write_text(json.dumps(payload, indent=2))
                    jobs.append(Job(job_id=job_id, target=target["target"], approach=approach,
                                    shs_seed=shs_seed, af3_seed=af3_seed,
                                    sequence=target["sequence"],
                                    seeded_pairs=target["pairs"], input_json=str(out_path)))

    with ws.manifest.open("w") as fh:
        for job in jobs:
            fh.write(json.dumps(asdict(job)) + "\n")
    logging.info("Wrote %d jobs to %s (%d failed to plan)", len(jobs), ws.manifest, failures)
    return 1 if failures and not jobs else 0


def pending_fold(ws: Workspace, jobs: List[Job], force: bool) -> List[Job]:
    return [j for j in jobs if force or find_model_cif(ws.found_af3_dir(j)) is None]


def stage_fold(ws: Workspace, config: Config, args) -> int:
    jobs = ws.read_manifest()
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
                                           job_output_dir=str(out_dir), seed=job.af3_seed))
        task_file.write_text("\n".join(lines) + "\n")
        submit_slurm_arrays(config, task_file, len(lines), ws.root,
                            dry_run=args.dry_run, cwd=ROOT_DIR)
        logging.info("Task list: %s", task_file)
        # Submission returns immediately, so downstream stages would run against results that do
        # not exist yet. SUBMITTED tells the 'all' driver to stop here.
        return SUBMITTED

    failures = 0
    for index, job in enumerate(todo, 1):
        out_dir = ws.af3_dir(job)
        out_dir.mkdir(parents=True, exist_ok=True)
        command = config.af3.render(input_json=job.input_json, input_dir=str(ws.inputs),
                                    output_dir=str(ws.af3), job=job.job_id,
                                    job_output_dir=str(out_dir), seed=job.af3_seed)
        logging.info("[%d/%d] folding %s", index, len(todo), job.job_id)
        code = run_command(command, config.af3.timeout_s, dry_run=args.dry_run,
                           log_path=ws.logs / f"af3_{job.job_id}.log")
        if code != 0:
            failures += 1
    logging.info("Folding done, %d failures", failures)
    return 0


def stage_annotate(ws: Workspace, config: Config, args) -> int:
    jobs = ws.read_manifest()
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
                      log_path=ws.logs / f"dssr_{job.job_id}.log")
        done += bool(ok)
        failures += (not ok)
    logging.info("DSSR: %d annotated, %d failed, %d still missing a model CIF",
                 done, failures, missing)
    return 0


def load_ground_truth(gt_dir: Optional[Path], target: str) -> Optional[Tuple[str, List[Tuple[int, int]]]]:
    if gt_dir is None:
        return None
    for candidate in (f"{target.upper()}_gt_dssr.json", f"{target.lower()}_gt_dssr.json",
                      f"{target.upper()}.json"):
        path = gt_dir / candidate
        if path.is_file():
            try:
                return parse_pairs(path)
            except Exception as exc:
                logging.debug("ground truth unreadable for %s: %s", target, exc)
                return None
    return None


def stage_score(ws: Workspace, config: Config, args) -> int:
    jobs = ws.read_manifest()
    rows: List[Dict[str, Any]] = []
    skipped = mismatched = 0

    for job in jobs:
        dssr_path = ws.dssr_json(job)
        if not dssr_path.exists():
            skipped += 1
            continue
        try:
            sequence, predicted = parse_pairs(dssr_path, chain=args.chain)
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
            "job_id": job.job_id, "target": job.target, "approach": job.approach,
            "shs_seed": job.shs_seed, "af3_seed": job.af3_seed, "length": len(job.sequence),
        }
        closure = pair_metrics(predicted, job.seeded_pairs, len(job.sequence))
        row.update({f"closure_{k}": v for k, v in closure.items()})
        row["canonical_fraction"] = canonical_fraction(sequence, predicted)
        found = ws.found_af3_dir(job)
        row.update(read_confidences(find_summary_confidences(found) if found else None))

        truth = load_ground_truth(args.gt_dssr_dir, job.target)
        if truth is not None and len(truth[0]) == len(job.sequence):
            accuracy = pair_metrics(predicted, truth[1], len(job.sequence))
            row.update({f"accuracy_{k}": v for k, v in accuracy.items()})
            seeded = pair_metrics(job.seeded_pairs, truth[1], len(job.sequence))
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
    parser.add_argument("--work-dir", type=Path, default=Path("results/shs_eval"))
    parser.add_argument("--targets-pickle", type=Path,
                        default=ROOT_DIR / "synthetic_msa_selection_max_rna_len_200_min_rna_len_15"
                                           "_max_protein_len_200_rnaformer2025_pdb_id2pairs_mapping.pkl")
    parser.add_argument("--input-json-dir", type=Path,
                        help="Original AF3 data files, so protein chains survive into the run")
    parser.add_argument("--gt-dssr-dir", type=Path,
                        help="Ground-truth DSSR annotations, e.g. evaluation/dssr/gt")
    parser.add_argument("--approaches", nargs="+", default=["potts", "watson_crick_cov"])
    parser.add_argument("--shs-seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--generator-flags", nargs=argparse.REMAINDER, default=[],
                        help="Everything after this flag is passed through to shs_generator.py")
    parser.add_argument("--chain", help="DSSR chain to score; defaults to the longest")
    parser.add_argument("--limit", type=int, help="Only plan the first N targets")
    parser.add_argument("--force", action="store_true", help="Redo work that already has outputs")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the AF3 and DSSR commands instead of running them")
    args = parser.parse_args(argv)

    try:
        config = Config.load(args.config)
    except ConfigError as exc:
        logging.error("%s", exc)
        return 2

    ws = Workspace(args.work_dir)
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
                "  python -m shs_eval.run annotate --config %s --work-dir %s\n"
                "  python -m shs_eval.run score    --config %s --work-dir %s",
                args.config, args.work_dir, args.config, args.work_dir)
            return 0
        if code != 0:
            return code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
