"""Cluster configuration for the SHS evaluation pipeline: AlphaFold 3 and DSSR are treated as
drop-in commands rendered from templates, so no site-specific paths are baked into the code."""

import json
import logging
import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class ConfigError(RuntimeError):
    pass


@dataclass
class CommandSpec:
    """A shell command template plus how long to allow it."""
    command: str
    timeout_s: int = 3600

    def render(self, **fields: Any) -> str:
        try:
            return self.command.format(**fields)
        except KeyError as exc:
            raise ConfigError(
                f"Command template references {exc} but the pipeline only provides "
                f"{sorted(fields)}. Fix the template in your config file."
            ) from exc


@dataclass
class SchedulerSpec:
    kind: str = "none"                       # none | slurm
    sbatch_args: List[str] = field(default_factory=list)
    array_throttle: int = 16                 # %N concurrency cap on the job array
    job_name: str = "shs-af3"
    # Lines run on the compute node before the AF3 command: `module load`, `conda activate`,
    # `export`, anything the node does not inherit from the submitting shell.
    setup: List[str] = field(default_factory=list)
    # Slurm's MaxArraySize rejects any array larger than this; work is split across several
    # submissions when it exceeds the limit. `scontrol show config | grep MaxArraySize` to check.
    max_array_size: int = 1000
    # Folds per array task. The binding limit is usually not MaxArraySize but the QOS cap on
    # submitted jobs (`sacctmgr show qos format=Name,MaxSubmitJobsPerUser`), which counts every
    # array task: 9190 one-fold tasks against a 1500 cap is rejected outright, while 10 folds per
    # task is 919 and fits. Walltime must cover the whole batch, not one fold.
    jobs_per_task: int = 1


@dataclass
class Config:
    af3: CommandSpec
    dssr: CommandSpec
    scheduler: SchedulerSpec = field(default_factory=SchedulerSpec)
    af3_seeds: List[int] = field(default_factory=lambda: [1])
    python: str = "python"

    @classmethod
    def load(cls, path: Path) -> "Config":
        text = Path(path).read_text()
        if path.suffix in (".yaml", ".yml"):
            data = _load_yaml(text, path)
        else:
            data = json.loads(text)
        for key in ("af3", "dssr"):
            if key not in data:
                raise ConfigError(f"Config {path} is missing the '{key}' section.")
        return cls(
            af3=CommandSpec(**data["af3"]),
            dssr=CommandSpec(**data["dssr"]),
            scheduler=SchedulerSpec(**data.get("scheduler", {})),
            af3_seeds=data.get("af3_seeds", [1]),
            python=data.get("python", "python"),
        )


def _load_yaml(text: str, path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise ConfigError(
            f"{path} is YAML but PyYAML is not installed. Either 'pip install pyyaml' or "
            f"write the config as JSON instead."
        ) from exc
    return yaml.safe_load(text)


def run_command(command: str, timeout_s: int, dry_run: bool = False,
                log_path: Optional[Path] = None, cwd: Optional[Path] = None) -> int:
    """Run one rendered command. Returns its exit code; -1 marks a timeout."""
    if dry_run:
        print(f"  [dry-run] {command}")
        return 0
    logging.debug("exec: %s", command)
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True,
                                timeout=timeout_s, cwd=str(cwd) if cwd else None)
    except subprocess.TimeoutExpired:
        logging.error("timed out after %ss: %s", timeout_s, command)
        if log_path:
            log_path.write_text(f"TIMEOUT after {timeout_s}s\n{command}\n")
        return -1
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(f"$ {command}\n\n--- stdout ---\n{result.stdout}\n"
                            f"--- stderr ---\n{result.stderr}\n")
    if result.returncode != 0:
        logging.error("exit %d: %s", result.returncode, command)
        logging.error("stderr tail: %s", result.stderr[-800:])
    return result.returncode


def submit_slurm_arrays(config: Config, task_file: Path, n_tasks: int, work_dir: Path,
                        dry_run: bool = False, cwd: Optional[Path] = None) -> List[str]:
    """Submit `n_tasks` commands from `task_file` as one or more SLURM arrays.

    Arrays are split at MaxArraySize; each chunk carries the offset that maps its task ids back
    onto line numbers in the shared task file.
    """
    log_dir = work_dir / "slurm_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    directives = "\n".join(f"#SBATCH {arg}" for arg in _pair_up(config.scheduler.sbatch_args))
    # Compute nodes inherit neither the submitting shell's directory nor its environment, so
    # anchor both explicitly.
    change_dir = f"cd {shlex.quote(str(cwd))}\n" if cwd else ""
    setup = "".join(f"{line}\n" for line in config.scheduler.setup)

    batch = max(1, config.scheduler.jobs_per_task)
    # Array tasks needed once each one runs `batch` folds.
    n_array_tasks = (n_tasks + batch - 1) // batch
    chunk = max(1, config.scheduler.max_array_size)
    job_ids: List[str] = []
    for part, start in enumerate(range(0, n_array_tasks, chunk)):
        size = min(chunk, n_array_tasks - start)
        suffix = "" if n_array_tasks <= chunk else f"_{part}"
        script_path = work_dir / f"af3_array{suffix}.sbatch"
        script = (
            "#!/usr/bin/env bash\n"
            f"#SBATCH --job-name={config.scheduler.job_name}{suffix}\n"
            f"#SBATCH --array=0-{size - 1}%{config.scheduler.array_throttle}\n"
            f"#SBATCH --output={log_dir}/%A_%a.out\n"
            f"#SBATCH --error={log_dir}/%A_%a.err\n"
            f"{directives}\n"
            # -u is enabled only after setup: `conda activate` and some module scripts read
            # unset variables and would abort the task under nounset.
            "set -eo pipefail\n\n"
            f"{change_dir}{setup}"
            "set -u\n\n"
            # One array task runs `batch` consecutive lines of the task file. -e is dropped for
            # the loop so a single failed fold does not discard the rest of the batch; failures
            # are counted and reported in the task's exit status instead.
            f'BATCH={batch}\n'
            f'FIRST=$(( (SLURM_ARRAY_TASK_ID + {start}) * BATCH + 1 ))\n'
            'LAST=$(( FIRST + BATCH - 1 ))\n'
            f'TOTAL={n_tasks}\n'
            'if [ "$LAST" -gt "$TOTAL" ]; then LAST=$TOTAL; fi\n'
            'set +e\n'
            'failed=0\n'
            'for LINE in $(seq "$FIRST" "$LAST"); do\n'
            f'  CMD=$(sed -n "${{LINE}}p" {task_file})\n'
            '  if [ -z "$CMD" ]; then continue; fi\n'
            '  echo "--- line $LINE / $TOTAL ---"\n'
            '  echo "$CMD"\n'
            '  eval "$CMD" || { echo "FOLD FAILED at line $LINE"; failed=$((failed + 1)); }\n'
            'done\n'
            'echo "batch done: lines $FIRST-$LAST, $failed failed"\n'
            'exit $(( failed > 0 ))\n'
        )
        script_path.write_text(script)
        script_path.chmod(0o755)
        logging.info("Wrote array script: %d array tasks x %d folds each (task-file lines %d-%d) -> %s",
                     size, batch, start * batch + 1, min((start + size) * batch, n_tasks), script_path)
        if dry_run:
            print(f"  [dry-run] sbatch {script_path}")
            continue
        result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True)
        if result.returncode != 0:
            raise ConfigError(f"sbatch failed: {result.stderr}")
        job_id = result.stdout.strip().split()[-1]
        logging.info("Submitted SLURM array %s (%d tasks)", job_id, size)
        job_ids.append(job_id)
    return job_ids


def _pair_up(args: List[str]) -> List[str]:
    """Turn ['-p', 'gpu', '--gres=gpu:1'] into ['-p gpu', '--gres=gpu:1'] for #SBATCH lines."""
    out, i = [], 0
    while i < len(args):
        if args[i].startswith("-") and i + 1 < len(args) and not args[i + 1].startswith("-"):
            out.append(f"{args[i]} {args[i + 1]}")
            i += 2
        else:
            out.append(args[i])
            i += 1
    return out
