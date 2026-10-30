#!/usr/bin/env python3
"""Run all analysis modules against an exported HPO result snapshot."""

import argparse
import importlib
import logging
import shutil
import sys
from pathlib import Path
from types import ModuleType
from typing import Optional

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.results_loader import ResultsContext, load_results


LOGGER = logging.getLogger("analyze_results")
MODULES_DIR = Path(__file__).resolve().parent / "modules"


class _ColorFormatter(logging.Formatter):
    COLORS = {
        logging.INFO: "\033[32m",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[1;31m",
    }
    BOLD = "\033[1m"
    RESET = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        message = super().format(record)
        logger_name = (
            f"{self.BOLD}[ORCHESTRATOR]{self.RESET}"
            if record.name == LOGGER.name
            else f"[{record.name}]"
        )
        if not sys.stderr.isatty():
            return f"{record.levelname} {logger_name} {message}"
        level = record.levelname
        color = self.COLORS.get(record.levelno)
        if color:
            level = f"{color}{level}{self.RESET}"
        return f"{level} {logger_name} {message}"


def _configure_logging() -> None:
    handler = logging.StreamHandler()
    handler.setFormatter(_ColorFormatter("%(message)s"))
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.handlers.clear()
    root_logger.addHandler(handler)


def _confirm_recreate(output_dir: Path) -> None:
    answer = input(
        f"Analysis output already exists at {output_dir}. "
        "Recreate it? [Y/n] "
    ).strip().lower()
    if answer not in {"", "y", "yes"}:
        raise SystemExit("Analysis cancelled.")
    shutil.rmtree(output_dir)


def prepare_output_directory(output_dir: Path) -> None:
    if output_dir.exists():
        if not output_dir.is_dir():
            raise NotADirectoryError(f"Analysis output is not a directory: {output_dir}")
        _confirm_recreate(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def _ensure_separate_directories(results_dir: Path, output_dir: Path) -> None:
    if (
        output_dir == results_dir
        or output_dir in results_dir.parents
        or results_dir in output_dir.parents
    ):
        raise ValueError(
            f"Analysis output must be separate from results: {output_dir}"
        )


def discover_modules() -> list[ModuleType]:
    modules = []
    for path in sorted(MODULES_DIR.glob("*.py")):
        if path.name.startswith("_"):
            continue
        module_name = f"analysis.modules.{path.stem}"
        modules.append(importlib.import_module(module_name))
    return modules


def run_analysis(
    results_dir: Path,
    output_dir: Path,
    analysis_config_path: Optional[Path] = None,
) -> None:
    _ensure_separate_directories(results_dir, output_dir)
    LOGGER.info("Loading exported results from %s", results_dir)
    context = load_results(results_dir, output_dir, analysis_config_path)
    LOGGER.info("Loaded %d trajectory configurations", len(context.config_ids))

    modules = discover_modules()
    if not modules:
        raise RuntimeError(f"No analysis modules found in {MODULES_DIR}")

    failures: list[tuple[str, Exception]] = []
    for module in modules:
        name = module.__name__.rsplit(".", maxsplit=1)[-1]
        module_output = output_dir / name
        LOGGER.info("[START] module: %s", name)
        try:
            run = module.run
            run(context, module_output)
        except Exception as exc:
            LOGGER.exception("[FAIL] module: %s", name)
            failures.append((name, exc))
        else:
            LOGGER.info("[DONE] module: %s", name)

    if failures:
        names = ", ".join(name for name, _ in failures)
        raise RuntimeError(f"Analysis modules failed: {names}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run all analyses for an exported HPO results directory."
    )
    parser.add_argument("results", type=Path, help="Exported results directory")
    parser.add_argument(
        "--output",
        type=Path,
        help="Analysis output directory (default: sibling analysis directory)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="JSON analysis configuration file, used by configurable modules",
    )
    args = parser.parse_args()
    _configure_logging()

    results_dir = args.results.expanduser().resolve()
    output_dir = (
        args.output.expanduser().resolve()
        if args.output is not None
        else results_dir.parent / "analysis"
    )
    LOGGER.info("Preparing analysis output directory: %s", output_dir)
    prepare_output_directory(output_dir)
    analysis_config_path = (
        args.config.expanduser().resolve() if args.config is not None else None
    )
    run_analysis(results_dir, output_dir, analysis_config_path)
    LOGGER.info("All analyses completed successfully.")


if __name__ == "__main__":
    main()
