"""Execute one AlphaFold 3 prediction in the active AF3 environment."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def run_af3_task(
    input_path: Path,
    output_directory: Path,
    model_directory: Path,
) -> Path:
    """Run one AF3 prediction, validate its model, and publish a marker."""

    input_path = input_path.resolve()
    output_directory = output_directory.resolve()
    model_directory = model_directory.expanduser().resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"AF3 input does not exist: {input_path}")

    alphafold_bin_directory = os.environ.get("ALPHAFOLD_BIN_DIR")
    database_directory = os.environ.get("ALPHAFOLD_DATABASES")

    if not alphafold_bin_directory:
        raise RuntimeError(
            "ALPHAFOLD_BIN_DIR is not set; was the AF3 module loaded?"
        )

    if not database_directory:
        raise RuntimeError(
            "ALPHAFOLD_DATABASES is not set; was the AF3 module loaded?"
        )

    alphafold_script = (
        Path(alphafold_bin_directory) / "run_alphafold.py"
    )
    if not alphafold_script.is_file():
        raise FileNotFoundError(
            f"AlphaFold entry point does not exist: {alphafold_script}"
        )

    output_directory.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            sys.executable,
            str(alphafold_script),
            f"--json_path={input_path}",
            f"--output_dir={output_directory}",
            f"--model_dir={model_directory}",
            f"--db_dir={database_directory}",
            "--norun_data_pipeline",
        ],
        check=True,
    )

    models = list(output_directory.glob("*/*_model.cif"))
    if len(models) != 1 or models[0].stat().st_size == 0:
        raise RuntimeError(
            "Expected exactly one nonempty AF3 model, "
            f"found {models}"
        )

    model = models[0].resolve()
    marker = output_directory / "af3_success.json"
    temporary = marker.with_suffix(".tmp")
    temporary.write_text(
        json.dumps({"model": str(model)}, indent=2)
    )
    temporary.replace(marker)

    return model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--model-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_af3_task(
        input_path=args.input,
        output_directory=args.output,
        model_directory=args.model_dir,
    )


if __name__ == "__main__":
    main()