#!/usr/bin/env python3

import argparse
import json
import re
import sys
from pathlib import Path


def pairs_to_dotbracket(pairs, size):
    """Convert base pairs to a dot-bracket string."""
    structure = ["."] * size

    for a, b in pairs:
        a, b = sorted((int(a), int(b)))

        if a < 0 or b >= size:
            raise ValueError(
                f"Pair ({a}, {b}) is outside structure length {size}"
            )

        if a == b:
            raise ValueError(f"Invalid self-pair: ({a}, {b})")

        if structure[a] != "." or structure[b] != ".":
            raise ValueError(
                f"Overlapping base pair: ({a}, {b})"
            )

        structure[a] = "("
        structure[b] = ")"

    return "".join(structure)


def infer_structure_length(*pair_lists):
    """Infer sequence length from the largest paired index."""
    indices = [
        int(index)
        for pairs in pair_lists
        for pair in pairs
        for index in pair
    ]

    if not indices:
        raise ValueError("Cannot infer structure length: no base pairs found")

    return max(indices) + 1


def load_trial_result(path):
    """Load and validate a trial_result.json file."""
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)

    if "scores" not in data:
        raise ValueError(f"Missing 'scores' in {path}")

    return data["scores"]


def format_trial_result(path):
    """Format target and predicted structures from one trial result."""
    scores = load_trial_result(path)
    output = []

    for score in scores:
        sequence_id = score["id"]
        target_pairs = score["target_pairs"]
        predicted_pairs = score["predicted_pairs"]

        size = infer_structure_length(
            target_pairs,
            predicted_pairs,
        )

        target = pairs_to_dotbracket(target_pairs, size)
        predicted = pairs_to_dotbracket(predicted_pairs, size)

        output.append(
            f"{sequence_id}\n"
            f"input:       {target}\n"
            f"prediction:  {predicted}\n"
        )

    return "\n".join(output)


def parse_config_ids(trajectory_path):
    """Extract Config IDs in their original order."""
    text = trajectory_path.read_text(encoding="utf-8")

    config_ids = re.findall(
        r"^\s*Config ID:\s*(\d+)\s*$",
        text,
        flags=re.MULTILINE,
    )

    if not config_ids:
        raise ValueError(
            f"No Config IDs found in {trajectory_path}"
        )

    return config_ids


def format_trajectory(neps_dir, trajectory_path):
    """Format structures across all configurations in the trajectory."""
    config_ids = parse_config_ids(trajectory_path)

    # Load all trial results, failing immediately if any is missing.
    results = {}

    for config_id in config_ids:
        trial_path = (
            neps_dir
            / "configs"
            / f"config_{config_id}"
            / "artifacts"
            / "trial_result.json"
        )

        if not trial_path.is_file():
            raise FileNotFoundError(
                f"Missing trial result for Config ID {config_id}: "
                f"{trial_path}"
            )

        results[config_id] = load_trial_result(trial_path)

    # Use the first configuration as the reference for sequence IDs
    # and target structures.
    reference_scores = results[config_ids[0]]

    reference_ids = [score["id"] for score in reference_scores]

    # Validate sequence IDs and order across configurations.
    for config_id in config_ids[1:]:
        ids = [score["id"] for score in results[config_id]]

        if ids != reference_ids:
            raise ValueError(
                f"Sequence IDs or order differ in Config ID {config_id}"
            )

    # Determine label width for aligned output.
    labels = ["input:"] + [
        f"CONFIG_ID-{config_id}:"
        for config_id in config_ids
    ]
    label_width = max(map(len, labels)) + 1

    output = []

    for sequence_index, reference in enumerate(reference_scores):
        sequence_id = reference["id"]
        target_pairs = reference["target_pairs"]

        all_pairs = [target_pairs]

        for config_id in config_ids:
            score = results[config_id][sequence_index]
            all_pairs.append(score["predicted_pairs"])

        size = infer_structure_length(*all_pairs)

        target = pairs_to_dotbracket(target_pairs, size)

        output.append(f"{sequence_id}:")

        output.append(
            f"{'input:':<{label_width}}{target}"
        )

        for config_id in config_ids:
            score = results[config_id][sequence_index]

            predicted = pairs_to_dotbracket(
                score["predicted_pairs"],
                size,
            )

            label = f"CONFIG_ID-{config_id}:"

            output.append(
                f"{label:<{label_width}}{predicted}"
            )

        output.append("")

    return "\n".join(output)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract dot-bracket structures from trial results "
            "or a configuration trajectory."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input JSON or trajectory TXT file",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output directory",
    )

    args = parser.parse_args()

    input_path = args.input.resolve()
    output_dir = args.output.resolve()

    if not input_path.is_file():
        parser.error(f"Input file does not exist: {input_path}")

    try:
        if input_path.name == "trial_result.json":
            content = format_trial_result(input_path)
            output_name = "dotbracket.txt"

        elif input_path.name == "best_config_trajectory.txt":
            # Expected layout:
            # neps/summary/best_config_trajectory.txt
            neps_dir = input_path.parent.parent

            content = format_trajectory(
                neps_dir,
                input_path,
            )
            output_name = "dotbracket_trajectory.txt"

        else:
            parser.error(
                "Input filename must be either "
                "'trial_result.json' or "
                "'best_config_trajectory.txt'"
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / output_name

        output_path.write_text(
            content.rstrip() + "\n",
            encoding="utf-8",
        )

        print(f"Wrote: {output_path}")

    except (ValueError, KeyError, json.JSONDecodeError, OSError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()