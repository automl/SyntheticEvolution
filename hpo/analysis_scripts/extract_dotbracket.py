#!/usr/bin/env python3

import argparse
import json
import logging
import re
import sys
from pathlib import Path


logger = logging.getLogger("extract_dotbracket")


def configure_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def pairs_to_dotbracket(pairs, size):
    """
    Convert base-pair indices to a best-effort dot-bracket string.

    If a nucleotide occurs in multiple pairs, the first pair is kept
    and subsequent conflicting pairs are skipped.

    Crossing pairs are also detected and skipped.

    Returns:
        structure, warnings
    """
    structure = ["."] * size
    paired_with = {}
    accepted_pairs = []
    warnings = []

    for pair in pairs:
        if len(pair) != 2:
            warnings.append(
                f"invalid pair {pair}; expected [a, b]"
            )
            continue

        a, b = sorted(pair)

        if a < 0 or b >= size:
            warnings.append(
                f"pair ({a}, {b}) is outside structure length {size}"
            )
            continue

        if a == b:
            warnings.append(
                f"self-pair ({a}, {b})"
            )
            continue

        # Check whether either nucleotide is already paired.
        conflicts = []

        if a in paired_with:
            conflicts.append(
                f"{a}-{paired_with[a]}"
            )

        if b in paired_with:
            conflicts.append(
                f"{b}-{paired_with[b]}"
            )

        if conflicts:
            warnings.append(
                f"overlapping pair {a}-{b} "
                f"(conflicts with {', '.join(conflicts)})"
            )
            continue

        # Check for crossing pairs.
        #
        # Existing pair (c, d) crosses (a, b) if:
        #
        #     a < c < b < d
        #
        # or:
        #
        #     c < a < d < b
        #
        crossing_pair = None

        for c, d in accepted_pairs:
            if (a < c < b < d) or (c < a < d < b):
                crossing_pair = (c, d)
                break

        if crossing_pair is not None:
            warnings.append(
                f"crossing pair {a}-{b} "
                f"(crosses {crossing_pair[0]}-{crossing_pair[1]})"
            )
            continue

        # Accept the pair.
        paired_with[a] = b
        paired_with[b] = a
        accepted_pairs.append((a, b))

        structure[a] = "("
        structure[b] = ")"

    return "".join(structure), warnings


def infer_structure_length(*pair_lists):
    """
    Infer structure length from the largest nucleotide index.

    Note:
        If the sequence has unpaired nucleotides after the final
        paired nucleotide, they cannot be inferred from pairs alone.
    """
    indices = [
        index
        for pairs in pair_lists
        for pair in pairs
        for index in pair
    ]

    if not indices:
        raise ValueError(
            "Cannot infer structure length: no base pairs found"
        )

    return max(indices) + 1


def load_trial_result(path):
    """
    Load a trial_result.json.

    Returns:
        {
            sequence_id: {
                "target": str,
                "prediction": str,
                "prediction_warnings": list[str],
                "target_warnings": list[str],
            }
        }
    """
    logger.info("Reading trial result: %s", path)

    try:
        with path.open() as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(
            f"Invalid JSON in {path}: {e}"
        ) from e

    if "scores" not in data:
        raise ValueError(
            f"Missing 'scores' field in {path}"
        )

    results = {}

    logger.info(
        "Found %d sequences in %s",
        len(data["scores"]),
        path,
    )

    for index, score in enumerate(data["scores"], start=1):
        sequence_id = score.get("id")

        if sequence_id is None:
            raise ValueError(
                f"{path}: score #{index} has no 'id'"
            )

        logger.debug(
            "Processing sequence %s from %s",
            sequence_id,
            path,
        )

        try:
            target_pairs = score["target_pairs"]
            predicted_pairs = score["predicted_pairs"]
        except KeyError as e:
            raise ValueError(
                f"{path}: sequence {sequence_id} is missing "
                f"{e.args[0]!r}"
            ) from e

        size = infer_structure_length(
            target_pairs,
            predicted_pairs,
        )

        target, target_warnings = pairs_to_dotbracket(
            target_pairs,
            size,
        )

        prediction, prediction_warnings = pairs_to_dotbracket(
            predicted_pairs,
            size,
        )

        for warning in target_warnings:
            logger.warning(
                "%s, target: %s",
                sequence_id,
                warning,
            )

        for warning in prediction_warnings:
            logger.warning(
                "%s, prediction: %s",
                sequence_id,
                warning,
            )

        results[sequence_id] = {
            "target": target,
            "prediction": prediction,
            "prediction_warnings": prediction_warnings,
            "target_warnings": target_warnings,
        }

    return results


def format_prediction(prediction, warnings):
    """Append INVALID(...) to a prediction when necessary."""
    if not warnings:
        return prediction

    return (
        f"{prediction} INVALID "
        f"({'; '.join(warnings)})"
    )


def write_trial_result(input_path, output_path):
    logger.info("Mode: trial_result")
    logger.info("Input: %s", input_path)
    logger.info("Output: %s", output_path)

    results = load_trial_result(input_path)

    with output_path.open("w") as f:
        for sequence_id, structures in results.items():
            f.write(f"{sequence_id}\n")

            f.write(
                f"input:       {structures['target']}\n"
            )

            prediction = format_prediction(
                structures["prediction"],
                structures["prediction_warnings"],
            )

            f.write(
                f"prediction:  {prediction}\n"
            )

            f.write("\n")

    logger.info(
        "Successfully wrote %d sequences to %s",
        len(results),
        output_path,
    )


def parse_config_ids(trajectory_path):
    """Return Config IDs in the order they appear."""
    logger.info(
        "Reading configuration trajectory: %s",
        trajectory_path,
    )

    text = trajectory_path.read_text()

    config_ids = [
        int(match)
        for match in re.findall(
            r"^Config ID:\s*(\d+)\s*$",
            text,
            re.MULTILINE,
        )
    ]

    logger.info(
        "Found %d Config IDs in trajectory: %s",
        len(config_ids),
        config_ids,
    )

    if not config_ids:
        raise ValueError(
            f"No Config IDs found in {trajectory_path}"
        )

    return config_ids


def write_trajectory(input_path, output_path):
    logger.info("Mode: trajectory")
    logger.info("Input: %s", input_path)
    logger.info("Output: %s", output_path)

    # best_config_trajectory.txt:
    #
    # neps/
    # ├── summary/
    # │   └── best_config_trajectory.txt
    # └── configs/
    #
    neps_dir = input_path.parent.parent
    configs_dir = neps_dir / "configs"

    logger.info(
        "NEPS directory: %s",
        neps_dir,
    )

    logger.info(
        "Configs directory: %s",
        configs_dir,
    )

    config_ids = parse_config_ids(input_path)

    config_results = []

    # Process configurations strictly in trajectory order.
    for config_number, config_id in enumerate(
        config_ids,
        start=1,
    ):
        trial_path = (
            configs_dir
            / f"config_{config_id}"
            / "artifacts"
            / "trial_result.json"
        )

        logger.info(
            "Processing trajectory entry %d/%d: "
            "Config ID %d",
            config_number,
            len(config_ids),
            config_id,
        )

        logger.debug(
            "Expected trial result: %s",
            trial_path,
        )

        # Missing configuration result is a hard failure.
        if not trial_path.is_file():
            logger.error(
                "Missing trial_result.json for Config ID %d: %s",
                config_id,
                trial_path,
            )

            raise FileNotFoundError(
                f"Missing trial_result.json for Config ID "
                f"{config_id}: {trial_path}"
            )

        try:
            results = load_trial_result(trial_path)
        except Exception:
            logger.exception(
                "Failed processing Config ID %d: %s",
                config_id,
                trial_path,
            )
            raise

        # Log invalid predictions with configuration context.
        for sequence_id, result in results.items():
            for warning in result["prediction_warnings"]:
                logger.warning(
                    "Config ID %d, sequence %s: %s",
                    config_id,
                    sequence_id,
                    warning,
                )

        config_results.append(
            (config_id, results)
        )

    # Use the first configuration for sequence ordering.
    first_config_id, first_results = config_results[0]

    sequence_ids = list(first_results.keys())

    logger.info(
        "Using Config ID %d as sequence reference",
        first_config_id,
    )

    logger.info(
        "Found %d sequences",
        len(sequence_ids),
    )

    expected_ids = set(sequence_ids)

    # Every configuration must contain exactly the same sequences.
    for config_id, results in config_results:
        actual_ids = set(results)

        if actual_ids != expected_ids:
            missing = expected_ids - actual_ids
            extra = actual_ids - expected_ids

            logger.error(
                "Sequence mismatch in Config ID %d. "
                "Missing: %s; Extra: %s",
                config_id,
                sorted(missing),
                sorted(extra),
            )

            raise ValueError(
                f"Sequence IDs in Config ID {config_id} do not "
                f"match those in Config ID {first_config_id}"
            )

    # Align all labels.
    labels = [
        "input:",
        *[
            f"CONFIG_ID-{config_id}:"
            for config_id in config_ids
        ],
    ]

    label_width = max(
        map(len, labels)
    )

    logger.debug(
        "Using label width %d",
        label_width,
    )

    with output_path.open("w") as f:
        for sequence_id in sequence_ids:
            target = first_results[sequence_id]["target"]

            f.write(
                f"{sequence_id}:\n"
            )

            f.write(
                f"{'input:':<{label_width}} "
                f"{target}\n"
            )

            for config_id, results in config_results:
                result = results[sequence_id]

                prediction = format_prediction(
                    result["prediction"],
                    result["prediction_warnings"],
                )

                label = f"CONFIG_ID-{config_id}:"

                f.write(
                    f"{label:<{label_width}} "
                    f"{prediction}\n"
                )

            f.write("\n")

    logger.info(
        "Successfully wrote trajectory to %s",
        output_path,
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract dot-bracket structures from NEPS "
            "trial results."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help=(
            "Input file: trial_result.json or "
            "best_config_trajectory.txt"
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output directory",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG logging",
    )

    args = parser.parse_args()

    configure_logging(
        args.verbose
    )

    input_path = args.input.resolve()
    output_dir = args.output.resolve()

    logger.info(
        "Starting dot-bracket extraction"
    )

    logger.info(
        "Input file: %s",
        input_path,
    )

    logger.info(
        "Output directory: %s",
        output_dir,
    )

    if not input_path.is_file():
        logger.error(
            "Input file does not exist: %s",
            input_path,
        )
        sys.exit(1)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:
        # Mode is determined exclusively from the input filename.
        if input_path.name == "trial_result.json":
            output_path = (
                output_dir / "dotbracket.txt"
            )

            write_trial_result(
                input_path,
                output_path,
            )

        elif input_path.name == "best_config_trajectory.txt":
            output_path = (
                output_dir / "dotbracket_trajectory.txt"
            )

            write_trajectory(
                input_path,
                output_path,
            )

        else:
            logger.error(
                "Unknown input filename: %s",
                input_path.name,
            )

            logger.error(
                "Expected trial_result.json or "
                "best_config_trajectory.txt"
            )

            sys.exit(1)

    except Exception:
        logger.exception(
            "Dot-bracket extraction failed"
        )
        sys.exit(1)

    logger.info(
        "Finished successfully"
    )


if __name__ == "__main__":
    main()

