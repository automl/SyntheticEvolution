#!/usr/bin/env python3

import json
import logging
import re
from collections import Counter
from pathlib import Path


logger = logging.getLogger("analysis.dotbracket")


def run(context, output_dir):
    """Create dot-bracket output for every exported trajectory configuration."""
    trajectory_path = context.results_dir / "summary" / "best_config_trajectory.txt"
    output_dir.mkdir(parents=True, exist_ok=True)
    warning_count = write_trajectory(
        trajectory_path,
        output_dir / "dotbracket_trajectory.txt",
        config_results={
            config_id: context.trial_results[config_id]
            for config_id in context.config_ids
        },
    )
    if warning_count:
        logger.warning(
            "Dot-bracket conversion skipped %d invalid or conflicting pairs.",
            warning_count,
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


def load_trial_result_from_data(data, source):
    """Convert an already loaded trial-result object to dot-bracket records."""
    if "scores" not in data:
        raise ValueError(
            f"Missing 'scores' field in {source}"
        )

    results = {}

    logger.info("Loaded %d sequences for Config ID %s", len(data["scores"]), source)

    for index, score in enumerate(data["scores"], start=1):
        sequence_id = score.get("id")

        if sequence_id is None:
            raise ValueError(
                f"{source}: score #{index} has no 'id'"
            )

        logger.debug(
            "Processing sequence %s from %s",
            sequence_id,
            source,
        )

        try:
            target_pairs = score["target_pairs"]
            predicted_pairs = score["predicted_pairs"]
        except KeyError as e:
            raise ValueError(
                f"{source}: sequence {sequence_id} is missing "
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


def write_trajectory(input_path, output_path, config_results=None):
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
    if config_results is None:
        raise ValueError("config_results are required for adapter execution")

    config_ids = parse_config_ids(input_path)

    trajectory_results = []
    warning_count = Counter()

    # Process configurations strictly in trajectory order.
    for config_number, config_id in enumerate(
        config_ids,
        start=1,
    ):
        logger.info(
            "Processing trajectory entry %d/%d: "
            "Config ID %d",
            config_number,
            len(config_ids),
            config_id,
        )

        if config_id not in config_results:
            raise KeyError(f"Missing trial result for Config ID {config_id}")

        try:
            results = load_trial_result_from_data(config_results[config_id], config_id)
        except Exception:
            logger.exception("Failed processing Config ID %d", config_id)
            raise

        # Log invalid predictions with configuration context.
        for sequence_id, result in results.items():
            for warning in result["prediction_warnings"]:
                warning_count[warning] += 1
            for warning in result["target_warnings"]:
                warning_count[warning] += 1

        trajectory_results.append(
            (config_id, results)
        )

    # Use the first configuration for sequence ordering.
    first_config_id, first_results = trajectory_results[0]

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
    for config_id, results in trajectory_results:
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

            for config_id, results in trajectory_results:
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
    return sum(warning_count.values())
