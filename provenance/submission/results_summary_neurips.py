from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable
from pprint import pformat

import numpy as np
import pandas as pd


RNAPROTEIN_METRICS = [
    "Dockq",
    "ICS",
    "QS_global",
    "Complex_TM",
    "Complex_RMSD",
    "Complex_LDDT_allAtom",
    # "Complex_RMSD_allAtom",
    # "RNA_LDDT_allAtom",
]

RNAONLY_METRICS = [
    "RNA_TM",
    "RNA_RMSD",
    "RNA_LDDT_allAtom",
]

# Metrics that are emitted to LaTeX. Keep these near the top so paper-table
# subsets can be changed without touching the summarization logic.
LATEX_RNAPROTEIN_METRICS = [
    "Dockq",
    "ICS",
    "QS_global",
    "Complex_TM",
    "Complex_RMSD",
    "Complex_LDDT_allAtom",
    # "Complex LDDT",
    # "Complex_RMSD_allAtom",
    # "RNA_LDDT_allAtom",
]

LATEX_RNAONLY_METRICS = [
    "RNA_TM",
    "RNA_RMSD",
    "RNA_LDDT_allAtom",
]

# Reduced metric subsets that are also emitted to LaTeX for every experiment.
REDUCED_LATEX_RNAPROTEIN_METRICS = [
    "Dockq",
    "Complex_RMSD",
    "Complex_LDDT_allAtom",
]

REDUCED_LATEX_RNAONLY_METRICS = [
    "RNA_TM",
    "RNA_RMSD",
    "RNA_LDDT_allAtom",
]

RNAPROTEIN_ALGORITHMS = [
    "alphafold",
    # "noMSA",
    # "rnafoldN100",
    "spotrnaN100",
    "rnaformerN100",
    "rnaformerN5000",
    "rnaformerN10000",
    "rnaformerN20000",
    # "dssrN100",
    # "rnaformerN100SHSandMSA",
    # "dssrN100SHSandMSA",
]

RNAONLY_ALGORITHMS = [
    "alphafold",
    # "rnafoldN100",
    # "spotrnaN100",
    "rnaformerN100",
    # "rnaformerN20SHSandMSA",
    # "rnaformerN50SHSandMSA",
    # "rnaformerN100SHSandMSA",
    # "rnaformerN200SHSandMSA",
    # "rnaformerN500SHSandMSA",
    # "rnaformerN1000SHSandMSA",
    # "dssrN100SHSandMSA",
    # "rnaformerN100S0",
    # "rnaformerN100S42",
    # "rnaformerN100S137",
    # "noMSA",
    # "rnaformerN100ShuffleandMSA",
    # "rnaformerN100Shuffle",
    # "rnaformerN100RandomandMSA",
    # "rnaformerN100Random",
    "rnaformerN5000",
    "rnaformerN10000",
    "rnaformerN20000",
    # "rnafoldN100Denoise",
    # "dssrN100",
    # "spotrnaN100Denoise",
    # "rnaformerN100Denoise",
    # "rhofold",
    # "rnafoldN100Rhofold",
    # "spotrnaN100Rhofold",
    # "rnaformerN100Rhofold",    
    # "rnaformerN100SHSandMSARhofold",
    # "rnafoldN100DenoiseRhofold",
    # "dssrN100Rhofold",
    # "spotrnaN100DenoiseRhofold",
    # "rnaformerN100DenoiseRhofold",
    ]

DEFAULT_ALGORITHMS = sorted(set(RNAPROTEIN_ALGORITHMS + RNAONLY_ALGORITHMS))

LATEX_METHOD_NAMES = {
    "alphafold": "MSA",
    "spotrnaN100": "SHS (SPOT-RNA)",
    "rnafoldN100": "SHS (RNAfold)",
    "rnaformerN100": "SHS (RNAformer)",
    "rnaformerN100S0": "SHS (RNAformer) + MSA (AF seed=0)",
    "rnaformerN100S42": "SHS (RNAformer) + MSA (AF seed=42)",
    "rnaformerN100S137": "SHS (RNAformer) + MSA (AF seed=137)",
    "rnaformerN100SHSandMSA": "SHS (RNAformer) + MSA",
    "dssrN100SHSandMSA": "SHS (GT) + MSA",
    "dssrN100": "SHS (GT)",
    "noMSA": "MSA removed baseline",
    "rnaformerN100ShuffleandMSA": "SHS (Shuffled) + MSA",
    "rnaformerN100Shuffle": "SHS (Shuffled)",
    "rnaformerN100RandomandMSA": "SHS (Random) + MSA",
    "rnaformerN100Random": "SHS (Random)",
    "dssrN100Rhofold": "RF+ SHS (GT)",
    "rhofold": "RF+ MSA",
    "rnaformerN100DenoiseRhofold": "RF+ SHS (RNAformer) Denoise",
    "spotrnaN100DenoiseRhofold": "RF+ SHS (SPOT-RNA Denoise)",
    "rnafoldN100DenoiseRhofold": "RF+ SHS (RNAfold Denoise)",
    "rnafoldN100Rhofold": "RF+ SHS (RNAfold)",
    "rnaformerN100Rhofold": "RF+ SHS (RNAformer)",
    "spotrnaN100Rhofold": "RF+ SHS (SPOT-RNA)",
    "rnafoldN100Denoise": "SHS (RNAfold Denoise)",
    "rnaformerN100Denoise": "SHS (RNAformer Denoise)",
    "spotrnaN100Denoise": "SHS (SPOT-RNA Denoise)",
    "rnaformerN20SHSandMSA": "SHS (RNAformerN20) + MSA",
    "rnaformerN50SHSandMSA": "SHS (RNAformerN50) + MSA",
    "rnaformerN200SHSandMSA": "SHS (RNAformerN200) + MSA",
    "rnaformerN500SHSandMSA": "SHS (RNAformerN500) + MSA",
    "rnaformerN1000SHSandMSA": "SHS (RNAformerN1000) + MSA",
    "rnaformerN100SHSandMSARhofold": "RF+ SHS (RNAformer) + MSA"
}

LATEX_SUBSET_HEADERS = {
    "Full benchmark eval": "Full benchmark",
    "Full benchmark eval, orphan targets": "Full benchmark w/o Homologs",
    "Full benchmark eval, non-orphan targets": "Full benchmark w/ Homologs",
    "Full benchmark eval, before 2021": "Full benchmark before 2021",
    "Full benchmark eval, before 2021 orphan targets": (
        "Full benchmark before 2021 w/o Homologs"
    ),
    "Full benchmark eval, before 2021 non-orphan targets": (
        "Full benchmark before 2021 w/ Homologs"
    ),
    "Full benchmark eval, after 2021": "Full benchmark after 2021",
    "Full benchmark eval, after 2021 orphan targets": (
        "Full benchmark after 2021 w/o Homologs"
    ),
    "Full benchmark eval, after 2021 non-orphan targets": (
        "Full benchmark after 2021 w/ Homologs"
    ),
    "RNA-Protein eval": "RNA-protein complexes",
    "RNA-Protein eval, orphan targets": "RNA-protein complexes w/o Homologs",
    "RNA-Protein eval, non-orphan targets": "RNA-protein complexes w/ Homologs",
    "RNA-Protein eval, before 2021": "RNA-protein complexes before 2021",
    "RNA-Protein eval, before 2021 orphan targets": (
        "RNA-protein complexes before 2021 w/o Homologs"
    ),
    "RNA-Protein eval, before 2021 non-orphan targets": (
        "RNA-protein complexes before 2021 w/ Homologs"
    ),
    "RNA-Protein eval, after 2021": "RNA-protein complexes after 2021",
    "RNA-Protein eval, after 2021 orphan targets": (
        "RNA-protein complexes after 2021 w/o Homologs"
    ),
    "RNA-Protein eval, after 2021 non-orphan targets": (
        "RNA-protein complexes after 2021 w/ Homologs"
    ),
    "RNA-Protein eval, 1v1 targets": "RNA-protein 1v1 complexes",
    "RNA-Protein eval, multichain targets": "RNA-protein multichain complexes",
    "RNA-only eval": "RNA single chain",
    "RNA-only eval, orphan targets": "RNA single chain w/o Homologs",
    "RNA-only eval, non-orphan targets": "RNA single chain w/ Homologs",
    "RNA-only eval, before 2021": "RNA single chain before 2021",
    "RNA-only eval, before 2021 orphan targets": (
        "RNA single chain before 2021 w/o Homologs"
    ),
    "RNA-only eval, before 2021 non-orphan targets": (
        "RNA single chain before 2021 w/ Homologs"
    ),
    "RNA-only eval, after 2021": "RNA single chain after 2021",
    "RNA-only eval, after 2021 orphan targets": (
        "RNA single chain after 2021 w/o Homologs"
    ),
    "RNA-only eval, after 2021 non-orphan targets": (
        "RNA single chain after 2021 w/ Homologs"
    ),
}

LOWER_IS_BETTER = {
    "Complex_RMSD",
    "RNA_RMSD",
}

LATEX_METRIC_NAMES = {
    "Dockq": "DockQ",
    "ICS": "ICS",
    "QS_global": "QS Global",
    "Complex_TM": "TM Score",
    "Complex_RMSD": "cRMSD",
    "Complex_LDDT_allAtom": "cLDDT",
    "RNA_TM": "TM Score",
    "RNA_LDDT_allAtom": "LDDT",
    "RNA_RMSD": "RMSD",
}


SubsetFn = Callable[[pd.DataFrame], pd.DataFrame]


@dataclass(frozen=True)
class Experiment:
    title: str
    file_template: str
    metrics: list[str]
    baseline_algorithm: str = "alphafold"
    algorithms: list[str] | None = None
    subset: SubsetFn | None = None
    query: str | None = None


@dataclass(frozen=True)
class ExperimentSummary:
    table: pd.DataFrame
    final_ids: set[str]


def make_split_experiments(
    title_prefix: str,
    file_prefix: str,
    metrics: list[str],
) -> list[Experiment]:
    return [
        Experiment(
            title=title_prefix,
            file_template=f"{file_prefix}_{{algorithm}}.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, orphan targets",
            file_template=f"{file_prefix}_{{algorithm}}_orphan.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, non-orphan targets",
            file_template=f"{file_prefix}_{{algorithm}}_non_orphan.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, before 2021",
            file_template=f"{file_prefix}_b2021_{{algorithm}}.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, before 2021 orphan targets",
            file_template=f"{file_prefix}_b2021_{{algorithm}}_orphan.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, before 2021 non-orphan targets",
            file_template=f"{file_prefix}_b2021_{{algorithm}}_non_orphan.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, after 2021",
            file_template=f"{file_prefix}_a2021_{{algorithm}}.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, after 2021 orphan targets",
            file_template=f"{file_prefix}_a2021_{{algorithm}}_orphan.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
        Experiment(
            title=f"{title_prefix}, after 2021 non-orphan targets",
            file_template=f"{file_prefix}_a2021_{{algorithm}}_non_orphan.csv",
            metrics=metrics,
            # baseline_algorithm = "rhofold" if not "alphafold" in ,
        ),
    ]


# Edit this list to add new result tables.  Use either `query` for simple pandas
# query expressions or `subset` for arbitrary dataframe filtering.
EXPERIMENTS = [
    # *make_split_experiments("Full benchmark eval", "All", RNAPROTEIN_METRICS),
    *make_split_experiments("RNA-Protein eval", "RNA-Protein", RNAPROTEIN_METRICS),
    Experiment(
        title="RNA-Protein eval, 1v1 targets",
        file_template="RNA-Protein_{algorithm}_1v1.csv",
        metrics=RNAPROTEIN_METRICS,
        # baseline_algorithm = "rnaformerN100",
    ),
    Experiment(
        title="RNA-Protein eval, multichain targets",
        file_template="RNA-Protein_{algorithm}_multichain.csv",
        metrics=RNAPROTEIN_METRICS,
        # baseline_algorithm = "rnaformerN100",
    ),
    *make_split_experiments("RNA-only eval", "RNA_Monomers", RNAONLY_METRICS),
    # Example subset experiments for later:
    # Experiment(
    #     title="Full benchmark, orphan targets",
    #     file_template="All_{algorithm}.csv",
    #     metrics=RNAPROTEIN_METRICS,
    #     query="exp_db_id in @ORPHAN_IDS",
    # ),
    # Experiment(
    #     title="RNA-only, longer RNAs",
    #     file_template="RNA_Monomers_{algorithm}.csv",
    #     metrics=RNAONLY_METRICS,
    #     subset=lambda df: df[pd.to_numeric(df["RNALength"], errors="coerce") >= 100],
    # ),
]


def discover_algorithms(csv_dir: Path, experiments: list[Experiment]) -> list[str]:
    algorithms: set[str] = set()

    for experiment in experiments:
        if experiment.algorithms:
            algorithms.update(experiment.algorithms)
            continue

        if "{algorithm}" not in experiment.file_template:
            continue

        prefix, suffix = experiment.file_template.split("{algorithm}", maxsplit=1)
        for csv_path in csv_dir.glob(f"{prefix}*{suffix}"):
            algorithm = csv_path.name.removeprefix(prefix).removesuffix(suffix)
            if algorithm:
                algorithms.add(algorithm)

    return sorted(algorithms)


def is_rnaonly_experiment(experiment: Experiment) -> bool:
    return set(experiment.metrics).issubset(set(RNAONLY_METRICS))


def default_algorithms_for_experiment(experiment: Experiment) -> list[str]:
    if is_rnaonly_experiment(experiment):
        return RNAONLY_ALGORITHMS
    return RNAPROTEIN_ALGORITHMS


def algorithms_for_experiment(
    experiment: Experiment,
    cli_algorithms: list[str] | None = None,
) -> list[str]:
    return experiment.algorithms or cli_algorithms or default_algorithms_for_experiment(experiment)


def validate_experiment_algorithms(
    experiments: list[Experiment],
    cli_algorithms: list[str] | None = None,
) -> None:
    for experiment in experiments:
        algorithms = algorithms_for_experiment(experiment, cli_algorithms)
        if experiment.baseline_algorithm not in algorithms:
            raise ValueError(
                f"Baseline algorithm '{experiment.baseline_algorithm}' is not in the "
                f"algorithm list for experiment '{experiment.title}': {algorithms}"
            )


def coerce_metric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan)


def latex_metrics_for_experiment(experiment: Experiment) -> list[str]:
    if set(experiment.metrics).issubset(set(RNAONLY_METRICS)):
        configured_metrics = LATEX_RNAONLY_METRICS
    else:
        configured_metrics = LATEX_RNAPROTEIN_METRICS

    return [metric for metric in configured_metrics if metric in experiment.metrics]


def reduced_latex_metrics_for_experiment(experiment: Experiment) -> list[str]:
    if is_rnaonly_experiment(experiment):
        configured_metrics = REDUCED_LATEX_RNAONLY_METRICS
    else:
        configured_metrics = REDUCED_LATEX_RNAPROTEIN_METRICS

    return [metric for metric in configured_metrics if metric in experiment.metrics]


def should_write_reduced_latex_table(experiment: Experiment) -> bool:
    return (
        reduced_latex_metrics_for_experiment(experiment)
        != latex_metrics_for_experiment(experiment)
    )


def format_number(value: float) -> str:
    if pd.isna(value):
        return "nan"
    return f"{value:.3f}"


def format_mean_std(mean: float, std: float) -> str:
    return f"{format_number(mean)} +- {format_number(std)}"


def format_mean_std_latex(mean: float, std: float) -> str:
    return rf"{format_number(mean)} $\pm$ {format_number(std)}"


def apply_subset(df: pd.DataFrame, experiment: Experiment) -> pd.DataFrame:
    subset_df = df

    if experiment.query:
        subset_df = subset_df.query(experiment.query)

    if experiment.subset:
        subset_df = experiment.subset(subset_df)

    return subset_df.copy()


def summarize_metric(df: pd.DataFrame, metric: str) -> tuple[int, float, float, float]:
    values = coerce_metric(df[metric]).dropna()
    return (
        int(values.shape[0]),
        float(values.mean()) if not values.empty else np.nan,
        float(values.std(ddof=1)) if values.shape[0] > 1 else np.nan,
        float(values.median()) if not values.empty else np.nan,
    )


def summarize_experiment(
    csv_dir: Path,
    experiment: Experiment,
    algorithms: list[str],
) -> pd.DataFrame:
    return summarize_experiment_with_ids(csv_dir, experiment, algorithms).table


def summarize_experiment_with_ids(
    csv_dir: Path,
    experiment: Experiment,
    algorithms: list[str],
) -> ExperimentSummary:
    dataframes = load_experiment_dataframes(csv_dir, experiment, algorithms)
    complete_ids = complete_case_ids(dataframes, experiment.metrics)
    complete_dataframes = filter_dataframes_to_ids(dataframes, complete_ids)
    table = summarize_dataframes(complete_dataframes, experiment)
    return ExperimentSummary(table=table, final_ids=complete_ids)


def load_experiment_dataframes(
    csv_dir: Path,
    experiment: Experiment,
    algorithms: Iterable[str],
) -> dict[str, pd.DataFrame]:
    dataframes = {}

    for algorithm in algorithms:
        csv_path = csv_dir / experiment.file_template.format(algorithm=algorithm)
        if not csv_path.exists():
            continue

        dataframes[algorithm] = apply_subset(pd.read_csv(csv_path), experiment)

    return dataframes


def summarize_dataframes(
    dataframes: dict[str, pd.DataFrame],
    experiment: Experiment,
) -> pd.DataFrame:
    rows = []

    for algorithm, df in dataframes.items():
        missing_metrics = [metric for metric in experiment.metrics if metric not in df.columns]

        for metric in experiment.metrics:
            if metric in missing_metrics:
                rows.append(
                    {
                        "algorithm": algorithm,
                        "rows": len(df),
                        "metric": metric,
                        "n": 0,
                        "mean": np.nan,
                        "std": np.nan,
                        "median_value": np.nan,
                        "mean +- std": "missing",
                        "median": "missing",
                    }
                )
                continue

            n, mean, std, median = summarize_metric(df, metric)
            rows.append(
                {
                    "algorithm": algorithm,
                    "rows": len(df),
                    "metric": metric,
                    "n": n,
                    "mean": mean,
                    "std": std,
                    "median_value": median,
                    "mean +- std": format_mean_std(mean, std),
                    "median": format_number(median),
                }
            )

    return pd.DataFrame(rows)


def dataframe_ids(df: pd.DataFrame) -> set[str]:
    if "exp_db_id" not in df.columns:
        return set()
    return set(df["exp_db_id"].dropna().astype(str))


def filter_dataframes_to_ids(
    dataframes: dict[str, pd.DataFrame],
    ids: set[str],
) -> dict[str, pd.DataFrame]:
    filtered = {}
    for algorithm, df in dataframes.items():
        if "exp_db_id" not in df.columns:
            filtered[algorithm] = df.iloc[0:0].copy()
            continue

        filtered[algorithm] = df[df["exp_db_id"].astype(str).isin(ids)].copy()

    return filtered


def common_ids_with_baseline(
    dataframes: dict[str, pd.DataFrame],
    baseline_algorithm: str,
) -> set[str]:
    if baseline_algorithm not in dataframes:
        return set()

    id_sets = [dataframe_ids(dataframes[baseline_algorithm])]
    for algorithm, df in dataframes.items():
        if algorithm == baseline_algorithm:
            continue
        id_sets.append(dataframe_ids(df))

    if not id_sets or any(not ids for ids in id_sets):
        return set()

    return set.intersection(*id_sets)


def dataframes_have_matching_ids(
    dataframes: dict[str, pd.DataFrame],
    baseline_algorithm: str,
) -> bool:
    if baseline_algorithm not in dataframes:
        return True

    baseline_ids = dataframe_ids(dataframes[baseline_algorithm])
    return all(dataframe_ids(df) == baseline_ids for df in dataframes.values())


def sample_counts_by_algorithm(dataframes: dict[str, pd.DataFrame]) -> dict[str, int]:
    return {algorithm: int(len(df)) for algorithm, df in dataframes.items()}


def metric_missing_by_id(
    df: pd.DataFrame,
    metrics: list[str],
) -> dict[str, list[str]]:
    missing_by_id: dict[str, list[str]] = {}
    if "exp_db_id" not in df.columns:
        return missing_by_id

    ids = df["exp_db_id"].astype(str)
    for metric in metrics:
        if metric not in df.columns:
            missing_ids = ids.dropna().unique()
        else:
            missing_ids = ids[coerce_metric(df[metric]).isna()].dropna().unique()

        for exp_db_id in missing_ids:
            missing_by_id.setdefault(str(exp_db_id), []).append(metric)

    return missing_by_id


def complete_metric_ids(df: pd.DataFrame, metrics: list[str]) -> set[str]:
    if "exp_db_id" not in df.columns:
        return set()
    if any(metric not in df.columns for metric in metrics):
        return set()

    metric_values = df[metrics].apply(coerce_metric)
    complete_rows = metric_values.notna().all(axis=1)
    return set(df.loc[complete_rows, "exp_db_id"].dropna().astype(str))


def complete_case_ids(
    dataframes: dict[str, pd.DataFrame],
    metrics: list[str],
) -> set[str]:
    id_sets = [complete_metric_ids(df, metrics) for df in dataframes.values()]
    if not id_sets or any(not ids for ids in id_sets):
        return set()
    return set.intersection(*id_sets)


def format_id_metric_pairs(missing_by_id: dict[str, list[str]]) -> str:
    parts = []
    for exp_db_id in sorted(missing_by_id):
        metrics = ", ".join(sorted(set(missing_by_id[exp_db_id])))
        parts.append(f"{exp_db_id} ({metrics})")
    return "; ".join(parts)


def experiment_id_variable_name(experiment: Experiment) -> str | None:
    mapping = {
        "RNA-Protein eval": "protein_rna_ids",
        "RNA-Protein eval, orphan targets": "protein_rna_orphan_ids",
        "RNA-Protein eval, non-orphan targets": "protein_rna_non_orphan_ids",
        "RNA-Protein eval, before 2021": "protein_rna_b2021_ids",
        "RNA-Protein eval, before 2021 orphan targets": "protein_rna_b2021_orphan_ids",
        "RNA-Protein eval, before 2021 non-orphan targets": (
            "protein_rna_b2021_non_orphan_ids"
        ),
        "RNA-Protein eval, after 2021": "protein_rna_a2021_ids",
        "RNA-Protein eval, after 2021 orphan targets": "protein_rna_a2021_orphan_ids",
        "RNA-Protein eval, after 2021 non-orphan targets": (
            "protein_rna_a2021_non_orphan_ids"
        ),
        "RNA-Protein eval, 1v1 targets": "protein_rna_1v1_ids",
        "RNA-Protein eval, multichain targets": "protein_rna_multichain_ids",
        "RNA-only eval": "rna_monomers_ids",
        "RNA-only eval, orphan targets": "rna_monomers_orphan_ids",
        "RNA-only eval, non-orphan targets": "rna_monomers_non_orphan_ids",
        "RNA-only eval, before 2021": "rna_monomers_b2021_ids",
        "RNA-only eval, before 2021 orphan targets": "rna_monomers_b2021_orphan_ids",
        "RNA-only eval, before 2021 non-orphan targets": (
            "rna_monomers_b2021_non_orphan_ids"
        ),
        "RNA-only eval, after 2021": "rna_monomers_a2021_ids",
        "RNA-only eval, after 2021 orphan targets": "rna_monomers_a2021_orphan_ids",
        "RNA-only eval, after 2021 non-orphan targets": (
            "rna_monomers_a2021_non_orphan_ids"
        ),
    }
    return mapping.get(experiment.title)


def sorted_ids(ids: Iterable[str]) -> list[str]:
    return sorted(set(ids))


def add_union_id_lists(id_lists: dict[str, list[str]]) -> None:
    subset_suffixes = [
        "",
        "_orphan",
        "_non_orphan",
        "_b2021",
        "_b2021_orphan",
        "_b2021_non_orphan",
        "_a2021",
        "_a2021_orphan",
        "_a2021_non_orphan",
    ]

    for suffix in subset_suffixes:
        rna_key = f"rna_monomers{suffix}_ids"
        protein_key = f"protein_rna{suffix}_ids"
        all_key = f"all_data{suffix}_ids"
        rna_ids = id_lists.get(rna_key, [])
        protein_ids = id_lists.get(protein_key, [])
        if rna_ids or protein_ids:
            id_lists[all_key] = sorted_ids([*rna_ids, *protein_ids])


def make_neurips_id_lists(
    csv_dir: Path,
    experiments: list[Experiment],
    cli_algorithms: list[str] | None = None,
) -> dict[str, list[str]]:
    id_lists: dict[str, list[str]] = {}

    for experiment in experiments:
        variable_name = experiment_id_variable_name(experiment)
        if variable_name is None:
            continue

        algorithms = algorithms_for_experiment(experiment, cli_algorithms)
        summary = summarize_experiment_with_ids(csv_dir, experiment, algorithms)
        id_lists[variable_name] = sorted_ids(summary.final_ids)

    add_union_id_lists(id_lists)
    return dict(sorted(id_lists.items()))


def write_neurips_ids(
    output_path: Path,
    csv_dir: Path,
    experiments: list[Experiment],
    cli_algorithms: list[str] | None = None,
) -> None:
    id_lists = make_neurips_id_lists(csv_dir, experiments, cli_algorithms)
    lines = [
        '"""Final complete-case NeurIPS evaluation IDs.',
        "",
        "Generated by results_summary_neurips.py. These lists are intended to",
        "replace the manually maintained ID lists in preprocess_neurips_evaluation_data.py.",
        '"""',
        "",
    ]

    for name, ids in id_lists.items():
        lines.append(f"{name} = {pformat(ids, width=88)}")
        lines.append("")

    lines.append("__all__ = [")
    for name in id_lists:
        lines.append(f"    {name!r},")
    lines.append("]")
    lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def strict_complete_case_note(
    dataframes: dict[str, pd.DataFrame],
    experiment: Experiment,
    metrics: list[str],
    complete_ids: set[str],
) -> str:
    counts = sample_counts_by_algorithm(dataframes)
    count_parts = [
        f"{latex_method_name(algorithm)} n={count}"
        for algorithm, count in counts.items()
    ]
    note_parts = [
        rf"\footnotesize\emph{{Note:}} Baseline is "
        rf"{latex_escape(latex_method_name(experiment.baseline_algorithm))}; "
        rf"reported results use common complete-case \texttt{{exp\_db\_id}} targets "
        rf"across all methods and metrics (n={len(complete_ids)}). "
        "Original evaluated sample counts were "
        f"{latex_escape(', '.join(count_parts))}."
    ]

    if experiment.baseline_algorithm in dataframes:
        baseline_ids = dataframe_ids(dataframes[experiment.baseline_algorithm])
    else:
        baseline_ids = set()

    method_notes = []
    for algorithm, df in dataframes.items():
        ids = dataframe_ids(df)
        missing_ids = sorted(baseline_ids - ids)
        excluded_ids = sorted(ids - complete_ids)
        missing_metric_ids = metric_missing_by_id(df, metrics)

        details = []
        if missing_ids:
            details.append("missing IDs: " + ", ".join(missing_ids))
        if missing_metric_ids:
            details.append("missing metrics: " + format_id_metric_pairs(missing_metric_ids))

        # Include complete-case exclusions even when they are caused by another
        # method, so each method's row count can be traced back to the final set.
        other_exclusions = [
            exp_db_id
            for exp_db_id in excluded_ids
            if exp_db_id not in missing_metric_ids and exp_db_id not in missing_ids
        ]
        if other_exclusions:
            details.append("dropped for common set: " + ", ".join(other_exclusions))

        if details:
            method_notes.append(
                f"{latex_method_name(algorithm)}: " + "; ".join(details)
            )

    if method_notes:
        note_parts.append(latex_escape(" ".join(method_notes)))

    return " ".join(note_parts)


def print_experiment_table(
    title: str,
    table: pd.DataFrame,
) -> None:
    print()
    print("=" * len(title))
    print(title)
    print("=" * len(title))

    if table.empty:
        print("No matching CSV files found.")
        return

    for algorithm, algorithm_table in table.groupby("algorithm", sort=False):
        print()
        print(f"{algorithm} (rows: {int(algorithm_table['rows'].iloc[0])})")
        print(
            algorithm_table[["metric", "n", "mean +- std", "median"]].to_string(
                index=False,
                justify="left",
            )
        )


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def latex_bold(text: str) -> str:
    return rf"\textbf{{{text}}}"


def latex_result_cell(mean: float, std: float) -> str:
    return format_mean_std_latex(mean, std)


def best_mask(values: pd.Series, metric: str) -> pd.Series:
    numeric_values = pd.to_numeric(values, errors="coerce")
    if numeric_values.dropna().empty:
        return pd.Series(False, index=values.index)

    best = numeric_values.min() if metric in LOWER_IS_BETTER else numeric_values.max()
    return numeric_values.eq(best)


def latex_metrics_for_table(table: pd.DataFrame) -> list[str]:
    available_metrics = set(table["metric"]) if not table.empty else set()
    if available_metrics.issubset(set(RNAONLY_METRICS)):
        configured_metrics = LATEX_RNAONLY_METRICS
    else:
        configured_metrics = LATEX_RNAPROTEIN_METRICS

    return [metric for metric in configured_metrics if metric in available_metrics]


def latex_method_name(algorithm: str) -> str:
    return LATEX_METHOD_NAMES.get(algorithm, algorithm)


def title_case_label(label: str) -> str:
    replacements = {
        "RNA-protein": "RNA-Protein",
        "RNA single chain": "RNA Single Chain",
        "Full benchmark": "Full Benchmark",
        "complexes": "Complexes",
        "complex": "Complex",
        "before": "Before",
        "after": "After",
        "targets": "Targets",
        "without homologs": "Without Homologs",
        "with homologs": "With Homologs",
        "orphan": "Orphan",
        "non-orphan": "Non-Orphan",
        "multichain": "Multichain",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)
    return label


def latex_metric_name(metric: str) -> str:
    return LATEX_METRIC_NAMES.get(metric, metric.replace("_", " "))


def latex_metric_header(metric: str) -> str:
    direction = r"$\downarrow$" if metric in LOWER_IS_BETTER else r"$\uparrow$"
    return f"{latex_metric_name(metric)} ({direction})"


def latex_subset_header(title: str) -> str:
    for suffix in [
        " (all available targets)",
        " (matched targets)",
        " (reduced metrics)",
    ]:
        if title.endswith(suffix):
            title = title[: -len(suffix)]
    return LATEX_SUBSET_HEADERS.get(title, title)


def latex_display_header(title: str) -> str:
    header = latex_subset_header(title)
    header = (
        header.replace("w/o Homologs", "without homologs")
        .replace("w/ Homologs", "with homologs")
        .replace("before 2021", "Pre-Cutoff")
        .replace("after 2021", "Post-Cutoff")
    )
    return title_case_label(header)


def is_reduced_metrics_table(title: str) -> bool:
    return title.endswith(" (reduced metrics)")


def latex_metric_note(metrics: list[str]) -> str:
    note_parts = []
    if "Complex_RMSD" in metrics:
        note_parts.append("cRMSD denotes complex RMSD computed on backbone P atoms")
    if "Complex_LDDT_allAtom" in metrics:
        note_parts.append("cLDDT denotes all-atom complex LDDT")
    if "RNA_RMSD" in metrics:
        note_parts.append("RMSD is computed on RNA backbone P atoms")
    if "RNA_TM" in metrics:
        note_parts.append(r"TM Score is computed on RNA backbone C3$'$ atoms")
    if "RNA_LDDT_allAtom" in metrics:
        note_parts.append("LDDT denotes all-atom RNA LDDT")

    if not note_parts:
        return ""

    return " " + "; ".join(note_parts) + "."


def latex_table_caption(
    title: str,
    sample_count: int | None,
    metrics: list[str],
) -> str:
    subset_header = latex_escape(latex_display_header(title))
    count_text = (
        f" (dataset size: n={sample_count})" if sample_count is not None else ""
    )
    return (
        f"{subset_header} Results{count_text}. "
        "Values are mean $\\pm$ standard deviation; bold values indicate the best "
        "performance for each metric."
        f"{latex_metric_note(metrics)}"
    )


def latex_table_label(title: str) -> str:
    label = latex_subset_header(title).lower()
    replacements = {
        "w/o": "without",
        "w/": "with",
        "before 2021": "pre cutoff",
        "after 2021": "post cutoff",
        "rna-protein": "rna protein",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)

    label = re.sub(r"[^a-z0-9]+", "_", label).strip("_")
    if is_reduced_metrics_table(title):
        label = f"{label}_reduced"
    return f"tab:{label}"


def table_sample_count(table: pd.DataFrame) -> int | None:
    rows = pd.to_numeric(table["rows"], errors="coerce").dropna()
    if rows.empty:
        return None

    unique_counts = rows.astype(int).unique()
    if len(unique_counts) == 1:
        return int(unique_counts[0])

    return int(rows.min())


def latex_input_header(table: pd.DataFrame) -> str:
    algorithms = table["algorithm"].dropna().astype(str)
    if algorithms.str.lower().str.contains("rhofold").any():
        return "RhoFold+ Input"
    return "AlphaFold~3 Input"


def counts_note(
    dataframes: dict[str, pd.DataFrame],
    baseline_algorithm: str,
    common_count: int | None = None,
) -> str:
    counts = sample_counts_by_algorithm(dataframes)
    count_parts = [
        f"{latex_method_name(algorithm)} n={count}"
        for algorithm, count in counts.items()
    ]
    note = (
        rf"\footnotesize\emph{{Note:}} Baseline is "
        rf"{latex_escape(latex_method_name(baseline_algorithm))}; "
        "original evaluated sample counts are "
        f"{latex_escape(', '.join(count_parts))}."
    )
    if common_count is not None:
        note += rf" Matched table uses common \texttt{{exp\_db\_id}} targets (n={common_count})."
    return note


def table_to_latex(
    title: str,
    table: pd.DataFrame,
    metrics: list[str] | None = None,
) -> str:
    if table.empty:
        return "\n".join(
            [
                rf"\begin{{table}}[htbp]",
                r"\centering",
                rf"\caption{{{latex_escape(latex_subset_header(title))} results. No matching CSV files found.}}",
                rf"\label{{{latex_table_label(title)}}}",
                rf"\end{{table}}",
                "No matching CSV files found.",
                "",
            ]
        )

    table = table.copy()
    if metrics is None:
        metrics = latex_metrics_for_table(table)
    else:
        available_metrics = set(table["metric"])
        metrics = [metric for metric in metrics if metric in available_metrics]
    table = table[table["metric"].isin(metrics)]

    sample_count = table_sample_count(table)
    input_header = latex_input_header(table)
    subset_header = latex_display_header(title)
    subset_header_with_count = subset_header
    if sample_count is not None:
        subset_header_with_count = f"{subset_header} (n={sample_count})"

    for metric, metric_table in table.groupby("metric", sort=False):
        table.loc[metric_table.index, "is_best"] = best_mask(metric_table["mean"], metric)

    n_columns = len(metrics) + 1
    metric_columns = "c" * len(metrics)
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{latex_table_caption(title, sample_count, metrics)}}}",
        rf"\label{{{latex_table_label(title)}}}",
        r"\resizebox{\textwidth}{!}{%",
        rf"\begin{{tabular}}{{l{metric_columns}}}",
        r"\toprule",
        rf"\multirow{{2}}{{*}}{{{input_header}}} & \multicolumn{{{len(metrics)}}}{{c}}{{{latex_escape(subset_header_with_count)}}} \\",
        rf"\cmidrule(lr){{2-{n_columns}}}",
        " & " + " & ".join(latex_metric_header(metric) for metric in metrics) + r" \\",
        r"\midrule",
    ]

    for algorithm in table["algorithm"].drop_duplicates():
        algorithm_table = table[table["algorithm"] == algorithm].set_index("metric")
        cells = [latex_escape(latex_method_name(algorithm))]

        for metric in metrics:
            if metric not in algorithm_table.index:
                cells.append("--")
                continue

            row = algorithm_table.loc[metric]
            if pd.isna(row["mean"]):
                cells.append("missing")
                continue

            cell = latex_result_cell(row["mean"], row["std"])
            if bool(row.get("is_best", False)):
                cell = latex_bold(cell)
            cells.append(cell)

        lines.append(" & ".join(cells) + r" \\")

    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"}",
            r"\end{table}",
        ]
    )
    lines.append("")
    return "\n".join(lines)


def experiment_to_latex(
    csv_dir: Path,
    experiment: Experiment,
    algorithms: list[str],
    metrics: list[str] | None = None,
    title_suffix: str = "",
) -> str:
    dataframes = load_experiment_dataframes(csv_dir, experiment, algorithms)
    title = f"{experiment.title}{title_suffix}"
    if not dataframes:
        return table_to_latex(title, pd.DataFrame(), metrics=metrics)

    if metrics is None:
        metrics = latex_metrics_for_experiment(experiment)
    complete_ids = complete_case_ids(dataframes, metrics)
    complete_dataframes = filter_dataframes_to_ids(dataframes, complete_ids)
    complete_table = summarize_dataframes(complete_dataframes, experiment)

    return table_to_latex(title, complete_table, metrics=metrics)


def experiment_to_reduced_latex(
    csv_dir: Path,
    experiment: Experiment,
    algorithms: list[str],
) -> str:
    return experiment_to_latex(
        csv_dir,
        experiment,
        algorithms,
        metrics=reduced_latex_metrics_for_experiment(experiment),
        title_suffix=" (reduced metrics)",
    )


def latex_output_path(output_path: Path, suffix: str) -> Path:
    return output_path.with_name(f"{output_path.stem}_{suffix}{output_path.suffix}")


def latex_document_lines(title: str) -> list[str]:
    return [
        r"\documentclass{article}",
        r"\usepackage{booktabs}",
        r"\usepackage{graphicx}",
        r"\usepackage{multirow}",
        r"\usepackage[margin=0.7in]{geometry}",
        r"\begin{document}",
        rf"\section*{{{latex_escape(title)}}}",
        "",
    ]


def write_latex_summary(
    output_path: Path,
    csv_dir: Path,
    experiments: list[Experiment],
    cli_algorithms: list[str] | None = None,
) -> tuple[Path, Path]:
    main_output_path = latex_output_path(output_path, "main")
    appendix_output_path = latex_output_path(output_path, "appendix")
    main_lines = latex_document_lines("NeurIPS Results Summary: Main Tables")
    appendix_lines = latex_document_lines("NeurIPS Results Summary: Appendix Tables")

    for experiment in experiments:
        algorithms = algorithms_for_experiment(experiment, cli_algorithms)
        appendix_lines.append(experiment_to_latex(csv_dir, experiment, algorithms))
        if should_write_reduced_latex_table(experiment):
            main_lines.append(experiment_to_reduced_latex(csv_dir, experiment, algorithms))
        else:
            main_lines.append(experiment_to_latex(csv_dir, experiment, algorithms))

    main_lines.append(r"\end{document}")
    appendix_lines.append(r"\end{document}")
    main_output_path.write_text("\n".join(main_lines), encoding="utf-8")
    appendix_output_path.write_text("\n".join(appendix_lines), encoding="utf-8")
    return main_output_path, appendix_output_path


def print_interface_metric_sanity_check(csv_dir: Path, algorithms: list[str]) -> None:
    metrics = ["Dockq", "ICS", "QS_global"]
    rows = []

    for algorithm in algorithms:
        all_path = csv_dir / f"All_{algorithm}.csv"
        rnaprotein_path = csv_dir / f"RNA-Protein_{algorithm}.csv"
        if not all_path.exists() or not rnaprotein_path.exists():
            continue

        all_df = pd.read_csv(all_path)
        rnaprotein_df = pd.read_csv(rnaprotein_path)

        for metric in metrics:
            if metric not in all_df.columns or metric not in rnaprotein_df.columns:
                continue

            all_values = coerce_metric(all_df[metric]).dropna()
            rnaprotein_values = coerce_metric(rnaprotein_df[metric]).dropna()
            medians_match = np.isclose(all_values.median(), rnaprotein_values.median())
            counts_match = len(all_values) == len(rnaprotein_values)
            rows.append(
                {
                    "algorithm": algorithm,
                    "metric": metric,
                    "n": len(all_values),
                    "zeros": int((all_values == 0).sum()),
                    "nonzero": int((all_values != 0).sum()),
                    "median": format_number(float(all_values.median())),
                    "All/RNA-Protein match": "yes"
                    if medians_match and counts_match
                    else "check",
                }
            )

    if not rows:
        return

    print()
    print("=============================")
    print("Interface metric sanity check")
    print("=============================")
    print(
        pd.DataFrame(rows).to_string(
            index=False,
            justify="left",
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Print NeurIPS evaluation summary tables from preprocessed CSVs."
    )
    parser.add_argument(
        "--csv-dir",
        type=Path,
        default=Path("results/csvs_neurips"),
        help="Directory containing CSVs produced by preprocess_neurips_evaluation_data.py.",
    )
    parser.add_argument(
        "--algorithms",
        nargs="+",
        default=None,
        help=(
            "Override algorithms for all experiments. By default RNA-protein "
            "experiments use RNAPROTEIN_ALGORITHMS and RNA-only experiments use "
            "RNAONLY_ALGORITHMS."
        ),
    )
    parser.add_argument(
        "--tex-out",
        type=Path,
        default=Path("results_summary_neurips.tex"),
        help=(
            "Base path for LaTeX summaries. The script writes sibling files with "
            "_main and _appendix appended to the stem."
        ),
    )
    parser.add_argument(
        "--ids-out",
        type=Path,
        default=Path("neurips_ids.py"),
        help="Path for the generated Python file containing final complete-case IDs.",
    )
    parser.add_argument(
        "--write-ids",
        action="store_true",
        help="Write final complete-case IDs to --ids-out.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    csv_dir = args.csv_dir

    if not csv_dir.exists():
        raise FileNotFoundError(f"CSV directory does not exist: {csv_dir}")

    validate_experiment_algorithms(EXPERIMENTS, args.algorithms)

    if not any(
        (csv_dir / experiment.file_template.format(algorithm=algorithm)).exists()
        for experiment in EXPERIMENTS
        for algorithm in algorithms_for_experiment(experiment, args.algorithms)
    ):
        raise FileNotFoundError(f"No matching experiment CSV files found in {csv_dir}")

    print_interface_metric_sanity_check(
        csv_dir,
        algorithms_for_experiment(EXPERIMENTS[0], args.algorithms),
    )

    for experiment in EXPERIMENTS:
        algorithms = algorithms_for_experiment(experiment, args.algorithms)
        summary = summarize_experiment_with_ids(csv_dir, experiment, algorithms)
        print_experiment_table(experiment.title, summary.table)

    main_tex_path, appendix_tex_path = write_latex_summary(
        args.tex_out,
        csv_dir,
        EXPERIMENTS,
        args.algorithms,
    )
    print()
    print(f"Wrote main LaTeX summary to {main_tex_path}")
    print(f"Wrote appendix LaTeX summary to {appendix_tex_path}")
    if args.write_ids:
        write_neurips_ids(args.ids_out, csv_dir, EXPERIMENTS, args.algorithms)
        print(f"Wrote final complete-case IDs to {args.ids_out}")


if __name__ == "__main__":
    main()
