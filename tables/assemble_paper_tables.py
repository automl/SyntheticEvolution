#!/usr/bin/env python3
"""Assemble copy-paste-ready LaTeX for the manuscript from the generated tables (results/tables/*.tex).

  python tables/assemble_paper_tables.py --manuscript <paper.tex> --dest <dir>

<dest>/main/<name>.tex                      main-text tables, one file each
<dest>/appendix_additional_results.tex      the full "Additional Results" appendix section in manuscript order: section and
                                            subsection headers, all tables (generated), the figure environments of that
                                            section copied verbatim from the manuscript, and the current prose of each
                                            subsection as LaTeX comments (old numbers are never emitted as text)
<dest>/appendix_sections/<name>.tex         tables for appendix sections outside "Additional Results"
<dest>/README.md                            index: file -> label -> placement -> generator
Every table is copied byte-for-byte from results/tables (verified by tables/verify_tables.py).
"""
import argparse, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"

MAIN = [("rna_single_chain_controls", "Sec. 6.1 (Topology, not depth)"),
        ("order_rna_single_chain", "Sec. 6.2 (Better 2D inputs lead to better 3D predictions)"),
        ("gates_single", "Sec. 6.2 (When SHS help)"),
        ("denoise_main", "Sec. 6.3 (Controlled interventions, denoising)"),
        ("rna_protein_complexes_post_cutoff_without_homologs_reduced", "Sec. 7.2 (RNA-protein complexes)")]

# (header line, [items]); item = table name | ("fig", label) | ("prose", label-of-section-in-manuscript)
APPENDIX = [
    (r"\section{Additional Results}\label{app:additional_results}", [("prose", "app:additional_results")]),
    (r"\subsection{Performance on All Single-Chain RNAs}\label{performance_rna_only}",
     [("prose", "performance_rna_only"), "rna_single_chain", "rna_single_chain_without_homologs", "rna_single_chain_with_homologs",
      "rna_single_chain_pre_cutoff", "rna_single_chain_pre_cutoff_without_homologs", "rna_single_chain_pre_cutoff_with_homologs",
      "rna_single_chain_post_cutoff", "rna_single_chain_post_cutoff_without_homologs", "rna_single_chain_post_cutoff_with_homologs"]),
    (r"\subsection{Performance on RNA-Protein Complexes}\label{app:performance_rp}",
     [("prose", "app:performance_rp"), "rna_protein_complexes", "rna_protein_complexes_without_homologs", "rna_protein_complexes_with_homologs",
      "rna_protein_complexes_pre_cutoff", "rna_protein_complexes_pre_cutoff_without_homologs", "rna_protein_complexes_pre_cutoff_with_homologs",
      "rna_protein_complexes_post_cutoff", "rna_protein_complexes_post_cutoff_without_homologs", "rna_protein_complexes_post_cutoff_with_homologs",
      "rna_protein_1v1_complexes", "rna_protein_multichain_complexes"]),
    (r"\subsection{RhoFold+ Single-Chain Results}\label{app:rhofold}",
     [("prose", "app:rhofold"), "rf_rna_single_chain", "rf_rna_single_chain_without_homologs", "rf_rna_single_chain_with_homologs",
      "rf_rna_single_chain_pre_cutoff", "rf_rna_single_chain_post_cutoff", "rf_rna_single_chain_post_cutoff_without_homologs",
      "rf_rna_single_chain_other_splits", "rf_rna_single_chain_native_msa"]),
    (r"\subsection{SHS and MSA}\label{app:shs_and_msa}", [("prose", "app:shs_and_msa"), "shs_and_msa_depth_sweep"]),
    (r"\subsection{SHS Depth Scaling Study}\label{app:scaling_depth}",
     [("prose", "app:scaling_depth"), "rna_single_chain_scaling", "rna_protein_complexes_scaling"]),
    (r"\subsection{Seed Robustness}\label{app:seed_robustness}", ["seed_robustness_table"]),
    (r"\subsection{Interventions}", []),
    (r"\subsubsection{Denoising}\label{app:rhofold_denoise}",
     [("prose", "app:rhofold_denoise"), "denoise_main_rf", ("fig", "fig:denoise_viual")]),
    (r"\subsubsection{Controlled Construction}\label{app:construction}",
     [("prose", "app:construction"), "construction_stage_mcc", ("fig", "fig:construction_sweep"), ("fig", "fig:construction_individual")]),
    (r"\subsubsection{Rescue}", [("prose", "Rescue"), ("fig", "fig:sup_rescue")]),
]

SECTIONS = [("gates_complex1v1", "App. When to trust SHS (app:trust), complexes", "tab:gates_complex"),
            ("orphan_tests", "App. paired orphan tests (app:orphan_tests)", "tab:orphan_tests"),
            ("near_floor", "App. near-floor orphans (app:near_floor)", "tab:near_floor"),
            ("long_rna", "App. RNAs longer than 200 nt (app:long_rna)", "tab:long_rna"),
            ("confidence_ptm", "App. AF3 confidence (app:confidence)", "tab:confidence"),
            ("novelty", "App. post-cutoff novelty (app:novelty)", "tab:novelty")]


def env_blocks(text, env):
    return [m for m in re.finditer(r"\\begin\{" + env + r"\*?\}.*?\\end\{" + env + r"\*?\}", text, re.S)]


def figure(text, label):
    for m in env_blocks(text, "figure"):
        if f"\\label{{{label}}}" in m.group(0): return m.group(0)
    raise SystemExit(f"figure {label} not found in manuscript")


def prose(text, key):
    """current prose of the (sub)section whose header contains `key`, as LaTeX comments."""
    hdr = re.search(r"\\(?:sub)*section\*?\{[^\n]*" + re.escape(key) + r"[^\n]*", text)
    if not hdr: return ""
    rest = text[hdr.end():]
    nxt = re.search(r"\\(?:sub)*section\*?\{", rest); seg = rest[: nxt.start()] if nxt else rest
    for env in ("table", "figure"):
        seg = re.sub(r"\\begin\{" + env + r"\*?\}.*?\\end\{" + env + r"\*?\}", "", seg, flags=re.S)
    lines = [l for l in seg.splitlines() if l.strip() and not l.lstrip().startswith("%")]
    if not lines: return ""
    return ("% --- current manuscript prose of this section (commented out: re-check every number before using) ---\n"
            + "\n".join("% " + l for l in lines) + "\n% ---\n")


def tab(name):
    return (TAB / f"{name}.tex").read_text()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manuscript", required=True); ap.add_argument("--dest", required=True)
    a = ap.parse_args(); text = Path(a.manuscript).read_text(); dest = Path(a.dest)
    (dest / "main").mkdir(parents=True, exist_ok=True); (dest / "appendix_sections").mkdir(parents=True, exist_ok=True)
    idx = ["| File | Label | Placement | Generator |", "|---|---|---|---|"]
    lab = lambda t: re.search(r"\\label\{([^}]*)\}", t).group(1)
    gen = lambda n: "make_extra_tables.py" if n in EXTRA else "make_crc_tables.py"
    for n, where in MAIN:
        t = tab(n); (dest / "main" / f"{n}.tex").write_text(t); idx.append(f"| main/{n}.tex | `{lab(t)}` | {where} | {gen(n)} |")
    out = ["% Additional Results appendix, generated by tables/assemble_paper_tables.py (do not edit numbers by hand).", ""]
    for header, items in APPENDIX:
        out += [header, ""]
        for it in items:
            if isinstance(it, tuple) and it[0] == "fig": out += [figure(text, it[1]), ""]
            elif isinstance(it, tuple) and it[0] == "prose":
                p = prose(text, it[1])
                if p: out += [p]
            else:
                t = tab(it); out += [t]
                idx.append(f"| appendix_additional_results.tex | `{lab(t)}` | {header.split('{')[1].split('}')[0]} | {gen(it)} |")
    (dest / "appendix_additional_results.tex").write_text("\n".join(out))
    for n, where, _ in SECTIONS:
        t = tab(n); (dest / "appendix_sections" / f"{n}.tex").write_text(t)
        idx.append(f"| appendix_sections/{n}.tex | `{lab(t)}` | {where} | {gen(n)} |")
    readme = ["# CRC tables, copy-paste ready", "",
              "Generated from the SyntheticEvolution repo (branch reproduce-neurips): `tables/make_crc_tables.py`, "
              "`tables/make_extra_tables.py`, assembled by `tables/assemble_paper_tables.py`; checked by `tables/verify_tables.py` "
              "(see VERIFY_REPORT.md). Do not edit numbers by hand: regenerate instead.", "",
              "* `main/`: one file per main-text table (paste the whole table environment).",
              "* `appendix_additional_results.tex`: the complete Additional Results section (replace the section in the manuscript "
              "from `\\section{Additional Results}` up to the paper checklist). Figure environments are copied from the current "
              "manuscript; the current prose is included as comments.",
              "* `appendix_sections/`: tables for appendix sections outside Additional Results (sections still to be written).",
              "* Removed relative to the submission: the per-seed SHS+MSA table (`tab:shs_and_msa`), replaced by `tab:seed_robustness`.",
              "", *idx, ""]
    (dest / "README.md").write_text("\n".join(readme))
    print(f"assembled -> {dest}")


EXTRA = {"gates_single", "gates_complex1v1", "seed_robustness_table", "long_rna", "orphan_tests", "order_rna_single_chain",
         "rf_rna_single_chain_other_splits", "confidence_ptm", "near_floor", "novelty", "construction_stage_mcc"}

if __name__ == "__main__":
    main()
