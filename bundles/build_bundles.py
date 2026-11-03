#!/usr/bin/env python3
"""Assemble the two download bundles (outputs zipped slim: see bundles/zip_bundle.sh) as symlink trees (zip follows symlinks, so no data is duplicated).

  2026_neurips_shs_inputs/   AF3 input JSONs (MSAs/SHS/templates embedded) + reference mmCIFs
  2026_neurips_shs_outputs/  raw AF3 outputs (all seeds/samples, confidences) + RhoFold+ outputs + metric DBs

Only targets of the paper (91 single chains, 377 complexes, 25 long RNAs, the case studies) are included.
Input JSONs are given once per condition with their multiseed `modelSeeds`; the seed-1 variants differ only in
`modelSeeds` (fold/prepare_inputs.py --seeds 1). Each bundle gets MANIFEST.tsv (bundle path -> source path).

  python bundles/build_bundles.py --which inputs|outputs --dest <dir>

Authors only: reads the original layout (AF_INPUT, AF_OUTPUT, GT_*, NOMSA_LONG_INPUT, RHOFOLD_DIR from config.local.sh).
"""
import argparse, csv, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
from paths import AF_INPUT, AF_OUTPUT, CRC, DB_DIR, DB_LONG, DB_RHOFOLD, GT_COMPLEX, GT_LONG, GT_SINGLE, IDS_DIR  # noqa: E402

if not (AF_INPUT and AF_OUTPUT):
    sys.exit("set AF_INPUT/AF_OUTPUT (original layout) in config.local.sh")
IDS = {k: {l.strip().lower() for l in open(IDS_DIR / f) if l.strip()}
       for k, f in (("single", "single_91.txt"), ("complex", "complex_377.txt"), ("long", "long_25.txt"))}
# multiseed output dir -> input dir it was folded from (bind mounts of the original fold scripts)
MS_INPUT = {
    "nature_methods_full_evaluation_data_files_af3_rerun_multiseeds_predictions": "nature_methods_full_evaluation_data_files_af3_rerun_multiseeds",
    "all_data_nature_methods_synMSA_config80_lsHPO_empty_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_empty_multiseed/config80",
    "all_data_nature_methods_synMSA_config80_lsHPO_rnafoldN100_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_N100_rnafold_multiseeds/config80",
    "all_data_nature_methods_synMSA_config80_lsHPO_spotrnaN100_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_N100_spotrna_multiseeds/config80",
    "all_data_nature_methods_synMSA_config80_lsHPO_rnaformerN100_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_N100_rnaformer_multiseed/config80",
    "all_data_nature_methods_synMSA_config80_lsHPO_dssrN100_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_N100_dssr_multiseeds/config80",
    "rnamonomers_denoise_data_nature_methods_synMSA_config80_lsHPO_N100_rnafold_multiseed": "rnamonomers_denoise_data_nature_methods_synMSA_config80_lsHPO_N100_rnafold_multiseed/config80",
    **{f"rnaonly_data_nature_methods_{k}_config80_lsHPO_N100_rnaformer_multiseed": f"rnaonly_data_nature_methods_{k}_config80_lsHPO_N100_rnaformer_multiseed/config80"
       for k in ("randomSHS", "shuffleSHS", "nat_and_random", "nat_and_shuffle")},
    **{f"rnaonly_data_nature_methods_nat_and_synMSA_config80_lsHPO_N{n}_rnaformer_multiseed": f"rnaonly_data_nature_methods_nat_and_synMSA_config80_lsHPO_N{n}_rnaformer_multiseed/config80"
       for n in (20, 50, 100, 200, 500, 1000)},
    **{f"rnaonly_data_nature_methods_synMSA_config80_lsHPO_N{n}_rnaformer_multiseed": f"rnaonly_data_nature_methods_synMSA_config80_lsHPO_N{n}_rnaformer_multiseed/config80"
       for n in (5000, 10000, 15000, 20000)},
    "crc2026_complex_N5000_rnaformer_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_N5000_rnaformer_multiseeds",
    "crc2026_complex_N20000_rnaformer_multiseed": "all_data_nature_methods_synMSA_config80_lsHPO_N20000_rnaformer_multiseeds",
    "crc2026_complex_N10000_rnaformer_multiseed": "crc2026_complex_N10000_rnaformer_multiseed",
    # depth completion (fold/runs/submit_crc_depth_folds.sh)
    **{f"all_data_nature_methods_synMSA_config80_lsHPO_rnaformerN{n}_multiseed": f"all_data_nature_methods_synMSA_config80_lsHPO_N{n}_rnaformer_multiseeds" for n in (20, 500)},
    "crc2026_N1000_rnaformer_multiseed": "crc2026_N1000_rnaformer_multiseed",
    "crc2026_complex_N15000_rnaformer_multiseed": "crc2026_complex_N15000_rnaformer_multiseed",
}
LONG = {  # condition -> (seed-1 output, multiseed output, multiseed input)
    "af3msa": ("a2021_rnaonly_predictions+b2021_rnaonly_predictions", "neurips_af3msa_rna_long_runs_multiseed", f"{AF_INPUT}/neurips_rna_long_runs_af3msa_multiseed/config80"),
    "dssr": ("neurips_dssr_rna_long_runs", "neurips_dssr_rna_long_runs_multiseed", f"{AF_INPUT}/neurips_rna_long_runs_dssr_multiseed/config80"),
    "nomsa": ("neurips_nomsa_rna_long_runs", "neurips_nomsa_rna_long_runs_multiseed", os.environ.get("NOMSA_LONG_INPUT", "")),
    "random": ("neurips_random_rna_long_runs", "neurips_random_rna_long_runs_multiseed", f"{AF_INPUT}/neurips_rna_long_runs_random_dssr_multiseed/config80"),
    "shuffle": ("neurips_shuffle_rna_long_runs", "neurips_shuffle_rna_long_runs_multiseed", f"{AF_INPUT}/neurips_rna_long_runs_shuffle_dssr_multiseed/config80"),
}
CASE_OUT = ["4P8Z_intervention", "6WJS_intervention", "4P8Z_intervention_full_sweep", "6WJS_intervention_full_sweep",
            "4P8Z_intervention_full_sweep_multiseed", "6WJS_intervention_full_sweep_multiseed", "construct91_multiseed",
            "construct91_N100_multiseed"]
CASE_IN = ["4P8Z_intervention", "6WJS_intervention", "4P8Z_intervention_data_config80_4p8z_intervention",
           "4P8Z_intervention_data_config80_4p8z_intervention_multiseed", "6WJS_intervention_scaling_sweep_config80_6wjs_intervention",
           "6WJS_intervention_scaling_sweep_config80_6wjs_intervention_multiseed", "construct91_multiseed",
           "construct91_N100_multiseed"]
RHOFOLD = os.environ.get("RHOFOLD_DIR", "")


class Tree:
    def __init__(self, dest):
        self.dest = Path(dest); self.rows = []
        self.dest.mkdir(parents=True, exist_ok=True)

    def link(self, rel, src):
        src = Path(src); p = self.dest / rel
        if not src.exists(): print("  missing source:", src); return
        p.parent.mkdir(parents=True, exist_ok=True)
        if not p.is_symlink(): os.symlink(src, p)
        self.rows.append((rel, str(src)))

    def targets(self, rel, src_dir, ids):
        """link every entry of src_dir whose name starts with a target id (AF3 output dirs / input subdirs)."""
        n = 0
        for e in sorted(Path(src_dir).iterdir()) if Path(src_dir).is_dir() else []:
            if e.name[:4].lower() in ids and (e.is_dir() and (e.name[4:5] in ("_", "") )):
                self.link(f"{rel}/{e.name}", e); n += 1
        return n

    def inputs(self, rel, src_dir, ids):
        """input roots are <idx>/<PDBID>...json or <PDBID>/<json>: link the json as <rel>/<PDBID>/<name>."""
        n = 0
        for j in sorted(Path(src_dir).glob("*/*.json")):
            if j.name[:4].lower() in ids:
                self.link(f"{rel}/{j.name[:4].upper()}/{j.name}", j); n += 1
        return n

    def close(self, readme):
        with open(self.dest / "MANIFEST.tsv", "w", newline="") as f:
            csv.writer(f, delimiter="\t", lineterminator="\n").writerows([("bundle_path", "source")] + self.rows)
        (self.dest / "README.txt").write_text(readme)
        print(f"{self.dest}: {len(self.rows)} entries")


def outputs(dest):
    t = Tree(dest)
    for c in csv.DictReader(open(CRC / "eval" / "conditions.tsv"), delimiter="\t"):
        name = f'{c["tag"]}_{c["idset"]}'; ids = IDS[c["idset"]]
        if c["seed1_dir"] == c["multiseed_dir"]:
            n = t.targets(f"af3/{name}/seeds_1_0_42_137", f"{AF_OUTPUT}/{c['multiseed_dir']}", ids)
        else:
            n = t.targets(f"af3/{name}/seed_1", f"{AF_OUTPUT}/{c['seed1_dir']}", ids)
            n += t.targets(f"af3/{name}/seeds_0_42_137", f"{AF_OUTPUT}/{c['multiseed_dir']}", ids)
        print(f"  {name}: {n} target dirs")
    for cond, (s1, ms, _) in LONG.items():
        for d in s1.split("+"):
            t.targets(f"af3_long_rna/{cond}/seed_1", f"{AF_OUTPUT}/{d}", IDS["long"])
        t.targets(f"af3_long_rna/{cond}/seeds_0_42_137", f"{AF_OUTPUT}/{ms}", IDS["long"])
    for d in CASE_OUT:
        t.link(f"af3_case_studies/{d}", f"{AF_OUTPUT}/{d}")
    for d in sorted(Path(RHOFOLD).glob("predictions_rhofold_*")):   # 5 SHS variants + native MSA (stdmsa)
        t.link(f"rhofold/{d.name}", d)
    t.link("metric_dbs/benchmark", DB_DIR)
    t.link("metric_dbs/long_rna", DB_LONG)
    t.link("metric_dbs/rhofold", DB_RHOFOLD)
    t.close("""2026_neurips_shs_outputs -- raw predictions and metric databases for the NeurIPS 2026 SHS paper.
af3/<condition>_<single|complex>/{seed_1,seeds_0_42_137}/<target>/  AlphaFold 3 outputs: seed-S_sample-K/model.cif and
    summary_confidences.json (pTM/ipTM/ranking score) for every seed/sample, ranking_scores.csv. The full per-sample
    confidences.json (PAE / per-atom pLDDT) and the input copies (*_data.json; see the inputs bundle) are not included.
    AlphaFold 3 outputs are subject to the AlphaFold 3 Output Terms of Use (TERMS_OF_USE.md in each dir).
af3_long_rna/, af3_case_studies/   long-RNA (>200 nt) ladder, 4P8Z/6WJS interventions, 91-RNA construction sweep.
rhofold/   RhoFold+ predictions, one deterministic run per input: SHS from RNAfold / SPOT-RNA / RNAformer, RNAfold-denoised,
    reference topology (DSSR), and the native RhoFold+ MSA pipeline (stdmsa); relaxed_1000_model.pdb is scored.
metric_dbs/   one SQLite DB per <condition>_<idset>_s<seed>_m<sample> (tables pred_*/exp_*), produced by eval/run_eval_array.sh:
    benchmark/ (AF3, 4 seeds x 5 samples), long_rna/, rhofold/ (<tag>_single_s1_m0.db).
Condition tags are explained in eval/conditions.tsv of the code repository. MANIFEST.tsv lists the source of every entry.
""")


def inputs(dest):
    t = Tree(dest)
    for rel, src in (("ground_truth/single_91", GT_SINGLE), ("ground_truth/complex_377", GT_COMPLEX), ("ground_truth/long_25", GT_LONG)):
        for f in sorted(Path(src).glob("*.cif")):
            t.link(f"{rel}/{f.name}", f.resolve())
    for f in sorted((CRC / "data" / "gt_intervention").glob("*.cif")):
        t.link(f"ground_truth/case_studies/{f.name}", f)
    seen = set()
    for c in csv.DictReader(open(CRC / "eval" / "conditions.tsv"), delimiter="\t"):
        src = MS_INPUT[c["multiseed_dir"]]; name = f'{c["tag"]}_{c["idset"]}'
        n = t.inputs(f"af3_inputs/{name}", f"{AF_INPUT}/{src}", IDS[c["idset"]])
        print(f"  {name}: {n} input JSONs  <- {src}")
    for cond, (_, _, src) in LONG.items():
        t.inputs(f"af3_inputs_long_rna/{cond}", src, IDS["long"])
    for d in CASE_IN:
        t.link(f"af3_inputs_case_studies/{d}", f"{AF_INPUT}/{d}")
    for d in sorted(Path(RHOFOLD, "inputs").iterdir()) if RHOFOLD else []:   # RhoFold+ A3M/FASTA inputs
        t.link(f"rhofold_inputs/{d.name}", d)
    if RHOFOLD:   # how RhoFold+ was run (pipeline scripts + README of the RhoFold+ runs)
        t.link("rhofold_inputs/run_scripts", Path(RHOFOLD, "scripts_dominik"))
    t.close("""2026_neurips_shs_inputs -- AlphaFold 3 inputs and reference structures for the NeurIPS 2026 SHS paper.
af3_inputs/<condition>_<single|complex>/<PDBID>/*.json  AF3 input JSONs with the RNA MSA slot holding natural MSAs, SHS,
    or controls (protein MSAs/templates unchanged). `modelSeeds` holds the multiseed seeds; seed-1 variants differ only in
    `modelSeeds`. Fold with --run_data_pipeline=False.
af3_inputs_long_rna/, af3_inputs_case_studies/   long-RNA ladder, 4P8Z/6WJS interventions, construction sweep (N=1000, N=100).
rhofold_inputs/   RhoFold+ inputs: prepared_shs_inputs_from_<variant>_folder/{a3m,fasta} (same SHS as the AF3 inputs, A3M form)
    and prepared_db_msa_inputs/ (FASTA for the native RhoFold+ MSA pipeline); run_scripts/ = the RhoFold+ run scripts.
ground_truth/   reference mmCIF files (PDB-derived) used for scoring.
MANIFEST.tsv lists the source of every entry.
""")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--which", choices=["inputs", "outputs"], required=True)
    ap.add_argument("--dest", required=True); a = ap.parse_args()
    (inputs if a.which == "inputs" else outputs)(a.dest)
