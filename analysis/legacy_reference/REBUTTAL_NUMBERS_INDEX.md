# Rebuttal numbers -> reproducible script (NeurIPS 35782)

Every quantitative claim in `~/REBUTTAL_DRAFT.md` maps to a script here. All read from
`reproduce_rebuttal/results/benchmark_multiseed/*.db` (multiseed) + `results/paper_splits/` +
`results/long_run/`. Run in the `synEvo` conda env unless noted. SHS/MSA/no-MSA are ALL the multiseed ensemble now
(no-MSA = `ms('empty')` from `empty_*.db`, folded + eval'd 2026-07-28).

| rebuttal location | number(s) | script |
|---|---|---|
| General (1), uTY4-Q1, Vnro-W1 | orphan parity: single +0.019 CI[-0.024,+0.061]; RP +0.004 [-0.011,+0.018] | `multiseed_orphan_parity.py` |
| General (2), uTY4-Q2 | long-RNA ladder 0.816/0.599/0.232/0.220/0.129, 63% gap | `general2_long_run_ladder.py` |
| General (2), faN3-Q3 | runtime: search 87% of end-to-end, 73x; 8XTQ/7PTQ zero homologs | `general2_runtime.py` |
| General (2) | AF3 confidence ladder (pTM, single n=91): SHS 0.408 / MSA 0.421 / noMSA 0.337 / shuffle 0.302 / random 0.181; SHS-noMSA +0.071, SHS-shuffle +0.106, SHS-random +0.227 (all p<.001); pTM~RNA_TM spearman 0.90. Depth-controlled (shuffle/random BELOW noMSA). Confirmed on ranking_score + single-seed-matched | `confidence_ladder.py` (+ `confidence_ladder_confirm.py`) |
| General (3) | seed-robustness table, shift <=0.014 (SPOT-RNA single 0.014) | `robustness_table.py` |
| General (3) | across-sample stability | `across_samples.py` |
| uTY4-Q3 | denoise: AF3 +0.002 p=0.20, RhoFold+ -0.026 p<0.001, ref +0.082/+0.158 | `multiseed_claims_bcd.py` (B) |
| uTY4-Q4 | near-floor solvable orphans (0.45: MSA 0.611/SHS 0.590; 0.60: 0.638/0.661) | `uTY4_q4_near_floor.py` |
| uTY4-Q4, faN3-Q2 | 2D-MCC gate (GT) single+complex; self-consistency; ribozymes -- MULTISEED no-MSA floor: GT-gate no-MSA 0.594/0.521 +0.142(0.002)/+0.260; selfcons 0.501/0.467 +0.234/+0.235; single/complex bins; ribozyme SHS-noMSA +0.081, SHS-MSA +0.073 | `multiseed_nomsa_final.py` (+ `multiseed_nomsa_cis.py` for the CIs); superseded `faN3_q2_failure.py` (single-seed floor) |
| faN3-Q1 | family-leakage boost MULTISEED (+0.134 no-MSA / +0.226 MSA / +0.191 SHS); MSA-SHS gap +0.032->+0.067 | `multiseed_nomsa_final.py`; superseded `faN3_q1_leakage.py` (single-seed floor) |
| faN3-Q1 | sequence novelty (28/31 single, 111/128 complex); structure max-TM: single median 0.35 (9/31>0.5, 2/31>0.7); complex median 0.57 (94/128>0.5) -- complex structurally LESS novel (protein-fold-dominated), so we lean on sequence + leakage there | `faN3_q1_similarity.py` (RnaBench env for cd-hit; USalign for structure; results in results/faN3_q1/structure_maxTM_complex.json) |
| (supporting) | 4-seed vs published-seed benchmark gather | `gather_4seed_benchmark.py` |
| (supporting) | Wilcoxon on orphan deltas | `wilcoxon_orphans.py` |

Notes:
- `faN3_q1_similarity.py --mode structure --subtype complex` is the ~32k-pair US-align matrix (slow).
- no-MSA multiseed fold + eval DONE (2026-07-28); `multiseed_nomsa_final.py`/`_cis.py` reproduce all
  no-MSA-floor numbers now in the draft. The old `faN3_q2_failure.py`/`faN3_q1_leakage.py` still use the
  single-seed floor (kept for provenance; superseded).

## Construction / recall-not-precision (mechanism, camera-ready candidate)
| location | number(s) | script |
|---|---|---|
| 4P8Z/6WJS intervention (multiseed, 4 seeds) | recall->TM monotone->plateau; pk dip was seed artifact; SD<=0.041; 4P8Z .387/.452/.651/.663/.755/.781, 6WJS .737/.821/.842/.849/.855/.855 | `analysis/construction/eval_intervention_ms.py` -> `construction_recall_TM_multiseed.png` |
| 91-RNA construction (within-RNA) | TM(all)>TM(wc) 72/91 dTM+0.083 Wilcoxon p=9e-13; per-RNA Spearman median rho+0.40 (72/89); pooled rho 0.36 p=8e-13; nc/pk don't help alone, gain=topology completion | `analysis/construction/build_construct_91.py`+`build_inputs_91.sh`+`construct91_eval.py`+`construct91_within.py` -> `construct91_recall_TM.png` |
