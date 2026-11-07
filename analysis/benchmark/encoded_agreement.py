"""Agreement of AF3 single-chain predictions with the encoded RNAformer topology, per input condition.

Per target: mean over 4 seeds x 5 samples of the MCC between the RNAformer topology encoded in the SHS and the DSSR pairs
of the prediction (eval/selfcons, results/selfcons/<cond>_single_s*_m*.tsv.gz). Conditions: RNAformer SHS, no MSA,
equal-depth shuffled and random SHS, natural MSA. Paired Wilcoxon + 5000x bootstrap CI of the mean difference.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from paths import *  # noqa: F401,F403  (CRC locations, split_ids, db helpers)
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

CONDS = [("rnaformer_single", "SHS (RNAformer)"), ("empty_single", "no MSA"), ("shuffleSHS_single", "shuffled SHS"),
         ("randomSHS_single", "random SHS"), ("af3_single", "natural MSA")]
t = pd.DataFrame({lab: selfcons_means(tag) for tag, lab in CONDS})
rng = np.random.default_rng(0)
rows = []
for _, lab in CONDS:
    rows.append(dict(condition=lab, n=t[lab].notna().sum(), mean_mcc_vs_encoded=t[lab].mean()))
    if lab == "SHS (RNAformer)":
        continue
    x = t[["SHS (RNAformer)", lab]].dropna(); d = (x.iloc[:, 0] - x.iloc[:, 1]).values
    bs = [rng.choice(d, len(d)).mean() for _ in range(5000)]
    rows[-1].update(n_paired=len(d), d_SHS_minus=d.mean(), ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5),
                    wilcoxon_p=wilcoxon(x.iloc[:, 0], x.iloc[:, 1]).pvalue, n_SHS_higher=int((d > 0).sum()))
x = t[["shuffled SHS", "no MSA"]].dropna()
print(f"shuffled SHS vs no MSA: {x.iloc[:, 0].mean():.3f} vs {x.iloc[:, 1].mean():.3f} "
      f"(p={wilcoxon(x.iloc[:, 0], x.iloc[:, 1]).pvalue:.1e}, n={len(x)})")
r = pd.DataFrame(rows)
print(r.round(4).to_string(index=False))
r.to_csv(f"{OUT_ROOT}/encoded_agreement.csv", index=False)
t.to_csv(f"{OUT_ROOT}/encoded_agreement_per_target.csv", index_label="id")
