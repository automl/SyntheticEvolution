import pandas as pd
import numpy as np
import re
from scipy.stats import wilcoxon

# ---------- helpers ----------
def signif_stars(p):
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 5e-2: return '*'
    return 'ns'

def run_wilcoxon_with_means(df, label):
    d = df[['AF','RF']].dropna()
    n = len(d)
    if n < 1:
        return {'label': label, 'n': 0, 'af_mean': np.nan, 'rf_mean': np.nan,
                'stat': np.nan, 'p': np.nan, 'stars': 'ns', 'note': 'no paired data'}
    af_mean = float(d['AF'].mean()); rf_mean = float(d['RF'].mean())
    try:
        stat, p = wilcoxon(d['AF'], d['RF'], zero_method='pratt', alternative='two-sided')
        stars = signif_stars(p)
        return {'label': label, 'n': n, 'af_mean': af_mean, 'rf_mean': rf_mean,
                'stat': float(stat), 'p': float(p), 'stars': stars, 'note': ''}
    except ValueError as e:
        # e.g., all differences zero or too few non-zero diffs
        return {'label': label, 'n': n, 'af_mean': af_mean, 'rf_mean': rf_mean,
                'stat': 0.0, 'p': 1.0, 'stars': 'ns', 'note': str(e)}

def pretty_print_result(res, prefix='   '):
    note = f"  [{res['note']}]" if res.get('note') else ''
    print(
        f"{prefix}{res['label']:<22} n={res['n']:>4d}  "
        f"AF mean={res['af_mean']:.3f}  RF mean={res['rf_mean']:.3f}  "
        f"Wilcoxon: stat={res['stat']:.3f}  p={res['p']:.3e} {res['stars']}{note}"
    )

def build_paired(af_df, rf_df):
    """Inner-join AF and RF on exp_db_id; returns DataFrame with ['exp_db_id','AF','RF']"""
    left = af_df[['exp_db_id','Complex_RMSD']].rename(columns={'Complex_RMSD':'AF'})
    right = rf_df[['exp_db_id','Complex_RMSD']].rename(columns={'Complex_RMSD':'RF'})
    return left.merge(right, on='exp_db_id', how='inner')

# ---------- family status (robust, disjoint) ----------
_none_list_regex = re.compile(r'^\s*\[\s*"None"(?:\s*,\s*"None")*\s*\]\s*$')

def _is_no_family(val) -> bool:
    if pd.isna(val):
        return True
    s = str(val).strip()
    return s in {'[]','None','nan'} or bool(_none_list_regex.match(s))

def family_status_by_id(meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns one row per PDBId with boolean 'has_fam'.
    has_fam = True iff ANY row for that PDBId is a real family (not NaN/[]/["None",...]).
    """
    m = meta_df[['PDBId','RNAFamily']].copy()
    m['no_fam'] = m['RNAFamily'].apply(_is_no_family)
    status = m.groupby('PDBId', as_index=False)['no_fam'] \
              .apply(lambda s: ~s.all()) \
              .rename(columns={'no_fam':'has_fam'})
    return status

def split_and_test_by_family(pair_df: pd.DataFrame, meta_df: pd.DataFrame, label_prefix: str):
    """
    pair_df: columns ['exp_db_id','AF','RF']
    meta_df: table with PDBId, RNAFamily
    Prints all / has_fam / no_fam results and counts.
    """
    status = family_status_by_id(meta_df)
    pair = pair_df.merge(status, left_on='exp_db_id', right_on='PDBId', how='left') \
                  .drop(columns=['PDBId'])
    pair['has_fam'] = pair['has_fam'].fillna(False)

    # all
    res_all = run_wilcoxon_with_means(pair, f'{label_prefix} all')
    # has_fam
    res_has = run_wilcoxon_with_means(pair[pair['has_fam']], f'{label_prefix} has_fam')
    # no_fam
    res_no  = run_wilcoxon_with_means(pair[~pair['has_fam']], f'{label_prefix} no_fam')

    # sanity check
    n_all = len(pair); n_has = int(pair['has_fam'].sum()); n_no = n_all - n_has
    print(f"   counts: all={n_all}  has_fam={n_has}  no_fam={n_no}")

    pretty_print_result(res_all)
    pretty_print_result(res_has)
    pretty_print_result(res_no)

# ---------- load metadata (with RNAFamily) ----------
r_pdb = pd.read_csv('results/pdb_rna_rna_from_alphafold_eval.csv')
p_pdb = pd.read_csv('results/pdb_protein_rna_from_alphafold_eval.csv')

# ---------- load AF/RF result tables ----------
# All
af_all = pd.read_csv('results/csvs/All_alphafold.csv')
rf_all = pd.read_csv('results/csvs/All_rnaformer.csv')

# RNA / Protein selection for "all"
r_af_all = af_all[af_all['exp_db_id'].isin(r_pdb['PDBId'])]
r_rf_all = rf_all[rf_all['exp_db_id'].isin(r_pdb['PDBId'])]
p_af_all = af_all[af_all['exp_db_id'].isin(p_pdb['PDBId'])]
p_rf_all = rf_all[rf_all['exp_db_id'].isin(p_pdb['PDBId'])]

# a2021 subsets
r_af_a2021              = pd.read_csv('results/csvs/RNA_Monomers_a2021_alphafold.csv')
p_af_a2021              = pd.read_csv('results/csvs/RNA-Protein_a2021_alphafold.csv')
r_af_a2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_a2021_alphafold_orphan.csv')
p_af_a2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_a2021_alphafold_orphan.csv')
r_af_a2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_a2021_alphafold_non_orphan.csv')
p_af_a2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_a2021_alphafold_non_orphan.csv')

r_rf_a2021              = pd.read_csv('results/csvs/RNA_Monomers_a2021_rnaformer.csv')
p_rf_a2021              = pd.read_csv('results/csvs/RNA-Protein_a2021_rnaformer.csv')
r_rf_a2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_a2021_rnaformer_orphan.csv')
p_rf_a2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_a2021_rnaformer_orphan.csv')
r_rf_a2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_a2021_rnaformer_non_orphan.csv')
p_rf_a2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_a2021_rnaformer_non_orphan.csv')

# b2021 subsets
r_af_b2021              = pd.read_csv('results/csvs/RNA_Monomers_b2021_alphafold.csv')
p_af_b2021              = pd.read_csv('results/csvs/RNA-Protein_b2021_alphafold.csv')
r_af_b2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_b2021_alphafold_orphan.csv')
p_af_b2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_b2021_alphafold_orphan.csv')
r_af_b2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_b2021_alphafold_non_orphan.csv')
p_af_b2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_b2021_alphafold_non_orphan.csv')

r_rf_b2021              = pd.read_csv('results/csvs/RNA_Monomers_b2021_rnaformer.csv')
p_rf_b2021              = pd.read_csv('results/csvs/RNA-Protein_b2021_rnaformer.csv')
r_rf_b2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_b2021_rnaformer_orphan.csv')
p_rf_b2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_b2021_rnaformer_orphan.csv')
r_rf_b2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_b2021_rnaformer_non_orphan.csv')
p_rf_b2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_b2021_rnaformer_non_orphan.csv')

# ---------- paired models (includes b2021) ----------
paired_models = [
    ('all',               r_af_all,              r_rf_all,              p_af_all,              p_rf_all),
    ('a2021',             r_af_a2021,            r_rf_a2021,            p_af_a2021,            p_rf_a2021),
    ('a2021_orphan',      r_af_a2021_orphan,     r_rf_a2021_orphan,     p_af_a2021_orphan,     p_rf_a2021_orphan),
    ('a2021_non_orphan',  r_af_a2021_non_orphan, r_rf_a2021_non_orphan, p_af_a2021_non_orphan, p_rf_a2021_non_orphan),
    ('b2021',             r_af_b2021,            r_rf_b2021,            p_af_b2021,            p_rf_b2021),
    ('b2021_orphan',      r_af_b2021_orphan,     r_rf_b2021_orphan,     p_af_b2021_orphan,     p_rf_b2021_orphan),
    ('b2021_non_orphan',  r_af_b2021_non_orphan, r_rf_b2021_non_orphan, p_af_b2021_non_orphan, p_rf_b2021_non_orphan),
]

for subset_name, r_af_df, r_rf_df, p_af_df, p_rf_df in paired_models:
    print(f"\nSubset: {subset_name}")

    # RNA
    r_pair = build_paired(r_af_df, r_rf_df)
    print(" - RNA AF vs RF")
    split_and_test_by_family(r_pair, r_pdb, 'RNA')

    # Protein
    p_pair = build_paired(p_af_df, p_rf_df)
    print(" - Protein AF vs RF")
    split_and_test_by_family(p_pair, p_pdb, 'Protein')

# ---------- POOLED ANALYSIS: ALL a2021 (RNA + Protein together) ----------
print("\n=== Pooled analysis: ALL a2021 (RNA + Protein together) ===")
# build paired a2021 tables
r_pair_a2021 = build_paired(r_af_a2021, r_rf_a2021)
p_pair_a2021 = build_paired(p_af_a2021, p_rf_a2021)
# concatenate RNA and Protein pairs
pooled_a2021 = pd.concat([r_pair_a2021, p_pair_a2021], ignore_index=True)

# 1) overall (pooled) Wilcoxon + means
overall_res = run_wilcoxon_with_means(pooled_a2021, 'a2021 pooled (all)')
pretty_print_result(overall_res, prefix='   ')

# 2) split by has_fam / no_fam using combined metadata
combined_meta = pd.concat(
    [r_pdb[['PDBId','RNAFamily']], p_pdb[['PDBId','RNAFamily']]],
    ignore_index=True
)
print(" - Split by family status (pooled)")
split_and_test_by_family(pooled_a2021, combined_meta, 'a2021 pooled')

# ---------- list families & compare a2021 vs b2021 ----------
def normalize_family_values(series: pd.Series) -> pd.Series:
    """Return series with NaNs and '["None", ...]' normalized to NaN for clean unique() display."""
    def norm_one(x):
        if _is_no_family(x):
            return np.nan
        return str(x)
    return series.apply(norm_one)

def families_in_subset(subset_df: pd.DataFrame, meta_df: pd.DataFrame) -> set:
    """Families present among exp_db_id of subset_df, using meta_df[PDBId,RNAFamily]."""
    ids = set(subset_df['exp_db_id'])
    fam = meta_df[meta_df['PDBId'].isin(ids)]['RNAFamily']
    fam = normalize_family_values(fam).dropna().unique()
    return set(fam)

print("\n=== Unique families (overall, from metadata tables) ===")
rna_fams_overall = normalize_family_values(r_pdb['RNAFamily']).dropna().unique()
prot_fams_overall = normalize_family_values(p_pdb['RNAFamily']).dropna().unique()
print(f"RNA families overall   ({len(rna_fams_overall)}): {sorted(rna_fams_overall)}")
print(f"Protein families overall({len(prot_fams_overall)}): {sorted(prot_fams_overall)}")

print("\n=== Families present in subsets: a2021 vs b2021 (RNA) ===")
rna_fam_a2021 = families_in_subset(r_af_a2021, r_pdb)
rna_fam_b2021 = families_in_subset(r_af_b2021, r_pdb)
print(f"a2021 RNA families ({len(rna_fam_a2021)}): {sorted(rna_fam_a2021)}")
print(f"b2021 RNA families ({len(rna_fam_b2021)}): {sorted(rna_fam_b2021)}")
print(f"Missing in b2021 (vs a2021): {sorted(rna_fam_a2021 - rna_fam_b2021)}")
print(f"New in b2021 (vs a2021):     {sorted(rna_fam_b2021 - rna_fam_a2021)}")

print("\n=== Families present in subsets: a2021 vs b2021 (Protein) ===")
prot_fam_a2021 = families_in_subset(p_af_a2021, p_pdb)
prot_fam_b2021 = families_in_subset(p_af_b2021, p_pdb)
print(f"a2021 Protein families ({len(prot_fam_a2021)}): {sorted(prot_fam_a2021)}")
print(f"b2021 Protein families ({len(prot_fam_b2021)}): {sorted(prot_fam_b2021)}")
print(f"Missing in b2021 (vs a2021): {sorted(prot_fam_a2021 - prot_fam_b2021)}")
print(f"New in b2021 (vs a2021):     {sorted(prot_fam_b2021 - prot_fam_a2021)}")
