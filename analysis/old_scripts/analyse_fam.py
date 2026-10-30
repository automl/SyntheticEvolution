import pandas as pd
import numpy as np
import re
import ast
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
    # left = af_df[['exp_db_id','Complex_RMSD']].rename(columns={'Complex_RMSD':'AF'})
    # right = rf_df[['exp_db_id','Complex_RMSD']].rename(columns={'Complex_RMSD':'RF'})
    left = af_df[['exp_db_id','Complex_LDDT']].rename(columns={'Complex_LDDT':'AF'})
    right = rf_df[['exp_db_id','Complex_LDDT']].rename(columns={'Complex_LDDT':'RF'})

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

r_rf_a2021              = pd.read_csv('results/csvs/RNA_Monomers_a2021_rnaformerN100.csv')
p_rf_a2021              = pd.read_csv('results/csvs/RNA-Protein_a2021_rnaformerN100.csv')
r_rf_a2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_a2021_rnaformerN100_orphan.csv')
p_rf_a2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_a2021_rnaformerN100_orphan.csv')
r_rf_a2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_a2021_rnaformerN100_non_orphan.csv')
p_rf_a2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_a2021_rnaformerN100_non_orphan.csv')

# b2021 subsets
r_af_b2021              = pd.read_csv('results/csvs/RNA_Monomers_b2021_alphafold.csv')
p_af_b2021              = pd.read_csv('results/csvs/RNA-Protein_b2021_alphafold.csv')
r_af_b2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_b2021_alphafold_orphan.csv')
p_af_b2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_b2021_alphafold_orphan.csv')
r_af_b2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_b2021_alphafold_non_orphan.csv')
p_af_b2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_b2021_alphafold_non_orphan.csv')

r_rf_b2021              = pd.read_csv('results/csvs/RNA_Monomers_b2021_rnaformerN100.csv')
p_rf_b2021              = pd.read_csv('results/csvs/RNA-Protein_b2021_rnaformerN100.csv')
r_rf_b2021_orphan       = pd.read_csv('results/csvs/RNA_Monomers_b2021_rnaformerN100_orphan.csv')
p_rf_b2021_orphan       = pd.read_csv('results/csvs/RNA-Protein_b2021_rnaformerN100_orphan.csv')
r_rf_b2021_non_orphan   = pd.read_csv('results/csvs/RNA_Monomers_b2021_rnaformerN100_non_orphan.csv')
p_rf_b2021_non_orphan   = pd.read_csv('results/csvs/RNA-Protein_b2021_rnaformerN100_non_orphan.csv')

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

val_ids = ['3BWP', '255D', '6E80', '4FAQ', '3Q50', '2ZY6', '4E8V', '7KD1', '3AM1', '2Q1R', '3GCA', '3ND3', '3DHS', '3NPN', '7M5O', '4E8P', '4E8Q', '4RBQ', '1U9S', '4WJ4', '6UES', '4EN5', '4E8M', '4C40', '6TF3', '5C5W', '4CS1', '4E8N', '5DA6', '6TB7', '4P8Z', '2A2E', '6IV9', '2A64', '5HSW', '413D', '3R4F', '2DVI', '4GMA', '6TFE', '3D0M', '4DS6', '387D', '7D7W', '6TF1', '6UET', '6T3S', '6DTD', '6PQ7', '4AOB']

paired_models_copy = []

for subset_name, r_af_df, r_rf_df, p_af_df, p_rf_df in paired_models:
    # Filter out validation IDs
    r_af_df_filt = r_af_df[~r_af_df['exp_db_id'].isin(val_ids)]
    r_rf_df_filt = r_rf_df[~r_rf_df['exp_db_id'].isin(val_ids)]
    p_af_df_filt = p_af_df[~p_af_df['exp_db_id'].isin(val_ids)]
    p_rf_df_filt = p_rf_df[~p_rf_df['exp_db_id'].isin(val_ids)]

    paired_models_copy.append((subset_name, r_af_df_filt, r_rf_df_filt, p_af_df_filt, p_rf_df_filt))

paired_models = paired_models_copy

for subset_name, r_af_df, r_rf_df, p_af_df, p_rf_df in paired_models:
    print(f"\nSubset: {subset_name}")


    r_af_df = r_af_df[~r_af_df['exp_db_id'].isin(val_ids)]
    r_rf_df = r_rf_df[~r_rf_df['exp_db_id'].isin(val_ids)]
    p_af_df = p_af_df[~p_af_df['exp_db_id'].isin(val_ids)]
    p_rf_df = p_rf_df[~p_rf_df['exp_db_id'].isin(val_ids)]

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

def compare_af_rf_by_family(pair_df: pd.DataFrame,
                            meta_df: pd.DataFrame,
                            label_prefix: str):
    """
    For a given paired AF/RF table and metadata with RNAFamily,
    run AF vs RF Wilcoxon per RNAFamily and print the results.
    """
    # Normalize family values (re-use your helper)
    meta = meta_df[['PDBId', 'RNAFamily']].copy()
    meta['RNAFamily_norm'] = normalize_family_values(meta['RNAFamily'])

    # Keep only entries with a real family annotation
    meta = meta.dropna(subset=['RNAFamily_norm']).drop_duplicates(['PDBId', 'RNAFamily_norm'])

    merged = pair_df.merge(meta[['PDBId', 'RNAFamily_norm']],
                           left_on='exp_db_id',
                           right_on='PDBId',
                           how='left') \
                    .drop(columns=['PDBId'])

    merged = merged.dropna(subset=['RNAFamily_norm'])

    if merged.empty:
        print(f"   [No annotated families found for {label_prefix}]")
        return

    # Group by family and run Wilcoxon AF vs RF
    print(f"   Per-family AF vs RF comparisons for {label_prefix}:")
    # Sort families by decreasing sample size for nicer output
    grouped = merged.groupby('RNAFamily_norm')
    fam_sizes = grouped.size().sort_values(ascending=False)

    for fam_name in fam_sizes.index:
        fam_df = grouped.get_group(fam_name)
        label = f"{label_prefix} fam={fam_name}"
        res = run_wilcoxon_with_means(fam_df, label)
        pretty_print_result(res, prefix='      ')


print("\n=== Unique families (overall, from metadata tables) ===")
rna_fams_overall = normalize_family_values(r_pdb['RNAFamily']).dropna().unique()
prot_fams_overall = normalize_family_values(p_pdb['RNAFamily']).dropna().unique()
print(f"RNA families overall   ({len(rna_fams_overall)}): {sorted(rna_fams_overall)}")
print(f"Protein families overall({len(prot_fams_overall)}): {sorted(prot_fams_overall)}")

all_fam = pd.concat([r_pdb[['PDBId','RNAFamily']], p_pdb[['PDBId','RNAFamily']]], ignore_index=True)
all_fam = all_fam[all_fam['PDBId'].isin(pd.concat([paired_models[0][1], paired_models[0][3]])['exp_db_id'])]
rna_fams = all_fam[all_fam['PDBId'].isin(paired_models[0][1]['exp_db_id'])]
prot_fams = all_fam[all_fam['PDBId'].isin(paired_models[0][3]['exp_db_id'])]
all_fam['RNAFamily_norm'] = normalize_family_values(all_fam['RNAFamily'])

print(f"All families overall: {len(all_fam)}")
print('Individual family counts:')
fam_counts = all_fam.groupby('RNAFamily').size().sort_values(ascending=False)
print(all_fam['RNAFamily'])
print(all_fam.groupby('RNAFamily').count())
print('=== Families present in RNA & Protein and their counts ===')
for fam_name, count in fam_counts.items():
    print(f"   {fam_name}: {count}")

print("\n=== Families present in RNA and their counts ===")
for fam_name, df in rna_fams.groupby('RNAFamily'):
    af_mean, af_std = paired_models[0][1][paired_models[0][1]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].mean(), paired_models[0][1][paired_models[0][1]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].std()
    rf_mean, rf_std = paired_models[0][2][paired_models[0][2]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].mean(), paired_models[0][2][paired_models[0][2]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].std()
    print(f"   {fam_name}: {df.shape[0]}, AFmean: {af_mean:.3f}, AFstd: {af_std:.3f}, RFmean: {rf_mean:.3f}, RFstd: {rf_std:.3f}")

print("\n=== Families present in Protein and their counts ===")
for fam_name, df in prot_fams.groupby('RNAFamily'):
    af_mean, af_std = paired_models[0][3][paired_models[0][3]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].mean(), paired_models[0][3][paired_models[0][3]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].std()
    rf_mean, rf_std = paired_models[0][4][paired_models[0][4]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].mean(), paired_models[0][4][paired_models[0][4]['exp_db_id'].isin(df['PDBId'])]['Complex_LDDT'].std()
    print(f"   {fam_name}: {df.shape[0]}, AFmean: {af_mean:.3f}, AFstd: {af_std:.3f}, RFmean: {rf_mean:.3f}, RFstd: {rf_std:.3f}")

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

print("\n=== Per-family AF vs RF comparisons for each subset (RNA & Protein) ===")

for subset_name, r_af_df, r_rf_df, p_af_df, p_rf_df in paired_models:
    print(f"\nSubset: {subset_name}")

    # Use the same validation IDs filter as above
    val_ids = ['3BWP', '255D', '6E80', '4FAQ', '3Q50', '2ZY6', '4E8V', '7KD1',
               '3AM1', '2Q1R', '3GCA', '3ND3', '3DHS', '3NPN', '7M5O', '4E8P',
               '4E8Q', '4RBQ', '1U9S', '4WJ4', '6UES', '4EN5', '4E8M', '4C40',
               '6TF3', '5C5W', '4CS1', '4E8N', '5DA6', '6TB7', '4P8Z', '2A2E',
               '6IV9', '2A64', '5HSW', '413D', '3R4F', '2DVI', '4GMA', '6TFE',
               '3D0M', '4DS6', '387D', '7D7W', '6TF1', '6UET', '6T3S', '6DTD',
               '6PQ7', '4AOB']

    # Filter out validation IDs (same as in the main subset loop)
    r_af_sub = r_af_df[~r_af_df['exp_db_id'].isin(val_ids)]
    r_rf_sub = r_rf_df[~r_rf_df['exp_db_id'].isin(val_ids)]
    p_af_sub = p_af_df[~p_af_df['exp_db_id'].isin(val_ids)]
    p_rf_sub = p_rf_df[~p_rf_df['exp_db_id'].isin(val_ids)]

    # Build paired tables
    r_pair_sub = build_paired(r_af_sub, r_rf_sub)
    p_pair_sub = build_paired(p_af_sub, p_rf_sub)

    # RNA per-family comparison
    print(" - RNA per-family")
    compare_af_rf_by_family(r_pair_sub, r_pdb, f"{subset_name} RNA")

    # Protein per-family comparison
    print(" - Protein per-family")
    compare_af_rf_by_family(p_pair_sub, p_pdb, f"{subset_name} Protein")
