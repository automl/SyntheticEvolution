import sqlite3

db_path = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/database/rbpDatabase.db"

columns_to_keep = [
    "exp_db_id", "FileName", "ProteinSequence",
    "ProteinChainIDs", "ProteinLength",
    "NumberProteins", "AAMotif", "AAC", "AAPproteinRNA", "RNASequence", "RNAChainIDs",
    "RNALength", "NumberRNAs", "RNAMotif", "RNAMotifLength", "RNAModified", "ContactList", "ChainIDpairList_proteinRNA",
    "Hbond_proteinRNA", "vdWbond_proteinRNA", "ProteinRNAInterfaceArea", "ProteinRNAInterfaceRatio", "Free_energy",
    "Binding_affinity_kd", "Seed", "Model", "af3_protein_pTM", "af3_rna_pTM", "af3_protein_ipTM", "af3_rna_ipTM",
    "af3_protein_pLDDT_avg", "af3_rna_pLDDT_avg", "af3_global_pae_avg", "af3_chain_pair_pae_min",
    "af3_fraction_disordered",
    "af3_has_clash", "af3_ranking_score", "af3_protein_MSA", "af3_rna_MSA", "rna_prot_interface_atom_ids",
    "interface_rna_atom_ids	All_RMSD", "Protein_RMSD", "RNA_RMSD", "Protein_LDDT", "RNA_LDDT", "Protein_TM", "RNA_TM",
    "Protein_GDT_TS", "RNA_GDT_TS", "RNA_INF", "RNA_DI", "RNA_f1_score", "RNA_precision", "RNA_recall", "RNA_mcc",
    "RNA_wl", "RNA_GlycosidicBond", "RNA_SugarPucker", "RNA_GammaAngle", "RNAmotif_score", "RNA_Stems",
    "RNA_HairpinLoops", "RNA_InternalLoops", "RNA_MultibranchLoops", "RNA_DanglingEnds", "RNA_Pseudoknots",
    "RNA_isDuplex", "RNAMotif_isDuplex", "RNA_ElectrostaticPotential", "RNA_GC_Content", "RNA_SequenceComplexity",
    "RNA_Shannon_Entropy_K1", "RNA_Shannon_Entropy_K2", "RNA_GC_Skew", "RNA_AU_Skew", "RNA_AU_Content", "RNA_MFE_Value",
    "RNA_MFE_Structure", "RNA_Ensemble_Diversity", "RNA_Centroid_Energy", "RNA_msa_size", "Protein_msa_size"
]

conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Step 1: Create new table with only desired columns
columns_str = ", ".join(columns_to_keep)
cursor.execute(f"""
    CREATE TABLE pred_new AS
    SELECT {columns_str} FROM pred_protein_rna
""")

# Step 2: Drop old table
cursor.execute("DROP TABLE pred_protein_rna")

# Step 3: Rename new table
cursor.execute("ALTER TABLE pred_new RENAME TO pred_protein_rna")

conn.commit()
conn.close()
