import sqlite3
import ast

# === CONFIGURATION ===
DELETE_IF = {
    "NumberProteins_gt": None,      # Delete if NumberProteins > this value (set to None to disable)
    "NumberRNAs_not": None,         # Delete if NumberRNAs != this value (set to None to disable)
    "AAmatch_score_not": None,      # Delete if AAmatch_score != this value (set to None to disable)
    "Ions_empty": True,          # Delete if Ions is not '[]' (set to False to disable)
    "Ligands_empty": None,       # Delete if Ligands is not '[]' (set to False to disable)
}
DB_PATH = "/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/database/rbpDatabase.db"

# === END CONFIGURATION ===

conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()
exp_ids_to_delete = set()

# Build WHERE clauses for pred_protein_rna
pred_clauses = []
if DELETE_IF["NumberProteins_gt"] is not None:
    pred_clauses.append(f"NumberProteins > {DELETE_IF['NumberProteins_gt']}")
if DELETE_IF["NumberRNAs_not"] is not None:
    pred_clauses.append(f"NumberRNAs != {DELETE_IF['NumberRNAs_not']}")
if DELETE_IF["AAmatch_score_not"] is not None:
    pred_clauses.append(f"AAmatch_score != {DELETE_IF['AAmatch_score_not']}")

# Collect exp_db_ids to delete from pred_protein_rna
if pred_clauses:
    where_pred = " OR ".join(f"({c})" for c in pred_clauses)
    cursor.execute(f"SELECT exp_db_id FROM pred_rna_rna WHERE {where_pred}")
    exp_ids_to_delete.update(row[0] for row in cursor.fetchall())
    cursor.execute(f"DELETE FROM pred_rna_rna WHERE {where_pred}")

# Build WHERE clauses for exp_protein_rna
exp_clauses = []
if DELETE_IF["Ions_empty"]:
    exp_clauses.append("(Ions IS NULL OR TRIM(Ions) == '[]')")
if DELETE_IF["Ligands_empty"]:
    exp_clauses.append("(Ligands IS NULL OR TRIM(Ligands) == '[]')")

# Collect PDBIds to delete from exp_protein_rna
if exp_clauses:
    where_exp = " OR ".join(f"({c})" for c in exp_clauses)
    cursor.execute(f"SELECT PDBId FROM exp_rna_rna WHERE {where_exp}")
    exp_ids_to_delete.update(row[0] for row in cursor.fetchall())
    cursor.execute(f"DELETE FROM exp_rna_rna WHERE {where_exp}")

# Optionally, delete from pred_protein_rna and exp_protein_rna by exp_ids_to_delete
if exp_ids_to_delete:
    cursor.executemany("DELETE FROM pred_rna_rna WHERE exp_db_id = ?", [(i,) for i in exp_ids_to_delete])
    cursor.executemany("DELETE FROM exp_rna_rna WHERE PDBId = ?", [(i,) for i in exp_ids_to_delete])

conn.commit()
conn.close()
print(f"Deleted {len(exp_ids_to_delete)} samples based on the selected criteria.")


"export talbe"

# import pandas as pd
#
# # SQL JOIN query to extract relevant columns
# query = """
# SELECT
#     pred.AAmatch_score,
#     pred.AAMotif AS Pred_AAMotif,
#     pred.exp_db_id,
#     exp.AAMotif AS Exp_AAMotif,
#     exp.PDBId
# FROM pred_protein_rna AS pred
# JOIN exp_protein_rna AS exp
# ON pred.exp_db_id = exp.PDBId
# """
#
# # Load into a pandas DataFrame
# df = pd.read_sql_query(query, conn)
#
# # Export to CSV
# df.to_csv("/Users/Iris/Desktop/BachelorProject/AF3InterfaceEval/database/combined_motif_table.csv", index=False)
#
# # Clean up
# conn.close()

# from Bio.Align import PairwiseAligner
# import re
#
# # Input: 3-letter AAMotif with incorrect indices
# aamotif_3letter = "ASP(503), PRO(504), ASP(512), TYR(614), LEU(656), CYS(657), SER(658), ILE(687), GLU(688), MET(691), LEU(692), ARG(713), HIS(715), THR(770), LEU(773), ARG(774), GLN(775), ALA(776), HIS(791), GLY(793), LEU(794), TYR(799), HIS(801), THR(803), SER(804), ARG(807), ARG(808), GLN(852)"
# protein_sequence = "ALFSPHLAESALDLGVQNGTYLRGKLRVSETNCFFGEIRGQWKGHNFERVLLPGRTNLNRAIHGDIVTVELLPVASWRPLREESEGAALARGYTPVGRVVGITTMNRRPFCGSIDVEELNKLALTGTVSVLFQPKDNRIPRIRITTAHLGDLKDKRLSVIIDDWGEHSSFPVGHYVEVLGTIGDKDTEAKVILLENDIPHYDFSEAVYDCLPKGEWNVTEEELGNRLDLRDLCVVSVDPLGCRDIDDALHCRRVNGNHLEVGVHIADVTHFLKEGTAMDEEAAKRSTSVYLVDRRINMLPQLLTENLCSIVADEDRYAFSIMWEFDENYSVVREFFGKTVIRSRAALYYGDAQRMIDDPEDESEAAVSLRYLMQLSRHFRKRREKDGALFLCSQEFKFKVDNDHVNPTDMQAYQTFDSNSMIEEWMLFANAAAARRVYASFPRWTLLRRHQAPAENAFDTLNEAIRRKIGVKLDDTTSLALNESLEKCVDPSDPYFNRLIRTLVTRCLRQAQYFSSSEVSKDEFHHFGLAMPIYTHFTSPIRRYADVIVHRQLAAALGIMDVSEHMVSVKMEALASNLNYRHEQAQKAGRDSQNLFTGFYLRNFANQEIPSEDGYVVKLSETHVFVLVPKYGQEGKIAKETLVRVPNLLDKVKVGIEVRASLVFSIIGLMKG"
#
# # Convert 3-letter motif to 1-letter string
# aa3_to_1 = {
#     'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
#     'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
#     'TYR': 'Y', 'VAL': 'V', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'NEP': 'H', 'HIP': 'H', 'ALY': 'K', 'MLY': 'K',
#     'M3L': 'K', 'MLZ': 'K', '2MR': 'R', 'AGM': 'R', 'MCS': 'C', 'HYP': 'P', 'HY3': 'H', 'LYZ': 'K', 'AHB': 'A',
#     'P1L': 'P', 'SNN': 'S', 'SNC': 'C', 'TRF': 'W', 'KCR': 'K', 'CIR': 'R', 'YHA': 'Y'
# }
#
# # Parse residues and indices
# matches = re.findall(r'([A-Z]{3})\((\d+)\)', aamotif_3letter)
# residues = [aa3_to_1[r] for r, i in matches]
# indices = [int(i) for r, i in matches]
#
# # Insert gaps to preserve spacing
# gapped_motif = residues[0]
# for i in range(1, len(residues)):
#     gap = indices[i] - indices[i - 1] - 1
#     gapped_motif += '-' * gap + residues[i]
#
# # Align using Biopython's PairwiseAligner
# aligner = PairwiseAligner()
# aligner.mode = 'local'
# aligner.open_gap_score = -10
# aligner.extend_gap_score = -0.5
# alignment = aligner.align(protein_sequence, gapped_motif)[0]
# print(alignment)
#
# # Find aligned positions
# prot_start = alignment.aligned[0][0][0]
# motif_pos = 0
# corrected = []
#
# for i, aa in enumerate(gapped_motif):
#     if aa != '-':
#         corrected.append((matches[motif_pos][0], prot_start + i + 1))  # +1 for 1-based
#         motif_pos += 1
#
# # Print corrected AAMotif
# corrected_motif = ', '.join(f"{res}({idx})" for res, idx in corrected)
# print("✅ Corrected AAMotif:")
# print(corrected_motif)


"delete experiment x-ray"

# # Step 1: Get all exp_db_id values where Experiment = "X-RAY DIFFRACTION"
# cursor.execute("SELECT PDBId FROM exp_protein_rna WHERE Experiment = 'X-RAY DIFFRACTION'")
# exp_ids_to_delete = [row[0] for row in cursor.fetchall()]
#
# # Step 2: Delete from pred_protein_rna using those exp_db_id values
# cursor.executemany(
#     "DELETE FROM pred_protein_rna WHERE exp_db_id = ?",
#     [(i,) for i in exp_ids_to_delete]
# )
#
# # Step 3: Delete from exp_protein_rna
# cursor.executemany(
#     "DELETE FROM exp_protein_rna WHERE PDBId = ?",
#     [(i,) for i in exp_ids_to_delete]
# )
#
# # Commit and close
# conn.commit()
# conn.close()



