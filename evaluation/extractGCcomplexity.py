import ast
import json
import re
from collections import Counter
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import RNA as ViennaRNA  # conda install -c bioconda viennarna
from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig
from structures.rna import RNA
from structures.dna import DNA

config = StartConfig()

# Define RNA modifications and their base nucleotides
RNA_MODIFICATIONS = {
    '6MZ': 'A', 'PSU': 'U', '5MC': 'C', 'OMC': 'C', '4OC': 'C', '5MU': 'U', 'OMU': 'U',
    'UR3': 'U', 'A2M': 'A', 'MA6': 'A', '2MG': 'G', 'OMG': 'G', '7MG': 'G', 'RSQ': 'G',
    '5CM': 'C', 'C34': 'C', '5HC': 'C', '6OG': 'G', '5FC': 'C', '3DR': 'A', '2MR': 'A',
    'AGM': 'A'
}

def add_columns():
    """Add sequence feature columns if they don't exist"""
    try:
        db = DatabaseMethods()
        pred_table = config.get_predicted_table_name()
        
        # List of columns to add with their types
        rna_columns = {
            "RNA_GC_Content": "FLOAT",
            "RNA_SequenceComplexity": "FLOAT",
            # Michael Uhl's metrics:
            "RNA_Shannon_Entropy_K1": "FLOAT",
            "RNA_Shannon_Entropy_K2": "FLOAT",
            "RNA_GC_Skew": "FLOAT",
            "RNA_AU_Skew": "FLOAT",
            "RNA_AU_Content": "FLOAT",
            "RNA_MFE_Value": "FLOAT",
            "RNA_MFE_Structure": "TEXT",
            "RNA_Ensemble_Diversity": "FLOAT",
            "RNA_Centroid_Energy": "FLOAT"
        }
        
        # dna_columns = {
        #     "DNA_GC_Content": "FLOAT",
        #     "DNA_SequenceComplexity": "FLOAT",
        #     "DNA_Shannon_Entropy_K1": "FLOAT",
        #     "DNA_Shannon_Entropy_K2": "FLOAT",
        #     "DNA_GC_Skew": "FLOAT",
        #     "DNA_AT_Skew": "FLOAT",
        #     "DNA_AT_Content": "FLOAT"
        # }
        
        # Check if columns exist and add them if they don't
        columns = db.execute_query(f"PRAGMA table_info({pred_table})")
        column_names = [col[1] for col in columns]
        
        # Add RNA columns
        for col_name, col_type in rna_columns.items():
            if col_name not in column_names:
                db.execute_query(f"ALTER TABLE {pred_table} ADD COLUMN {col_name} {col_type}")
        
        # # Add DNA columns
        # for col_name, col_type in dna_columns.items():
        #     if col_name not in column_names:
        #         db.execute_query(f"ALTER TABLE pred_protein_rna ADD COLUMN {col_name} {col_type}")
            
    except Exception as e:
        print(f"ERROR: Failed to add columns: {str(e)}")

def get_longest_sequence(sequence_text):
    """Extract the longest sequence from a string that might be a list"""
    try:
        # Try to parse as JSON/list
        if sequence_text.startswith('['):
            sequences = json.loads(sequence_text)
            if isinstance(sequences, list):
                return max(sequences, key=len)
        return sequence_text
    except:
        # If parsing fails, return the original text
        return sequence_text

def clean_rna_sequence(sequence):
    """
    Clean RNA sequence by converting modifications to their base nucleotides
    and removing any remaining non-standard nucleotides.
    Returns only standard nucleotides (A, U, G, C).
    """
    # First, convert known modifications to their base nucleotides
    for mod, base in RNA_MODIFICATIONS.items():
        sequence = sequence.replace(f"({mod})", base)
    
    # Remove any remaining parentheses and their contents
    sequence = re.sub(r'\([^)]*\)', '', sequence)
    
    # Keep only standard nucleotides
    cleaned = ''.join(c for c in sequence if c in 'AUGC')
    
    return cleaned

def calculate_gc_content(sequence):
    """Calculate GC content as a ratio between 0 and 1"""
    if not sequence:
        return None
    
    # Clean sequence first
    cleaned_seq = clean_rna_sequence(sequence)
    if not cleaned_seq:
        return None
    
    gc_count = cleaned_seq.count('G') + cleaned_seq.count('C')
    total_length = len(cleaned_seq)
    return round(gc_count / total_length, 2) if total_length > 0 else 0

def calculate_linguistic_complexity(sequence):
    """
    Calculate linguistic sequence complexity using k-mer based approach.
    How diverse a sequence is compared to its possible unique substrings.
    LC = Number of observed k-mers/Number of possible k-mers

    Reference:
    Trifonov EN. Making sense of the human genome. In: Sarma RH, Sarma MH,
    editors. Structure and Methods. Volume 1: Human Genome Initiative and DNA
    Recombination. Adenine Press; 1990. p. 69-77.

    This method considers the observed vs possible k-mers for different window sizes.
    Returns a value between 0 and 1.
    """
    # Clean sequence first
    sequence = clean_rna_sequence(sequence)
    
    if not sequence or len(sequence) < 2:
        return 0.0
    
    seq_length = len(sequence)
    k_complexities = []
    
    # Calculate for different k-mer sizes up to sequence length
    for k in range(1, min(seq_length + 1, 6)):  # Limit to 5-mers for efficiency
        # Count observed k-mers
        observed = len(set(sequence[i:i+k] for i in range(seq_length - k + 1)))
        
        # For each k, maximum possible unique k-mers is min(4^k, len-k+1)
        max_possible = min(4**k, seq_length - k + 1)
        
        # Calculate complexity for this k-mer size
        k_complexity = observed / max_possible
        k_complexities.append(k_complexity)
    
    # Calculate average complexity across all k values
    complexity = sum(k_complexities) / len(k_complexities) if k_complexities else 0.0

    # print(f"Linguistic complexity - k-mer complexities: {[f'{x:.2f}' for x in k_complexities]}, final: {complexity:.2f}")
    return round(complexity, 2)

def calculate_shannon_entropy(sequence):
    """
    Calculate Shannon entropy of the sequence.
    Higher values indicate more complex/random sequences.

    Reference:
    Shannon CE. A Mathematical Theory of Communication.
    Bell System Technical Journal. 1948;27(3):379-423.
    DOI: 10.1002/j.1538-7305.1948.tb01338.x

    H = -sum(p_i * log2(p_i))
    where p_i is the probability of each nucleotide
    Returns a value between 0 and 1.
    """
    # Clean sequence first
    sequence = clean_rna_sequence(sequence)
    
    if not sequence:
        return 0.0
    
    # Count frequency of each nucleotide
    counts = Counter(sequence)
    length = len(sequence)
    
    # Calculate entropy
    entropy = 0
    for count in counts.values():
        probability = count / length
        entropy -= probability * math.log2(probability)
    
    # For RNA (4 nucleotides), maximum entropy is 2 bits (log2(4))
    # Normalize to [0,1] by dividing by 2
    normalized_entropy = entropy / 2.0

    # print(f"Shannon entropy - raw: {entropy:.2f}, normalized: {normalized_entropy:.2f}")
    return round(normalized_entropy, 2)

def calculate_sequence_complexity(sequence):
    """
    Calculate overall sequence complexity using multiple metrics.
    Returns a weighted average of different complexity measures.
    The return value is between 0 and 1.
    """
    if not sequence:
        return 0.0
    
    # Show original and cleaned sequence
    cleaned_seq = clean_rna_sequence(sequence)
    # print(f"\nOriginal sequence: {sequence[:50]}...")
    # print(f"Cleaned sequence: {cleaned_seq[:50]}...")
    
    # Calculate different complexity metrics
    linguistic = calculate_linguistic_complexity(cleaned_seq)
    entropy = calculate_shannon_entropy(cleaned_seq)
    
    # Debug prints with 2 decimal places
    # print(f"Final metrics - Linguistic: {linguistic:.2f}, Entropy: {entropy:.2f}")
    
    # Combine metrics (weighted average)
    complexity = round(0.6 * linguistic + 0.4 * entropy, 2)
    # print(f"Combined complexity: {complexity:.2f}\n")
    
    return complexity

def get_ntc_dic(seq, rna=True):
    """Get single nucleotide count dictionary"""
    ntc_dic = {'A': 0, 'C': 0, 'G': 0, 'U' if rna else 'T': 0}
    for nt in seq:
        if nt in ntc_dic:
            ntc_dic[nt] += 1
    return ntc_dic

def get_kmer_counts_dic(seq, k=2, rna=True):
    """Get k-mer count dictionary"""
    bases = ['A', 'C', 'G', 'U' if rna else 'T']
    kmer_dic = {}
    
    # Initialize dictionary with all possible k-mers
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if all(base in bases for base in kmer):
            kmer_dic[kmer] = kmer_dic.get(kmer, 0) + 1
            
    return kmer_dic

def calc_seq_entropy(seq_l, counts_dic, k=1):
    """Calculate Shannon entropy for k-mers"""
    entropy = 0
    total_kmers = sum(counts_dic.values())
    
    if total_kmers == 0:
        return 0
        
    for count in counts_dic.values():
        if count > 0:
            prob = count / total_kmers
            entropy -= prob * math.log2(prob)
            
    return round(entropy, 2)

def calc_seq_gc_skew(seq):
    """Calculate GC skew: (G-C)/(G+C)"""
    ntc_dic = get_ntc_dic(seq)
    g_count = ntc_dic['G']
    c_count = ntc_dic['C']
    if (g_count + c_count) == 0:
        return 0
    return round((g_count - c_count) / (g_count + c_count), 2)

def calc_seq_at_skew(seq, base):
    """Calculate AT/AU skew: (A-T)/(A+T) or (A-U)/(A+U) for RNA"""
    ntc_dic = get_ntc_dic(seq, rna=(base=='U'))
    a_count = ntc_dic['A']
    t_count = ntc_dic[base]
    if (a_count + t_count) == 0:
        return 0
    return round((a_count - t_count) / (a_count + t_count), 2)

def calc_seq_at_content(seq, base):
    """Calculate AT/AU content"""
    ntc_dic = get_ntc_dic(seq, rna=(base=='U'))
    at_count = ntc_dic['A'] + ntc_dic[base]
    total = sum(ntc_dic.values())
    if total == 0:
        return 0
    return round(at_count / total, 2)

def calculate_rna_structure_features(sequence):
    """Calculate RNA structure features using ViennaRNA"""
    try:
        # Create an RNA fold compound
        fc = ViennaRNA.fold_compound(sequence)
        
        # Compute MFE structure and energy
        mfe_structure, mfe_value = fc.mfe()
        
        # Calculate ensemble features
        fc.pf()  # Populate base pair probability matrix
        ensemble_diversity = fc.mean_bp_distance()
        centroid_structure, centroid_energy = fc.centroid()
        
        return {
            'mfe_value': round(mfe_value, 2),
            'mfe_structure': mfe_structure,
            'ensemble_diversity': round(ensemble_diversity, 2),
            'centroid_energy': round(centroid_energy, 2)
        }
    except Exception as e:
        print(f"Error calculating RNA structure features: {str(e)}")
        return {
            'mfe_value': None,
            'mfe_structure': None,
            'ensemble_diversity': None,
            'centroid_energy': None
        }

def process_sequences():
    """Process all sequences in the database"""
    try:
        db = DatabaseMethods()
        pred_table = config.get_predicted_table_name()
        
        # Ensure columns exist
        add_columns()
        
        # Get all sequences
        query = f"SELECT exp_db_id, RNASequence FROM {pred_table} WHERE RNASequence IS NOT NULL"
        sequences = db.execute_query(query)
        
        for pdb_id, seq_text in sequences:
            print(f"\nProcessing {pdb_id}...")
            
            # Check if this is a DNA or RNA sequence
            try:
                dna = DNA.get_dna_from_db(id=pdb_id, file_name=None)
                rna = RNA.get_rna_from_db(id=pdb_id, file_name=None)
                
                is_dna = dna.is_dna if dna else False
                is_rna = rna.is_rna if rna else False
                
                if not (is_dna or is_rna):
                    print(f"Skipping {pdb_id}: Neither DNA nor RNA")
                    continue
                    
            except Exception as e:
                print(f"Error checking DNA/RNA status for {pdb_id}: {str(e)}")
                continue
            
            # Get the longest sequence if it's a list
            sequence = get_longest_sequence(seq_text)
            if not sequence:
                continue
                
            # Clean sequence
            cleaned_seq = clean_rna_sequence(sequence)
            if not cleaned_seq:
                continue
            
            # Calculate common features
            gc_content = calculate_gc_content(cleaned_seq)
            seq_complexity = calculate_sequence_complexity(cleaned_seq)
            
            # Calculate additional sequence features
            mono_nt_counts = get_ntc_dic(cleaned_seq, rna=is_rna)
            di_nt_counts = get_kmer_counts_dic(cleaned_seq, k=2, rna=is_rna)
            
            shannon_k1 = calc_seq_entropy(len(cleaned_seq), mono_nt_counts)
            shannon_k2 = calc_seq_entropy(len(cleaned_seq), di_nt_counts, k=2)
            gc_skew = calc_seq_gc_skew(cleaned_seq)
            
            # Calculate RNA-specific features
            if is_rna:
                struct_features = calculate_rna_structure_features(cleaned_seq)
                au_skew = calc_seq_at_skew(cleaned_seq, 'U')
                au_content = calc_seq_at_content(cleaned_seq, 'U')
                
                # Update RNA features
                db.execute_query(f"""
                    UPDATE {pred_table} 
                    SET RNA_GC_Content = ?,
                        RNA_SequenceComplexity = ?,
                        RNA_Shannon_Entropy_K1 = ?,
                        RNA_Shannon_Entropy_K2 = ?,
                        RNA_GC_Skew = ?,
                        RNA_AU_Skew = ?,
                        RNA_AU_Content = ?,
                        RNA_MFE_Value = ?,
                        RNA_MFE_Structure = ?,
                        RNA_Ensemble_Diversity = ?,
                        RNA_Centroid_Energy = ?
                    WHERE exp_db_id = ?
                """, (gc_content, seq_complexity, shannon_k1, shannon_k2, 
                      gc_skew, au_skew, au_content,
                      struct_features['mfe_value'], struct_features['mfe_structure'],
                      struct_features['ensemble_diversity'], struct_features['centroid_energy'],
                      pdb_id))
            
            # Calculate DNA-specific features
            if is_dna:
                at_skew = calc_seq_at_skew(cleaned_seq, 'T')
                at_content = calc_seq_at_content(cleaned_seq, 'T')
                
                # Update DNA features
                db.execute_query(f"""
                    UPDATE {pred_table} 
                    SET DNA_GC_Content = ?,
                        DNA_SequenceComplexity = ?,
                        DNA_Shannon_Entropy_K1 = ?,
                        DNA_Shannon_Entropy_K2 = ?,
                        DNA_GC_Skew = ?,
                        DNA_AT_Skew = ?,
                        DNA_AT_Content = ?
                    WHERE exp_db_id = ?
                """, (gc_content, seq_complexity, shannon_k1, shannon_k2,
                      gc_skew, at_skew, at_content,
                      pdb_id))
            
            print(f"Updated features for {pdb_id} - {'RNA' if is_rna else 'DNA'}")
            
    except Exception as e:
        print(f"ERROR: Failed to process sequences: {str(e)}")
    finally:
        db.close_connection()

if __name__ == "__main__":
    process_sequences()
