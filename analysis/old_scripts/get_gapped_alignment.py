import json

from pathlib import Path
from matplotlib import pyplot as plt

def fetch_msa_string(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    for key, value in data.items():
        if key == "sequences":
            for chain in value:
                if "rna" in chain:
                    
                    msa = chain["rna"]["unpairedMsa"]
                    break
    return msa

def msa_string_to_list(msa_string):
    return msa_string.split('\n')

def get_sequences_from_msa(msa_string):
    msa_list = msa_string_to_list(msa_string)
    labels = []
    sequences = []
    for line in msa_list:
        if not line.startswith('>') and line.strip():
            sequences.append(line.strip())
        else:
            labels.append(line.strip())
    return labels, sequences

def pad_insertions(seqs):
    """
    Given a list of aligned RNA sequences where:
      - Uppercase letters and '-' are the consensus columns.
      - Lowercase letters are insertions relative to consensus.
    Returns a new list of sequences all the same length, with '-' 
    padding inserted so that every insertion column is represented 
    in every sequence.
    """
    cons_lists = []
    ins_lists = []

    # Step 1: Split each sequence into 'consensus' chars and insertion segments
    for seq in seqs:
        cons = []
        ins = []
        current_ins = ''
        for ch in seq:
            if ch.islower():
                # build up an insertion
                current_ins += ch
            else:
                # commit any insertion we just passed
                ins.append(current_ins)
                current_ins = ''
                # record the consensus column (A/C/G/U or '-')
                cons.append(ch)
        # tail insertion after the last consensus column
        ins.append(current_ins)
        cons_lists.append(cons)
        ins_lists.append(ins)

    # Step 2: sanity-check that every seq has the same number of consensus columns
    n_cons = len(cons_lists[0])
    if any(len(c) != n_cons for c in cons_lists):
        raise ValueError("Consensus column count differs between sequences")

    # Step 3: for each insertion slot (there are n_cons+1 of them), find the max insertion length
    max_ins = []
    for i in range(n_cons + 1):
        max_ins.append(max(len(ins_lists[s][i]) for s in range(len(seqs))))

    # Step 4: rebuild each sequence, padding its insertions with '-' to the max length
    padded = []
    for cons, ins in zip(cons_lists, ins_lists):
        out = []
        for i in range(n_cons):
            # pad this insertion segment, then add the consensus char
            out.append(ins[i].ljust(max_ins[i], '-'))
            out.append(cons[i])
        # final insertion after the last consensus char
        out.append(ins[n_cons].ljust(max_ins[n_cons], '-'))
        padded.append(''.join(out))

    return padded

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Process RNA MSA data.")
    parser.add_argument('--algorithm', type=str, help='Algorithm that produced the datafiles')
    parser.add_argument('--data_dir', type=str, help='Directory of the original json datafiles')
    parser.add_argument('--out_dir', type=str, help='Output directory')

    args = parser.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    json_files = list(Path(args.data_dir).glob("*.json"))
    
    for p in json_files:
        pdb_id = str(p.stem)[:4].upper()
        algorithm = args.algorithm
        out_dir = args.out_dir
        
        print(f"Processing {p} for PDB ID {pdb_id} with algorithm {algorithm}")
        
        # Construct the path to the JSON file
        json_path = Path(p)
        
        # Fetch MSA string from JSON file
        msa_string = fetch_msa_string(json_path)
        
        if msa_string is None:
            print(f"No MSA found in {p}")
            continue
        # elif len(msa_string_to_list(msa_string)) -1 <= 2:
        #     print(f"Skipping {p} due to insufficient sequences.")
        #     continue
        else:
            # Convert MSA string to list of sequences
            labels, sequences = get_sequences_from_msa(msa_string)
            
            # Pad insertions in sequences
            padded_sequences = pad_insertions(sequences)
            
            # Print results
            print('MSA from', algorithm, 'for', pdb_id, 'contains', len(sequences), 'sequences')
            
            Path(out_dir).mkdir(parents=True, exist_ok=True)
            
            # Write padded sequences to FASTA file
            print(f"Writing gapped RNA alignment to {out_dir}/{pdb_id}_gapped_rna_alignment_{algorithm}.fasta")
            with open(Path(out_dir, f"{pdb_id}_gapped_rna_alignment_{algorithm}.fasta"), 'w') as f:
                for label, sequence in zip(labels, padded_sequences):
                    f.write(f"{label}\n")
                    f.write(''.join(sequence.upper()) + '\n')
            
            print(f"Labels: {labels[:10]}")  # Print first 10 labels for brevity
            for s in padded_sequences[:10]:  # Print first 10 sequences for brevity
                print(''.join(s.upper()))
            
            print('Done writing gapped RNA alignment.')
    
