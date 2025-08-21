import json
import os
import torch
import sys
import collections
import subprocess
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import forgi.graph.bulge_graph as fgb

from typing import Dict, List, Tuple, Set, Optional, Union
from pathlib import Path

from database.databaseMethods import DatabaseMethods
from database.startConfig import StartConfig

config = StartConfig()

def db2pairs(structure, start_index=0):
    """
    Converts dot-bracket string into a list of pairs.

    Input:
      structure <string>: A sequence in dot-bracket format.
      start_index <int>: Starting index of first nucleotide (default zero-indexing).

    Returns:
      pairs <list>: A list of tuples of (index1, index2, pk_level).

    """
    level_stacks = collections.defaultdict(list)
    closing_partners = {')': '(', ']': '[', '}': '{', '>': '<'}
    levels = {')': 0, ']': 1, '}': 2, '>': 3}

    pairs = []

    for i, sym in enumerate(structure, start_index):
        if sym == '.':
            continue
        # high order pks are alphabetical characters
        if sym.isalpha():
            if sym.isupper():
                level_stacks[sym].append(i)
            else:
                try:  # in case we have invalid preditions, we continue with next bracket
                    op = level_stacks[sym.upper()].pop()
                    pairs.append((op, i,
                                  ord(sym.upper()) - 61))  # use asci code if letter is used to asign PKs, start with level 4 (A has asci code 65)
                except:
                    continue
        else:
            if sym in closing_partners.values():
                level_stacks[sym].append(i)
            else:
                try:  # in case we have invalid preditions, we continue with next bracket
                    op = level_stacks[closing_partners[sym]].pop()
                    pairs.append([op, i])  # , levels[sym]])
                except:
                    continue
    return sorted(pairs, key=lambda x: x[0])


def to_bpseq(pairs, sequence, outpath):
    bpseq_output = [[str(i), s, '0']
                    for i, s in enumerate(sequence, 1)]
    if pairs:
        pos1, pos2 = zip(*pairs)
        for p1, p2 in zip(pos1, pos2):
            bpseq_output[p1][2] = str(p2+1)
            bpseq_output[p2][2] = str(p1+1)
    with open(outpath, 'w+') as f:
        for line in bpseq_output:
            f.write('\t'.join(line)+'\n')

def annotate_pairs(pairs, sequence):
    watson_pairs, wobble_pairs, other_pairs = type_pairs(pairs, sequence)
    lone_pairs = lone_pair(pairs)
    multiplets = multiplets_pairs(pairs)
    return {
        'sequence': sequence,
        'pairs': pairs,
        'watson_pairs': watson_pairs,
        'wobble_pairs': wobble_pairs,
        'nc_pairs': other_pairs,
        'lone_pairs': lone_pairs,
        'multiplets': multiplets
    }

class BpRNA():
    def __init__(self, bprna_dir, bpseq_path):
        self.cwd = str(Path(bprna_dir).resolve())
        # self.working_dir = str(Path(working_dir).resolve())
        # self.bpseq_path = Path(self.working_dir, 'tmp.bpseq')
        self.bpseq_path = Path(bpseq_path)
        self.st_file = Path(bprna_dir, f"{self.bpseq_path.stem}.st")

    def run(self):
        p = subprocess.call(["perl", "bpRNA.pl", Path(self.bpseq_path.resolve())], cwd=self.cwd)  # , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Path(self.bpseq_path.resolve()).unlink()
        return p
    
    def parse_st_output(self):
        try:
            with open(self.st_file) as f:
                lines = f.readlines()

        except pd.errors.ParserError:
            print("# Read st ERROR: ParserError st file", self.st_file)
            return False

        lines_plain = []
        header_lines = []
        for line in lines:
            if not line.startswith('#'):
                lines_plain.append(line.strip())
            else:
                header_lines.append(line.strip())
        # lines = lines_plain
        # print(lines_plain)
        # print(header_lines)

        annotations = {
            'name': header_lines[0].split(': ')[-1],
            'length': int(''.join(c for c in header_lines[1] if c.isdigit()).split(': ')[-1]),
            'page_number': int(header_lines[2].split(': ')[-1]),
            'sequence': ''.join([a for a in lines_plain[0] if a != '\n']),
            'secondary_structure': ''.join([a for a in lines_plain[1] if a != '\n']),
            'annotation_string': ''.join([a for a in lines_plain[2] if a != '\n']),
            'pseudoknot_annotation_string': ''.join([a for a in lines_plain[3] if a != '\n']),
            'additional_annotations': [a for a in lines_plain[4:] if a != '\n']
        }

        Path(self.st_file.resolve()).unlink()

        return annotations

################################################################################
# following code snippets originate from SPOT-RNA repo at https://github.com/jaswindersingh2/SPOT-RNA
# we use them for plotting with varna in the visualization module
################################################################################


# copy-paste from SPOT-RNA2 source code
def lone_pair(pairs):
    lone_pairs = []
    pairs.sort()
    for i, I in enumerate(pairs):
        if ([I[0] - 1, I[1] + 1] not in pairs) and ([I[0] + 1, I[1] - 1] not in pairs):
            lone_pairs.append(I)

    return lone_pairs


# copy-paste from SPOT-RNA2 source code
def type_pairs(pairs, sequence):
    sequence = [i.upper() for i in sequence]
    # seq_pairs = [[sequence[i[0]],sequence[i[1]]] for i in pairs]

    AU_pair = []
    GC_pair = []
    GU_pair = []
    other_pairs = []
    for i in pairs:
        if [sequence[i[0]],sequence[i[1]]] in [["A","U"], ["U","A"]]:
            AU_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","C"], ["C","G"]]:
            GC_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","U"], ["U","G"]]:
            GU_pair.append(i)
        else:
            other_pairs.append(i)
    watson_pairs_t = AU_pair + GC_pair
    wobble_pairs_t = GU_pair
    other_pairs_t = other_pairs
    # print(watson_pairs_t)
    return watson_pairs_t, wobble_pairs_t, other_pairs_t



# copy-paste from SPOT-RNA2 source code
def multiplets_pairs(pred_pairs):

    pred_pair = [i[:2] for i in pred_pairs]
    temp_list = flatten(pred_pair)
    temp_list.sort()
    new_list = sorted(set(temp_list))
    dup_list = []
    for i in range(len(new_list)):
        if (temp_list.count(new_list[i]) > 1):
            dup_list.append(new_list[i])

    dub_pairs = []
    for e in pred_pair:
        if e[0] in dup_list:
            dub_pairs.append(e)
        elif e[1] in dup_list:
            dub_pairs.append(e)

    temp3 = []
    for i in dup_list:
        temp4 = []
        for k in dub_pairs:
            if i in k:
                temp4.append(k)
        temp3.append(temp4)

    return temp3

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result



class SecondaryStructureAnalyser:
    def __init__(self, dssr_file: str):
        self.dssr_file = dssr_file
        self.dssr_data = self._load_dssr_file(self.dssr_file)
        self.bg = self._create_bulge_graph()
    
    @staticmethod
    def _load_dssr_file(dssr_file: str) -> dict:
        try:
            with open(dssr_file) as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading DSSR file: {e}")
            return {}
            
    def _create_bulge_graph(self) -> fgb.BulgeGraph:
        """Create forgi BulgeGraph from dot-bracket notation"""
        if 'dbn' not in self.dssr_data:
            return None
        
        dbn_data = self.dssr_data['dbn']
        if 'all_chains' not in dbn_data:
            return None
        
        structure = dbn_data['all_chains']['sstr']
        sequence = dbn_data['all_chains']['bseq']
        self.sequence = sequence
        self.structure = structure

        try:
            # classmethod from_dotbracket(dotbracket_str, seq=None, name=None, dissolve_length_one_stems=False, remove_pseudoknots=False)
            bg = fgb.BulgeGraph.from_dotbracket(dotbracket_str=structure, 
                                              seq=sequence,
                                              remove_pseudoknots=False)
            return bg
        except Exception as e:
            # print(f"Error: {e}")
            return None

    # https: // viennarna.github.io / forgi / graph_tutorial.html
    def get_stems(self) -> List[Dict]:
        stems = []
        if self.bg:
            for stem in self.bg.stem_iterator():
                stem_data = {
                    'stem_id': stem,
                    'length': self.bg.stem_length(stem),
                    'pairs': list(self.bg.stem_bp_iterator(stem)),
                    'residues': [
                        list(self.bg.define_residue_num_iterator(stem, 0)),
                        list(self.bg.define_residue_num_iterator(stem, 1))
                    ]
                }
                stems.append(stem_data)
        return stems

    def get_hairpin_loops(self) -> List[Dict]:
        hairpins = []
        if self.bg:
            for element in self.bg.hloop_iterator():
                loop = {
                    'hairpin_id': element,
                    'length': len(list(self.bg.define_range_iterator(element))),
                    'closing_pair': self.bg.connections(element)[0],  # First connected stem
                    'residues': list(self.bg.define_residue_num_iterator(element))
                }
                hairpins.append(loop)
        return hairpins

    def get_internal_loops(self) -> List[Dict]:
        internal_loops = []
        if self.bg:
            for element in self.bg.iloop_iterator():
                internal = {
                    'internal_id': element,
                    'length': len(list(self.bg.define_range_iterator(element))),
                    'closing_pairs': self.bg.connections(element),
                    'residues': [
                        list(self.bg.define_residue_num_iterator(element, 0)),
                        list(self.bg.define_residue_num_iterator(element, 1))
                    ]
                }
                internal_loops.append(internal)
        return internal_loops

    def get_multibranch_loops(self) -> List[Dict]:
        """Get multibranch loop (junction) information"""
        junctions = []
        if self.bg:
            for junction in self.bg.mloop_iterator():
                # Get connected stems
                connected_stems = list(self.bg.connections(junction))
                
                # Get residues if any exist
                try:
                    residues = list(self.bg.define_residue_num_iterator(junction))
                except Exception:
                    residues = []
                
                junction_data = {
                    'junction_id': junction,
                    'degree': len(connected_stems),  # Number of stems meeting at this junction
                    'residues': residues,
                    'connected_stems': connected_stems  # Add this to see which stems are connected
                }
                junctions.append(junction_data)
        return junctions
    
    def get_element_string(self) -> str:
        # Get string representation of all elements in BulgeGraph
        if self.bg:
            return self.bg.to_element_string()
        return ""
            
    # def get_dangling_ends(self) -> Tuple[List[str], List[str]]:
    #     dangling_start = []
    #     dangling_end = []
    #     if self.bg:
    #         for end in self.bg.exterior_iterator():
    #             residues = list(self.bg.define_residue_num_iterator(end))
    #             if not residues:
    #                 continue
    #             if self.bg.define_a(end) == 0:  # 5' end
    #                 dangling_start.extend(residues)
    #             else:  # 3' end
    #                 dangling_end.extend(residues)
    #     return dangling_start, dangling_end
    
    
    def extract_dangling_ends(self):
        """
        Extract dangling ends (5', 3', and internal) using forgi functions.
        Ensures only actual unpaired extensions of stems are considered.
        """
        # Dictionary to store dangling end positions with structured output
        dangling_ends = {
            "5' end": [],      # [(start, end)]
            "3' end": [],      # [(start, end)]
            "internal": []     # [{"stem_id": ..., "side": ..., "residues": [(start, end)]}]
        }
    
        if not self.bg:
            return dangling_ends
    
        def get_contiguous_regions(residues):
            """Convert list of residues into contiguous (start, end) tuples."""
            if not residues:
                return []
            residues = sorted(residues)
            regions = []
            start = residues[0]
            prev = start
            for res in residues[1:]:
                if res != prev + 1:  # Gap detected
                    regions.append((start, prev))
                    start = res
                prev = res
            regions.append((start, prev))  # Add last region
            return regions
    
        # 1. Identify 5' and 3' terminal dangling ends
        floop_residues = []
        for floop in self.bg.floop_iterator():
            if self.bg.is_single_stranded(floop):
                floop_residues.extend(self.bg.define_residue_num_iterator(floop))
        dangling_ends["5' end"] = get_contiguous_regions(floop_residues)
    
        tloop_residues = []
        for tloop in self.bg.tloop_iterator():
            if self.bg.is_single_stranded(tloop):
                tloop_residues.extend(self.bg.define_residue_num_iterator(tloop))
        dangling_ends["3' end"] = get_contiguous_regions(tloop_residues)
    
        # 2. Identify internal dangling ends (true extensions of stems)
        # for stem in self.bg.stem_iterator():
        #     start, end = stem  # Stem positions
        #     left_flanking = []
        #     right_flanking = []
    
        #     # Get flanking regions at the 5' and 3' end of each stem
        #     for side in [0, 1]:  # 0 = 5' end of stem, 1 = 3' end
        #         flank_res = self.bg.flanking_nuc_at_stem_side(stem, side)
    
        #         # Collect contiguous unpaired residues extending beyond the stem
        #         temp_flanking = []
        #         while (
        #             flank_res and 
        #             flank_res != 0 and 
        #             flank_res != self.bg.seq_length + 1 and 
        #             self.bg.pairing_partner(flank_res) is None and
        #             self.bg.is_single_stranded(self.bg.get_elem(flank_res))  # Ensure it's truly unpaired
        #         ):
        #             temp_flanking.append(flank_res)
        #             flank_res = flank_res - 1 if side == 0 else flank_res + 1  # Move left or right
    
        #         # Validate: Ensure residues are not inside a hairpin or loop
        #         temp_flanking = [
        #             res for res in temp_flanking if self.bg.get_elem(res) not in self.bg.hloop_iterator()
        #         ]
    
        #         # Store contiguous regions only if they extend the stem
        #         if temp_flanking:
        #             if side == 0:
        #                 left_flanking = temp_flanking
        #             else:
        #                 right_flanking = temp_flanking
    
        #     # Ensure these flanking residues are **not** the same as the terminal dangling ends
        #     left_flanking = [
        #         (min(left_flanking), max(left_flanking))
        #         for region in get_contiguous_regions(left_flanking)
        #         if region not in dangling_ends["5' end"]
        #     ]
    
        #     right_flanking = [
        #         (min(right_flanking), max(right_flanking))
        #         for region in get_contiguous_regions(right_flanking)
        #         if region not in dangling_ends["3' end"]
        #     ]
    
        #     # Store structured information per stem in a dictionary format
        #     if left_flanking or right_flanking:
        #         stem_entry = {
        #             "stem_id": stem,
        #             "dangling_ends": []
        #         }
        #         if left_flanking:
        #             stem_entry["dangling_ends"].append({"side": "5'", "residues": left_flanking})
        #         if right_flanking:
        #             stem_entry["dangling_ends"].append({"side": "3'", "residues": right_flanking})
    
        #         dangling_ends["internal"].append(stem_entry)
    
        return dangling_ends
    
    
    def forms_duplex(self) -> bool:
        # Determine if RNA forms a duplex structure
        if not self.bg:
            return False
        has_stem = any(True for _ in self.bg.stem_iterator())
        has_hairpin = any(True for _ in self.bg.hloop_iterator())
        return has_stem or has_hairpin
    
    def get_all_features(self) -> Dict:
        stems = self.get_stems()
        hairpins = self.get_hairpin_loops()
        internal_loops = self.get_internal_loops()
        multibranch = self.get_multibranch_loops()
        # dangling_start, dangling_end = self.get_dangling_ends()
        dangling_ends = self.extract_dangling_ends()
        pseudoknots = self.bg.pseudoknotted_basepairs() if self.bg else []
        forgi_estring = self.get_element_string()
        
        # Count features
        feature_counts = [
            len(stems),
            sum(1 for loop in hairpins if loop['length'] > 0),
            len(internal_loops),
            sum(1 for loop in multibranch if loop['degree'] > 0),
            # len(dangling_start) + len(dangling_end),
            # len(pseudoknots)
        ]
        
        return {
            'annotation_string': forgi_estring,
            'stems': stems,
            'hairpin_loops': hairpins,
            'internal_loops': internal_loops,
            'multibranch_loops': multibranch,
            'dangling_ends': dangling_ends,
            # 'dangling_ends': (dangling_start, dangling_end),
            'pseudoknots': pseudoknots,
            'forms_duplex': self.forms_duplex(),
            'sequence': self.sequence,
            'structure': self.structure,
            'feature_counts': feature_counts
        }

def analyse_secondary_structure(dssr_file: str) -> Dict:
    # processor = SecondaryStructureProcessor(dssr_file)
    analyser = SecondaryStructureAnalyser(dssr_file)
    return analyser.get_all_features()

def process_dssr_files(dssr_dir: str, bprna_dir: Optional[Union[str, Path]] = None) -> None:
    db = DatabaseMethods()

    # Check and create columns if they don't exist
    pred_table = config.get_predicted_table_name()
    exp_table = config.get_reference_table_name()
    columns_to_check = [
        'RNA_Stems TEXT',
        'RNA_HairpinLoops TEXT',
        'RNA_InternalLoops TEXT',
        'RNA_MultibranchLoops TEXT',
        'RNA_DanglingEnds TEXT',
        'RNA_Pseudoknots TEXT',
        'RNA_isDuplex INTEGER'
    ]

    # Check and create columns for both tables
    for table in [pred_table, exp_table]:
        for column_def in columns_to_check:
            column_name = column_def.split()[0]
            try:
                db.execute_query(f"SELECT {column_name} FROM {table} LIMIT 1")
            except:
                print(f"Adding column {column_name} to {table}")
                db.execute_query(f"ALTER TABLE {table} ADD COLUMN {column_def}")
    results = []

    for filename in os.listdir(dssr_dir):
        if filename.endswith('_dssr.json'):
            dssr_path = os.path.join(dssr_dir, filename)
            features = analyse_secondary_structure(dssr_path)
            pdb_id = filename.split('_')[0].upper()

            # Extract data for database update
            stems_data = [stem['residues'][0] for stem in features.get('stems', [])]
            hairpin_data = [loop['residues'] for loop in features.get('hairpin_loops', [])]
            internal_data = [loop['residues'][0] for loop in features.get('internal_loops', [])]
            multibranch_data = [loop['residues'] for loop in features.get('multibranch_loops', [])]
            dangling_data = [
                features['dangling_ends'].get("5' end", []),
                features['dangling_ends'].get("3' end", [])
            ]
            pseudoknots = features.get('pseudoknots', [])
            is_duplex = 1 if features.get('forms_duplex', 0) == 1 else 0

            # Update database based on file type
            if '_af_' in filename:
                query = f"""
                UPDATE {pred_table}
                SET RNA_Stems = ?,
                    RNA_HairpinLoops = ?,
                    RNA_InternalLoops = ?,
                    RNA_MultibranchLoops = ?,
                    RNA_DanglingEnds = ?,
                    RNA_Pseudoknots = ?,
                    RNA_isDuplex = ?
                WHERE exp_db_id = ?
                """
            elif '_pdb_' in filename:
                query = f"""
                UPDATE {exp_table}
                SET RNA_Stems = ?,
                    RNA_HairpinLoops = ?,
                    RNA_InternalLoops = ?,
                    RNA_MultibranchLoops = ?,
                    RNA_DanglingEnds = ?,
                    RNA_Pseudoknots = ?,
                    RNA_isDuplex = ?
                WHERE PDBId = ?
                """
            else:
                continue

            # Update database
            values = (
                str(stems_data),
                str(hairpin_data),
                str(internal_data),
                str(multibranch_data),
                str(dangling_data),
                str(pseudoknots),
                is_duplex,
                pdb_id
            )
            db.execute_query(query, values)

            # bpRNA annotations
            annotations = {}
            if bprna_dir is not None:
                db = features['structure']
                sequence = features['sequence']
                
                pairs = db2pairs(structure=db)
                pair_annotiations = annotate_pairs(pairs, sequence)
                flat_multi_annotations = [p for mp in pair_annotiations['multiplets'] for p in mp]
                multiplet_free_pairs = [p for p in pairs if p not in flat_multi_annotations]
                outpath = f"{pdb_id}.bpseq"
                to_bpseq(multiplet_free_pairs, sequence, outpath)
                bprna = BpRNA(bprna_dir, outpath)
                p = bprna.run()
                if p == 0:
                    annotations = bprna.parse_st_output()
                    annotations.update(pair_annotiations)
            
            results.append({
                'pdb_id': pdb_id,
                # 'forgi_features': features,
                # 'bprna_features': annotations
                'features': features
            })

    # Write detailed analysis to file
    output_path = os.path.join(config.parent_folder, 'secondary_features.txt')
    with open(output_path, 'w') as f:
        for result in results:
            f.write(f"\nRNA Secondary Structure Analysis for {result['pdb_id']}")
            
            features = result['features']
            f.write(f"\nForms duplex: {features['forms_duplex']}\n")
            
            # Write stem information
            f.write("\nStems:")
            if features['stems']:
                for stem in features['stems']:
                    f.write(f"\n  - ID: {stem['stem_id']}")
                    f.write(f"\n    Length: {stem['length']}")
                    f.write(f"\n    Base pairs: {len(stem['pairs'])}")
                    f.write(f"\n    Residues: {stem['residues']}\n")
            else:
                f.write("\n  No stems found\n")
            
            # Write hairpin loops
            f.write("\nHairpin Loops:")
            if features['hairpin_loops']:
                for loop in features['hairpin_loops']:
                    f.write(f"\n  - ID: {loop['hairpin_id']}")
                    f.write(f"\n    Length: {loop['length']}")
                    f.write(f"\n    Residues: {loop['residues']}\n")
            else:
                f.write("\n  No hairpin loops found\n")
            
            # Write internal loops
            f.write("\nInternal Loops:")
            if features['internal_loops']:
                for loop in features['internal_loops']:
                    f.write(f"\n  - ID: {loop['internal_id']}")
                    f.write(f"\n    Length: {loop['length']}")
                    f.write(f"\n    Residues: {loop['residues']}\n")
            else:
                f.write("\n  No internal loops found\n")
            
            # Write multibranch loops
            f.write("\nMultibranch Loops:")
            if features['multibranch_loops']:
                for loop in features['multibranch_loops']:
                    f.write(f"\n  - ID: {loop['junction_id']}")
                    f.write(f"\n    Degree: {loop['degree']}")
                    f.write(f"\n    Connected stems: {loop['connected_stems']}")
                    f.write(f"\n    Unpaired residues: {loop['residues']}")
                    if not loop['residues']:
                        f.write(" (Direct stem junction)")
                    f.write("\n")
            else:
                f.write("\n  No multibranch loops found\n")
            
            # Write dangling ends
            f.write("\nDangling Ends:")
            dangling_ends = features['dangling_ends']
            f.write("\n  5' end: {}".format(dangling_ends.get("5' end", [])))
            f.write("\n  3' end: {}\n".format(dangling_ends.get("3' end", [])))
            
            # Write pseudoknots
            f.write("\nPseudoknots:")
            if features['pseudoknots']:
                f.write(f"\n  {features['pseudoknots']}\n")
            else:
                f.write("\n  No pseudoknots found\n")
            
            f.write("-"*50 + "\n")

    return results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <dssr_directory_path>")
        sys.exit(1)
    dssr_dir = sys.argv[1]
    # bprna_path = '/home/fred/research'
    bprna_dir = '/home/fred/current_projects/github/RnaBench/external_algorithms/bpRNA'  # optional
    # results = process_dssr_files(dssr_dir, bprna_dir=bprna_dir)  # bprna-dir is optional
    results = process_dssr_files(dssr_dir, bprna_dir=None)  # bprna-dir is optional
    # print(results)