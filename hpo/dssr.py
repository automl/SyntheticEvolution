"""Run DSSR in isolation; map its full pairs list through validated nucleotide IDs.

Supports the initial contract: one unmodified RNA, chain A, AF3 residue numbering
1..N. Reject missing/reordered/mismatched residues rather than silently shifting
positions. Verified against the supplied DSSR 2.9.3 output.
"""
import json
import subprocess
from pathlib import Path
from scoring import normalize_pairs


def parse_pairs(data, sequence):
    nts = data.get('nts', [])
    if len(nts) != len(sequence):
        raise ValueError('DSSR nucleotide count does not match input')
    mapping = {}
    for nt in nts:
        pos = nt['nt_resnum'] - 1
        if (nt['chain_name'] != 'A' or nt.get('nt_type') != 'RNA'
                or not 0 <= pos < len(sequence)
                or nt['nt_code'] != sequence[pos]
                or nt['nt_id'] in mapping):
            raise ValueError('Unexpected DSSR chain, sequence or residue identifiers')
        mapping[nt['nt_id']] = pos
    if set(mapping.values()) != set(range(len(sequence))):
        raise ValueError('Missing or duplicate DSSR residue positions')
    pairs = data.get('pairs', [])
    if data.get('num_pairs', len(pairs)) != len(pairs):
        raise ValueError('Inconsistent DSSR pair count')
    return normalize_pairs([(mapping[p['nt1']], mapping[p['nt2']]) for p in pairs], len(sequence))


def run_dssr(model, directory, sequence, executable='x3dna-dssr', timeout=300):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / 'dssr.json'
    with (directory / 'dssr.log').open('w') as log:
        subprocess.run([executable, f'-i={Path(model).resolve()}', '--json', f'-o={output}'],
                       cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=timeout)
    return parse_pairs(json.loads(output.read_text()), sequence)
