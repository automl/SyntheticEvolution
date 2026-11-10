"""Extraction step: turn a predicted structure into a zero-based pair list.

* RNA: run DSSR in isolation and map its full pairs list through validated
  nucleotide IDs (verified against DSSR 2.9.3 output).
* Protein: read C-alpha positions with Biopython, build a distance matrix and
  binary-mask it with a cutoff to obtain residue contacts.

Initial contract: one unmodified chain A with AF3 residue numbering 1..N.
Missing/reordered/mismatched residues are rejected rather than silently
shifting positions. Every extractor returns the same ``list[(i, j)]`` so scoring
does not depend on the molecule type.
"""
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Callable, Union

# Files an extractor publishes into the archive, one per chain type.
RESULT_FILES = ('dssr.json', 'contacts.json')


def parse_pairs(data: dict[str, Any], sequence: str) -> list[tuple[int, int]]:
    """Validate DSSR nucleotide mappings and return zero-based base pairs.

    DSSR uses one-based residue numbers and nucleotide identifiers in its pair
    records. The mapping is validated before pairs are converted so malformed,
    reordered, or mismatched structures cannot silently produce wrong indices.
    """
    nts = data.get('nts', [])
    if len(nts) != len(sequence):
        raise ValueError('DSSR nucleotide count does not match input')

    # Map DSSR residue identifiers to the corresponding zero-based sequence
    # position while checking the expected chain, residue type, and base code.
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

    # Preserve DSSR's pair order while translating nucleotide IDs to indices.
    pairs = data.get('pairs', [])
    if data.get('num_pairs', len(pairs)) != len(pairs):
        raise ValueError('Inconsistent DSSR pair count')
    return [(mapping[p['nt1']], mapping[p['nt2']]) for p in pairs]


def run_dssr(
    model: Union[str, Path],
    directory: Union[str, Path],
    sequence: str,
    executable: str = 'x3dna-dssr',
    timeout: int = 300,
    scratch: Union[str, Path, None] = None,
) -> list[tuple[int, int]]:
    """Run DSSR for one predicted model and parse its JSON pair output.

    Only ``dssr.json`` is written to ``directory``. ``dssr.log`` and DSSR's
    working directory go to ``scratch`` (default: ``directory``).
    """
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    scratch = directory if scratch is None else Path(scratch).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    model_path = Path(model).resolve()
    log_path = scratch / 'dssr.log'
    result_path = directory / 'dssr.json'

    # DSSR writes its side files (dssr-*.pdb/.ct/.bpseq/...) into its working
    # directory and has no flag to suppress them, so run it in a temporary
    # directory and keep only the JSON output.
    with tempfile.TemporaryDirectory(prefix='dssr-', dir=scratch) as work_directory:
        output_path = Path(work_directory) / 'dssr.json'
        dssr_command = [
            executable,
            f'-i={model_path}',
            '--json',
            f'-o={output_path}',
        ]

        # Keep DSSR's stdout and stderr in the log while allowing a non-zero
        # exit status to propagate to the pipeline.
        with log_path.open('w') as log:
            subprocess.run(
                dssr_command,
                cwd=work_directory,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=timeout,
            )

        # Parse only after DSSR has completed successfully and written its output.
        text = output_path.read_text()
        data = json.loads(text)

    # Publish the JSON for analysis before validating it, so it is also kept
    # when the residue mapping rejects it.
    temporary = result_path.with_suffix('.tmp')
    temporary.write_text(text)
    temporary.replace(result_path)
    return parse_pairs(data, sequence)


# Standard residues AF3 writes; own table because Biopython's key case varies by version.
THREE_TO_ONE = dict(zip(
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
    'ACDEFGHIKLMNPQRSTVWY',
))

PROTEIN_DEFAULTS = {'atom': 'CA', 'cutoff': 8.0, 'min_separation': 6}


def contact_pairs(
    model: Union[str, Path],
    sequence: str,
    *,
    chain: str = 'A',
    atom: str = 'CA',
    cutoff: float = 8.0,
    min_separation: int = 6,
) -> list[tuple[int, int]]:
    """Return zero-based residue contacts from one chain of an mmCIF model.

    Positions come from Biopython, the contact map is a numpy boolean mask over
    the pairwise distance matrix. Only i < j with ``j - i >= min_separation`` is
    kept so sequence neighbours are not reported as contacts.
    """
    try:
        import numpy as np
        from Bio.PDB.MMCIFParser import MMCIFParser
    except ImportError as exc:
        raise ImportError(
            'Protein extraction needs biopython and numpy in the controller environment'
        ) from exc

    structure = MMCIFParser(QUIET=True).get_structure('model', str(model))
    chains = {c.id: c for c in structure[0]}
    if chain not in chains:
        raise ValueError(f'Chain {chain} not found in {model}')

    residues = [r for r in chains[chain] if r.id[0] == ' ']
    if len(residues) != len(sequence):
        raise ValueError('Protein residue count does not match input')

    coordinates = []
    for residue, code in zip(residues, sequence):
        if THREE_TO_ONE.get(residue.get_resname().upper()) != code:
            raise ValueError('Unexpected protein residue order or identity')
        if atom not in residue:
            raise ValueError(f'Residue {residue.id[1]} has no {atom} atom')
        coordinates.append(residue[atom].coord)

    coordinates = np.asarray(coordinates, dtype=float)
    distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=-1)
    mask = np.triu(distances <= cutoff, k=min_separation)
    return [(int(i), int(j)) for i, j in np.argwhere(mask)]


def _extract_rna(model, directory, sequence, settings, scratch=None):
    return run_dssr(
        model, directory, sequence,
        settings['dssr'], settings['dssr_timeout_seconds'],
        scratch=scratch,
    )


def _extract_protein(model, directory, sequence, settings, scratch=None):
    options = {**PROTEIN_DEFAULTS, **settings.get('protein', {})}
    pairs = contact_pairs(model, sequence, **options)
    work_directory = Path(directory)
    work_directory.mkdir(parents=True, exist_ok=True)
    result_path = work_directory / 'contacts.json'
    temporary = result_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(dict(options=options, pairs=pairs)))
    temporary.replace(result_path)
    return pairs


EXTRACTORS: dict[str, Callable] = {
    'rna': _extract_rna,
    'protein': _extract_protein,
}


def extract_pairs(
    chain_type: str,
    model: Union[str, Path],
    directory: Union[str, Path],
    sequence: str,
    settings: dict[str, Any],
    scratch: Union[str, Path, None] = None,
) -> list[tuple[int, int]]:
    """Dispatch to the extractor for ``chain_type`` using the 'extraction' config.

    ``directory`` receives the published result (``dssr.json`` or
    ``contacts.json``); ``scratch`` receives logs and temporary files.
    """
    if chain_type not in EXTRACTORS:
        raise ValueError(f'No extractor for chain type: {chain_type}')
    return EXTRACTORS[chain_type](model, directory, sequence, settings, scratch)
