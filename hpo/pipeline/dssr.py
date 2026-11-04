"""Run DSSR in isolation; map its full pairs list through validated nucleotide IDs.

Supports the initial contract: one unmodified RNA, chain A, AF3 residue numbering
1..N. Reject missing/reordered/mismatched residues rather than silently shifting
positions. Verified against the supplied DSSR 2.9.3 output.
"""
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Union


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

DSSR
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
