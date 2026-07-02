"""Build AlphaFold3 input JSON."""
from typing import Any, Dict


def build_input_json(rna_seq: str, protein_seq: str) -> Dict[str, Any]:
    """Return an AlphaFold3-format input dict for an RNA (chain A) + protein
    (chain B) complex. The RNA MSA is seeded with the query row only; the name
    is a placeholder that callers typically overwrite."""
    rna_seq = rna_seq.upper()
    protein_seq = protein_seq.upper()
    unpaired_msa = f">query\n{rna_seq}"
    return {
        "name": "default_name",
        "modelSeeds": [1],
        "sequences": [
            {
                "rna": {
                    "sequence": rna_seq,
                    "modifications": [],
                    "unpairedMsa": unpaired_msa,
                    "id": "A",
                }
            },
            {
                "protein": {
                    "sequence": protein_seq,
                    "modifications": [],
                    "id": "B",
                }
            },
        ],
        "dialect": "alphafold3",
        "version": 1,
    }
