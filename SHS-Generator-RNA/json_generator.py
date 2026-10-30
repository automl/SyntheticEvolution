"""Build native AlphaFold 3 JSON and attach an MSA without mutating input data."""
from copy import deepcopy
from typing import Any, Dict, Sequence


def msa_to_a3m(msa: Sequence[str]) -> str:
    if not msa:
        raise ValueError('MSA must contain a query')
    return '\n'.join(['>query\n' + msa[0]] +
                     [f'>sample_{i}\n{row}' for i, row in enumerate(msa[1:])])


def build_input_json(rna_seq: str, msa=None, *, name='default_name',
                     model_seeds=None) -> Dict[str, Any]:
    """Return an RNA-only native AF3 input (chain A).

    The existing build_input_json(rna_seq) call still creates a query-only MSA.
    The generator validates sequences and alignments before this serialization.
    """
    rna_seq = rna_seq.upper()
    return dict(
        name=name, 
        modelSeeds=[1] if model_seeds is None else list(model_seeds),
        sequences=[{'rna': dict(sequence=rna_seq, modifications=[], 
                    id='A',
                    unpairedMsa=msa_to_a3m([rna_seq] if msa is None else msa))
                    }],
        dialect='alphafold3', version=1,
    )


def attach_rna_msas(data, msas, *, name=None, model_seeds=None):
    """Copy AF3 input and replace MSAs by sequence-entry index.

    Keep protein/DNA/ligand entries and all unrelated metadata. Embedded MSA and
    unpairedMsaPath are alternatives, so remove the old path on updated RNA chains.
    """
    result = deepcopy(data)
    for index, msa in msas.items():
        rna = result['sequences'][index]['rna']
        rna['sequence'] = msa[0]
        rna['unpairedMsa'] = msa_to_a3m(msa)
        rna.pop('unpairedMsaPath', None)
    if name is not None:
        result['name'] = name
    if model_seeds is not None:
        result['modelSeeds'] = list(model_seeds)
    return result
