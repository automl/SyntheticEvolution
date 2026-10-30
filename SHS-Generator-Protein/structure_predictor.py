"""Secondary-structure prediction backends for shs_generator.

Kept separate so shs_generator does not eagerly import RnaBench / tensorflow /
torch: those are only needed when --structure_predictor is used, so
shs_generator imports this module lazily and each backend imports its heavy
dependencies inside the function.

Each predictor returns *raw* predicted pairs in whatever shape its backend
produces (a dict, a list of [i, j], or a list of [i, j, extra]); the caller is
responsible for normalising/converting them.

Note: the spotrna and rnaformer backends use paths relative to the current
working directory (external_algorithms/, data/, RNAformer/), so shs_generator
must be run from the repository root when those predictors are selected.
"""
from pathlib import Path


def predict_rnafold(seq):
    from RnaBench.lib.rna_folding_algorithms.rnafold import RNAFold
    pred_pairs, _ = RNAFold()(seq)
    return pred_pairs


def predict_spotrna(seq):
    from RnaBench.lib.rna_folding_algorithms.DL.spotrna import SpotRna
    return SpotRna()(seq)


def predict_rnaformer(seq, pdb_id):
    import pickle
    import subprocess
    import pandas as pd

    df = pd.read_pickle('data/rnaformer_predictions.pkl')
    d = df[df['pdb_id'].str.lower() == pdb_id]
    if d.empty and df[df['sequence'] == seq]['pairs'].empty:
        sample = {'Id': pdb_id.upper(), 'sequence': seq}
        df = pd.DataFrame([sample])

        with open(Path(f'RNAformer/datasets/{pdb_id.upper()}.pkl'), 'wb') as f:
            pickle.dump(df, f)

        subprocess.run(['./RNAformer/run_cpu.sh', 'infer_RNAformer.py', '-c', '1',
                        '-m', 'models/af3_like_finetune', '-p', f'datasets/{pdb_id.upper()}.pkl'])

        pred = pd.read_pickle(f'RNAformer/datasets/{pdb_id.upper()}_processed.pkl')

        pred_pairs = sorted(pred['pairs'].values[0], key=lambda x: x[0])
    elif not d.empty:
        pred_pairs = df[df['pdb_id'].str.lower() == pdb_id]['pairs'].values[0]
    else:
        pred_pairs = df[df['sequence'] == seq]['pairs'].values[0]
    return [[p1, p2, 0] for p1, p2 in pred_pairs]


def predict(predictor, seq, pdb_id=None):
    """Dispatch to a named predictor and return its raw predicted pairs.
    pdb_id is only used by rnaformer (already lowercased/truncated by the caller)."""
    name = predictor.lower()
    if name == "rnafold":
        return predict_rnafold(seq)
    if name == "spotrna":
        return predict_spotrna(seq)
    if name == "rnaformer":
        return predict_rnaformer(seq, pdb_id)
    raise ValueError(f"Unknown structure predictor: {predictor}")
