"""The protein generator must describe the chain to AF3 as a protein, not RNA."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_protein_af3_input_has_protein_chain(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / 'SHS-Generator-Protein'))
    sys.modules.pop('json_generator', None)
    try:
        import json_generator
        data = json_generator.build_input_json('MKTA', ['MKTA', 'MKtTA'], name='p', model_seeds=[3])
    finally:
        sys.modules.pop('json_generator', None)
    assert data['modelSeeds'] == [3] and data['name'] == 'p'
    (chain,) = data['sequences']
    assert list(chain) == ['protein']
    protein = chain['protein']
    assert protein['id'] == 'A' and protein['sequence'] == 'MKTA'
    assert protein['unpairedMsa'].startswith('>query\nMKTA')
    assert protein['pairedMsa'] == '' and protein['templates'] == []
