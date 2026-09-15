"""Contract tests for the NEW shs_generator.py request-file interface.

Place in SyntheticEvolution/hpo/tests/. Run with pytest; subprocesses use
SHS_PYTHON when set, otherwise the Python running pytest. Optional SHS_GENERATOR
can point at a checkout elsewhere. No Slurm, AF3 inference or DSSR is called.
Includes request/legacy CLI checks, pure API checks, and seeded regression fixtures.
"""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

SEQUENCE = 'AGUAGUAGUAGUAGA'


@pytest.fixture(scope='module')
def generator_cli():
    """Require the real generator and new CLI before testing error cases.

    This prevents an old CLI rejecting every request from making negative tests
    appear to pass. Missing dependencies also produce a useful setup failure.
    """
    script = Path(os.environ.get(
        'SHS_GENERATOR',
        str(Path(__file__).resolve().parents[2] / 'SHS-Generator/shs_generator.py'),
    )).resolve()
    python = os.environ.get('SHS_PYTHON', sys.executable)
    assert script.is_file(), f'Generator not found: {script}; set SHS_GENERATOR if needed'
    cli = [python, str(script)]
    result = subprocess.run(cli + ['--help'], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert '--request-json' in result.stdout and '--output-json' in result.stdout, (
        'Implement --request-json and --output-json in shs_generator.py first.'
    )
    return cli


@pytest.fixture
def request_data():
    """Small RNA-only input; indels disabled to isolate mutation behavior."""
    return dict(
        id='dataset-row-1', task_name='rna_00000', sequence=SEQUENCE,
        pairs=[[0,14], [1,13], [7,12]], shs_seed=42, af3_seed=17,
        parameters=dict(
            N=6, mutation_rate_paired=0.2, mutation_rate_unpaired=0.2,
            pair_mutation_approach='watson_crick_cov',
            stem_single_insertion_prob=0.0, stem_long_insertion_prob=0.0,
            stem_single_deletion_prob=0.0, stem_pair_deletion_prob=0.0,
            loop_single_insertion_prob=0.0, loop_long_insertion_prob=0.0,
            loop_single_deletion_prob=0.0, loop_long_deletion_prob=0.0,
        ),
    )


@pytest.fixture
def invoke(generator_cli, tmp_path):
    """Launch the real CLI from an unrelated cwd, with paths containing spaces."""
    def run(request, case='case'):
        folder = tmp_path / case / 'task with spaces'
        folder.mkdir(parents=True)
        request_path = folder / 'shs_generation_request.json'
        output_path = folder / 'af3_input.json'
        text = request if isinstance(request, str) else json.dumps(request)
        request_path.write_text(text)
        result = subprocess.run(
            generator_cli + ['--request-json', str(request_path),
                             '--output-json', str(output_path)],
            cwd=tmp_path, capture_output=True, text=True, timeout=60,
        )
        # Generation must not modify the pipeline's request record.
        assert request_path.read_text() == text
        return result, output_path
    return run


def read_success(result, path):
    assert result.returncode == 0, result.stdout + result.stderr
    assert path.is_file(), f'Generator did not write requested output: {path}'
    return json.loads(path.read_text())


def msa_sequences(data):
    """Parse FASTA/A3M records, allowing sequences wrapped across multiple lines."""
    records = []
    for line in data['sequences'][0]['rna']['unpairedMsa'].splitlines():
        if not line.strip():
            continue
        if line.startswith('>'):
            records.append('')
        else:
            assert records, 'MSA sequence before first FASTA header'
            records[-1] += line.strip()
    assert records and all(records), 'MSA contains no sequences or empty records'
    return records


def test_request_cli_writes_rna_only_af3_input(invoke, request_data):
    data = read_success(*invoke(request_data))
    assert data['name'] == request_data['task_name']
    assert data['modelSeeds'] == [17]  # AF3 seed must not be replaced by SHS seed 42.
    assert data['dialect'] == 'alphafold3'
    assert data['version'] == 1
    assert len(data['sequences']) == 1
    assert set(data['sequences'][0]) == {'rna'}
    rna = data['sequences'][0]['rna']
    assert rna['id'] == 'A'
    assert rna['sequence'] == SEQUENCE
    assert rna['modifications'] == []
    records = msa_sequences(data)
    assert len(records) == request_data['parameters']['N']
    assert records[0] == SEQUENCE
    for row in records:
        assert set(row) <= set('ACGUacgu-')
        assert len(''.join(c for c in row if not c.islower())) == len(SEQUENCE)


def test_empty_pair_list_is_valid(invoke, request_data):
    request_data['pairs'] = []
    data = read_success(*invoke(request_data))
    assert msa_sequences(data)[0] == SEQUENCE


def test_n_one_produces_query_only(invoke, request_data):
    request_data['parameters']['N'] = 1
    assert msa_sequences(read_success(*invoke(request_data))) == [SEQUENCE]


def test_zero_mutation_rates_preserve_all_sequences(invoke, request_data):
    request_data['parameters'].update(mutation_rate_paired=0.0, mutation_rate_unpaired=0.0)
    assert msa_sequences(read_success(*invoke(request_data))) == [SEQUENCE] * 6


def test_pair_list_controls_which_positions_mutate(invoke, request_data):
    # With approach none and rates 1/0, every paired position changes and every
    # unpaired position stays identical. This detects ignored or misparsed pairs.
    request_data['parameters'].update(
        pair_mutation_approach='none', mutation_rate_paired=1.0,
        mutation_rate_unpaired=0.0,
    )
    paired = {pos for pair in request_data['pairs'] for pos in pair}
    records = msa_sequences(read_success(*invoke(request_data)))
    for row in records[1:]:
        assert len(row) == len(SEQUENCE)
        for pos, (original, generated) in enumerate(zip(SEQUENCE, row)):
            assert (generated != original) == (pos in paired)


def test_same_request_and_seeds_are_reproducible(invoke, request_data):
    first = read_success(*invoke(request_data, 'first'))
    second = read_success(*invoke(copy.deepcopy(request_data), 'second'))
    assert first == second


@pytest.mark.parametrize('problem', ['missing_sequence', 'unknown_parameter',
                                    'invalid_pair', 'invalid_n', 'malformed_json'])
def test_invalid_request_fails_without_output(invoke, request_data, problem):
    if problem == 'missing_sequence':
        del request_data['sequence']
    elif problem == 'unknown_parameter':
        request_data['parameters']['mutation_rate_typo'] = 0.2
    elif problem == 'invalid_pair':
        request_data['pairs'] = [[0, len(SEQUENCE)]]
    elif problem == 'invalid_n':
        request_data['parameters']['N'] = 0
    else:
        request_data = '{not valid JSON'
    result, output = invoke(request_data)
    assert result.returncode != 0, 'Invalid request was accepted'
    assert not output.exists(), 'Failure left an output that could look successful'
    assert (result.stdout + result.stderr).strip(), 'Failure needs a diagnostic'

# Tests below cover the pure Python API and coexistence with the legacy CLI.
# All generator code still runs in SHS_PYTHON, not in pytest's environment.

@pytest.fixture
def run_api(generator_cli, tmp_path):
    def run(code, *arguments):
        preamble = (
            'import sys, json, random\n'
            'sys.path.insert(0, ' + repr(str(Path(generator_cli[1]).parent)) + ')\n'
            'from shs_generator import MsaGenerator\n'
            'from generator_config import MutationParameters, build_pair_map\n'
        )
        return subprocess.run([generator_cli[0], '-c', preamble + code, *arguments],
                              cwd=tmp_path, text=True, capture_output=True, timeout=60)
    return run


def test_python_api_has_no_file_or_global_rng_side_effects(run_api, tmp_path):
    result = run_api('''
sequence = 'AGUAGUAGUAGUAGA'
parameters = MutationParameters(N=10)
pairs = build_pair_map(sequence, [[0,14],[1,13]], .3, .1)
state = random.getstate()
first = MsaGenerator(parameters, seed=42).generate(sequence, pairs)
assert random.getstate() == state
assert first == MsaGenerator(parameters, seed=42).generate(sequence, pairs)
assert len(first) == 10 and first[0] == sequence
''')
    assert result.returncode == 0, result.stdout + result.stderr
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('case_index', range(10))
def test_seeded_output_matches_original_generator(run_api, case_index):
    fixture = Path(__file__).with_name('fixtures') / 'generator_regression.json'
    result = run_api('''
from pathlib import Path
case = json.loads(Path(sys.argv[1]).read_text())[int(sys.argv[2])]
parameters = MutationParameters(**case['parameters'])
pairs = build_pair_map(case['sequence'], case['pairs'],
                      case['mutation_rate_paired'], case['mutation_rate_unpaired'])
msa = MsaGenerator(parameters, seed=case['seed']).generate(case['sequence'], pairs)
assert msa == case['msa'], case['id']
''', str(fixture.resolve()), str(case_index))
    assert result.returncode == 0, result.stdout + result.stderr


def test_per_position_rates_and_parameter_validation(run_api):
    result = run_api('''
parameters = MutationParameters(N=4, pair_mutation_approach='none',
    stem_single_insertion_prob=0, stem_long_insertion_prob=0,
    stem_single_deletion_prob=0, stem_pair_deletion_prob=0,
    loop_single_insertion_prob=0, loop_long_insertion_prob=0,
    loop_single_deletion_prob=0, loop_long_deletion_prob=0)
sequence = 'ACGU'
pairs = build_pair_map(sequence, [], mutation_rates=[0,1,0,1])
msa = MsaGenerator(parameters, seed=3).generate(sequence, pairs)
for row in msa[1:]:
    assert [a != b for a,b in zip(sequence,row)] == [False,True,False,True]
for kwargs in ({'N':0}, {'N':True}, {'wobble_prob':float('nan')},
               {'max_insertion_fraction':-1}, {'pair_mutation_approach':'typo'}):
    try: MutationParameters(**kwargs)
    except ValueError: pass
    else: raise AssertionError(kwargs)
''')
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('input_args', [
    ['--structure', '((.....(....)))'],
    ['--structure', '[[0,14],[1,13],[7,12]]'],
    ['--interactions', '[]'],
])
def test_legacy_direct_input_with_explicit_output(generator_cli, tmp_path, input_args):
    path = tmp_path / 'nested' / 'short.json'
    result = subprocess.run(generator_cli + ['--rna-seq', SEQUENCE, '-N', '1',
        '--seed', '42', '--output-json', str(path)] + input_args,
        cwd=tmp_path, capture_output=True, text=True, timeout=60)
    data = read_success(result, path)
    assert data['name'] == 'short'
    assert msa_sequences(data) == [SEQUENCE]


def test_legacy_output_directory(generator_cli, tmp_path):
    folder = tmp_path / 'old output'
    result = subprocess.run(generator_cli + ['--rna-seq', SEQUENCE, '--interactions', '[]',
        '-N', '1', '--output-json-dir', str(folder)],
        cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    files = list(folder.glob('*.json'))
    assert len(files) == 1
    assert '_custom_rnamsa_N1_' in files[0].name


def test_existing_af3_input_preserves_other_entities(generator_cli, tmp_path):
    data = dict(name='old', dialect='alphafold3',version=1, modelSeeds=[9],
        sequences=[{'rna':dict(id='A',sequence=SEQUENCE,modifications=[],
                               unpairedMsaPath='unused.a3m')},
                   {'ligand':dict(id='B',ccdCodes=['MG'])}])
    source, output = tmp_path / 'source.json', tmp_path / 'out.json'
    source.write_text(json.dumps(data))
    result = subprocess.run(generator_cli + ['--input-json-path',str(source),
        '--interactions','[]','-N','1','--output-json',str(output)],
        cwd=tmp_path,capture_output=True,text=True,timeout=60)
    out = read_success(result,output)
    assert out['sequences'][1] == data['sequences'][1]
    assert out['modelSeeds'] == [9]
    assert 'unpairedMsaPath' not in out['sequences'][0]['rna']
    assert json.loads(source.read_text()) == data


@pytest.mark.parametrize('extra', [
    ['--rna-seq', SEQUENCE], ['-N', '3'], ['--seed', '9'],
    ['--output_json_dir', 'other-output'],
])
def test_request_rejects_conflicting_cli_options(generator_cli, tmp_path, request_data, extra):
    request = tmp_path / 'request.json'
    output = tmp_path / 'out.json'
    request.write_text(json.dumps(request_data))
    result = subprocess.run(generator_cli + ['--request-json',str(request),
        '--output-json',str(output)] + extra, cwd=tmp_path,
        capture_output=True,text=True,timeout=60)
    assert result.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize('structure', ['((...', '((.............', '..............)', '[[0,15]]'])
def test_invalid_structure_fails_before_writing(generator_cli, tmp_path, structure):
    output = tmp_path / 'out.json'
    result = subprocess.run(generator_cli + ['--rna-seq', SEQUENCE,
        '--structure', structure,'--output-json',str(output)],
        cwd=tmp_path,capture_output=True,text=True,timeout=60)
    assert result.returncode != 0
    assert not output.exists()
