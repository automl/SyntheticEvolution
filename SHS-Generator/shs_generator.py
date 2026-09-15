#!/usr/bin/env python
"""Generate synthetic RNA MSAs with a small explicit Python API.

MsaGenerator only operates on sequences, PairMap and MutationParameters.
CLI helpers below the class handle requests, AF3 JSON and files.
"""
import argparse
from dataclasses import fields
import json
import logging
import random
import sys
from pathlib import Path
from typing import Dict, List
import numpy as np
from typing import Optional

# Retain the repository import path for optional structure predictors.
ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from pair_map import PairMap
from generator_config import (
    MutationParameters, APPROACHES, build_pair_map, split_parameters,
    validate_seed, validate_sequence,
)

PAIR_MUTATION_PROBABILITIES = {
    "A": {"U": 1.0},
    "C": {"G": 1.0},
    "G": {"C": 0.75, "U": 0.25},
    "U": {"A": 0.75, "G": 0.25},
}

PAIR_MUTATIONS: Dict[str, List[str]] = {
    'CU': ['GU', 'AU', 'CG'],
    'CA': ['UA', 'CG'],
    'GA': ['UA', 'GC', 'GU'],
    'CC': ['CG', 'GC'],
    'AA': ['AU', 'UA'],
    'UU': ['AU', 'UA'],
    'GG': ['GC', 'CG', 'GU', 'UG'],
    'GC': ['AU', 'CG', 'GC'],
    'CG': ['GC', 'AU'],
    'AU': ['UA', 'GC'],
    'UA': ['AU', 'CG'],
    'GU': ['GC', 'AU', 'UG'],
    'UG': ['GC', 'AU', 'GU']
}

WC_PAIRS: List[str] = ['AU', 'UA', 'GC', 'CG']

class MsaGenerator:
    """Generate MSAs with an instance-local random stream and no file I/O.

    Make a new instance with the same seed to reproduce the same MSA. Consecutive
    generate() calls on one instance advance its stream (also useful for chains).
    PairMap contains per-position rates; parameters contains the remaining settings.
    """
    def __init__(self, parameters: Optional[MutationParameters] = None, seed=1):
        self.parameters = MutationParameters() if parameters is None else parameters
        if not isinstance(self.parameters, MutationParameters):
            raise TypeError('parameters must be MutationParameters, not argparse.Namespace')
        validate_seed(seed)
        self.rng = random.Random(None if seed is None else int(seed))
        self.pair_map = None

    def generate(self, rna_sequence: str, pair_map: PairMap) -> List[str]:
        """Return query followed by N-1 synthetic A3M sequence rows (no headers)."""
        sequence = validate_sequence(rna_sequence)
        if not isinstance(pair_map, PairMap):
            raise TypeError('pair_map must be a PairMap')
        matrix = pair_map._pairs_mat
        if matrix.shape != (len(sequence), len(sequence)):
            raise ValueError('PairMap size must match RNA length')
        if not np.all(np.isfinite(matrix)) or np.any(matrix < 0) or np.any(matrix > 1):
            raise ValueError('PairMap rates and interactions must be finite and in [0,1]')
        if not np.array_equal(matrix, matrix.T):
            raise ValueError('PairMap must be symmetric')
        self.pair_map = pair_map
        self.max_insertion_length = max(int(len(sequence) * self.parameters.max_insertion_fraction), 2)
        self.max_deletion_length = max(int(len(sequence) * self.parameters.max_deletion_fraction), 2)
        msa = [sequence] + [self.mutate_sequence(sequence) for _ in range(self.parameters.N - 1)]
        validate_msa(msa, sequence, self.parameters.N)
        return msa

    def loop_insertion(self) -> str:
        """Long insertions take priority, then single insertions. If all long insertions are disregarded
        the single insertion rate can be recovered accurately. Any insertion is randomly selected"""
        if self.rng.random() < self.parameters.loop_long_insertion_prob:
            insertion_len = self.rng.randint(2, self.max_insertion_length)
            return ''.join(self.rng.choice('augc') for _ in range(insertion_len))
        if self.rng.random() < self.parameters.loop_single_insertion_prob:
            return self.rng.choice('acgu')
        return ""

    def stem_insertion(self) -> str:
        """Long insertions take priority, then single insertions. If all long insertions are disregarded
        the single insertion rate can be recovered accurately. Any insertion is randomly selected"""
        if self.rng.random() < self.parameters.stem_long_insertion_prob:
            insertion_len = self.rng.randint(2, self.max_insertion_length)
            return ''.join(self.rng.choice('augc') for _ in range(insertion_len))
        if self.rng.random() < self.parameters.stem_single_insertion_prob:
            return self.rng.choice('acgu')
        return ""

    def paired_deletion(self, i, mutated_partner: str, partner_index: float) -> bool:
        single_del = self.parameters.stem_single_deletion_prob
        pair_del =  self.parameters.stem_pair_deletion_prob * self.pair_map.interaction(i, partner_index)
        if self.pair_map.is_multiplet_member(i):
            return self.rng.random() < single_del
        if i < partner_index: # P(-)
            return self.rng.random() < single_del + pair_del
        if mutated_partner == "-": # P(-|-)
            return self.rng.random() * (single_del + pair_del) < pair_del
        return self.rng.random() * (1 - single_del - pair_del) < single_del # P(-|!-)

    def mutate_random(self, nt: str, prob: float):
        if self.rng.random() < prob:
            return self.rng.choice([c for c in 'AUGC' if c != nt])
        return nt
    
    def mutate_unpaired(self, nt: str, loop_long_del_len: int, mutation_rate: float) -> tuple[str, int]:
        """Long deletions take priority, then single deletions, then mutations. Therefore
        the mutation rate is can be recovered accurately if deletions are disregarded."""
        if self.rng.random() < self.parameters.loop_long_deletion_prob:
            loop_long_del_len = self.rng.randint(2, self.max_deletion_length)
        if loop_long_del_len > 0:
            return "-", loop_long_del_len - 1
        if self.rng.random() < self.parameters.loop_single_deletion_prob:
            return "-", 0
        return self.mutate_random(nt, mutation_rate), 0

    def mutate_cov(self, nt: str, partners_original: np.ndarray, partners_mutated: np.ndarray,
                    interactions: np.ndarray, partner_mutation_rates: np.ndarray, mutation_rate: float) -> str:
        """Guarantees that all partners get mutated together but randomly and therefore only increases covariance."""
        opts = list(zip((partners_original != partners_mutated), interactions, partner_mutation_rates))
        m = mutation_rate
        prob = m
        for mut, inter, pm in opts:
            if mut:
                prob += ((m * pm) + inter * m * (1 - pm)) / max(pm, 0.0001)
            else:
                prob += ((m * (1 - pm)) - inter * m * (1 - pm)) / (1 - min(pm, 0.9999))
        return self.mutate_random(nt, prob / (len(opts)+1))

    def mutate_wc(self, nt: str, partners_original: np.ndarray, partners_mutated: np.ndarray, 
                  interactions: np.ndarray, partner_mutation_rates: np.ndarray, mutation_rate: float, increase_cov: bool) -> str:
        """Try to maximize the number of watson crick base pairs while also trying to leave no partner.
        This uses the probabilities in PAIR_MUTATION_PROBABILITIES."""
        # first in the multiplet is random with mutation rate
        if len(partners_original) == 0:
            return self.mutate_random(nt, mutation_rate)
        # use mutate_cov to decide if it should mutate for high covariance or decide randomly if nt should mutate for low covariance
        if increase_cov:
            mutate = self.mutate_cov(nt, partners_original, partners_mutated, interactions, partner_mutation_rates, mutation_rate) != nt
        else:
            mutate = self.mutate_random(nt, mutation_rate) != nt
        if not mutate:
            return nt
        # chose with interactions, always mutate
        options = {"A": 0.00001, "U": 0.00001, "G": 0.00001, "C": 0.00001}
        for j, mut in enumerate(partners_mutated):
            for opt, prob in PAIR_MUTATION_PROBABILITIES.get(mut, {}).items():
                options[opt] += prob * interactions[j]
        options.pop(nt)
        return self.rng.choices(list(options.keys()), list(options.values()))[0]
        

    def mutate_pair_original(self, p_nt: str, nt: str, mutation_rate: float, partner_mutation_rate: float) -> str:
        if not self.rng.random() < (mutation_rate + partner_mutation_rate) / 2:
            return p_nt + nt
        if self.rng.random() < self.parameters.wobble_prob:
            return self.rng.choice(['GU', 'UG'])
        candidates = PAIR_MUTATIONS.get(p_nt + nt, WC_PAIRS)
        return self.rng.choice(candidates)

    def mutate_sequence(self, seq: str) -> str:
        approach = self.parameters.pair_mutation_approach
        if approach not in ["covariance", "watson_crick", "watson_crick_cov", "original", "none"]:
            logging.error("Unknown input for --pair-mutation-approach, '%s', please use one of the provided options", self.parameters.pair_mutation_approach)
            raise ValueError()
        seq = np.array(list(seq))
        new_seq = np.empty(len(seq), str)
        insertions = []
        loop_long_del_len = 0
        for i, nt in enumerate(seq):
            new_nt = nt
            new_insertion = ""
            mutation_rate = self.pair_map.mutation_rate(i)

            if self.pair_map.is_unpaired(i):
                new_nt, loop_long_del_len = self.mutate_unpaired(nt, loop_long_del_len, mutation_rate)
            
            if self.pair_map.is_paired(i):
                loop_long_del_len = 0
                partners =  np.array(self.pair_map.partners(i), int)
                prev = partners[partners < i]
                interactions = np.array([self.pair_map.interaction(i, j) for j in prev])
                partner_mutation_rates = np.array([self.pair_map.mutation_rate(j) for j in prev])
                if approach == "none":
                    new_nt = self.mutate_random(nt, self.pair_map.mutation_rate(i))
                if approach == "covariance":
                    new_nt = self.mutate_cov(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates,mutation_rate)
                if approach == "watson_crick":
                    new_nt = self.mutate_wc(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates, mutation_rate,  False)
                if approach == "watson_crick_cov":
                    new_nt = self.mutate_wc(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates, mutation_rate, True)
                if approach == "original" :
                    if self.pair_map.is_basic_pair(i) and i > partners[0]:
                        new_seq[partners[0]], new_nt = self.mutate_pair_original(seq[partners[0]], nt, mutation_rate, partner_mutation_rates[0])
                    if self.pair_map.is_multiplet_member(i):
                        new_nt = self.mutate_wc(nt, seq[prev], new_seq[prev], interactions, partner_mutation_rates, mutation_rate, True)
                # optionally override mutation if a deletion happens
                if self.paired_deletion(i, new_seq[partners[0]], partners[0]):
                    new_nt = "-"

            if self.pair_map.is_paired(i-1) and self.pair_map.is_paired(i):
                new_insertion = self.stem_insertion()
            else:
                new_insertion = self.loop_insertion()

            new_seq[i] = new_nt
            insertions.append(new_insertion)
        return ''.join(np.ravel(list(zip(insertions, new_seq)))) + self.loop_insertion()


# Input/output helpers: used by main(), never by MsaGenerator.

def validate_msa(msa, sequence, count):
    """Validate the emitted A3M representation before writing an AF3 input."""
    if len(msa) != count or not msa or msa[0] != sequence:
        raise ValueError('MSA query or sequence count mismatch')
    for row in msa:
        if set(row) - set('ACGUacgu-'):
            raise ValueError('MSA contains invalid characters')
        if len(''.join(c for c in row if not c.islower())) != len(sequence):
            raise ValueError('MSA aligned lengths differ')


def parse_json(text):
    """Accept JSON arrays and the legacy tuple-style CLI syntax."""
    return json.loads(text.replace('(', '[').replace(')', ']'))


def parse_args(argv=None):
    # SUPPRESS lets us distinguish explicit options from absent defaults.
    parser = argparse.ArgumentParser(
        description='Generate an RNA MSA and native AF3 input JSON.',
        argument_default=argparse.SUPPRESS,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--request-json', help='Pipeline request with sequence, pairs, parameters and seeds')
    source.add_argument('--rna-seq', help='RNA sequence; supply structure/interactions/predictor too')
    source.add_argument('--input-json-path', help='Existing native AF3 JSON to adapt')
    output = parser.add_mutually_exclusive_group()
    output.add_argument('--output-json', help='Exact output path; parent directories are created')
    output.add_argument('--output-json-dir', help='Output directory; defaults to custom_msa_json_output')
    parser.add_argument('--structure', help='Dot-bracket or zero-based JSON pair list')
    parser.add_argument('--interactions', help='JSON pairs with optional interaction strengths')
    parser.add_argument('--structure-predictor', choices=('rnafold', 'spotrna', 'rnaformer', 'dssr'))
    parser.add_argument('--mutation-rate-paired', type=float)
    parser.add_argument('--mutation-rate-unpaired', type=float)
    parser.add_argument('--mutation-rates', help='JSON per-position mutation rates')
    for field in fields(MutationParameters):
        name = field.name
        if name == 'N':
            parser.add_argument('-N', type=int)
        elif name == 'pair_mutation_approach':
            parser.add_argument('--pair-mutation-approach', choices=APPROACHES)
        else:
            flag = name.replace('_', '-')
            parser.add_argument('--' + flag, type=float)
    parser.add_argument('--seed', type=int, help='SHS seed; default None for legacy CLI')
    parser.add_argument('--af3-seed', type=int, help='Override AF3 model seed; default 1 for new inputs')
    parser.add_argument('--task-name', help='Short output/model name; overrides legacy parameter-based name')
    parser.add_argument('--pdb-id')
    parser.add_argument('--protein-seq', help='Deprecated; ignored (new direct inputs are RNA-only)')
    parser.add_argument('--pair-mutation', help='Deprecated unused option; retained for CLI compatibility')
    parser.add_argument('--max-chains', type=int)
    
    parser.add_argument('--plot', action='store_true')      # deprecated; will be moved to analyse
    parser.add_argument('--print-msa', action='store_true') # deprecated; will be moved to analyse
    parser.add_argument('--show-plot', action='store_true') # deprecated; will be moved to analyse
    args = parser.parse_args(argv)
    if hasattr(args, 'request_json'):
        permitted = {'request_json', 'output_json', 'output_json_dir'}
        conflicts = set(vars(args)) - permitted
        if conflicts:
            parser.error('--request-json cannot be combined with generation options: ' + ', '.join(sorted(conflicts)))
    return args


def valid_task_name(name):
    # This name may become a filename when output-json is omitted.
    if not isinstance(name, str) or not name or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in name):
        raise ValueError('task_name must contain only letters, digits, underscore or hyphen')
    return name


def load_request(path):
    """Validate the pipeline request and produce explicit generator inputs."""
    request = json.loads(Path(path).read_text())
    if not isinstance(request, dict):
        raise ValueError('Request must be a JSON object')
    required = {'task_name', 'sequence', 'pairs', 'parameters', 'shs_seed', 'af3_seed'}
    if required - request.keys():
        raise ValueError(f'Missing request fields: {sorted(required - request.keys())}')
    if request.keys() - required - {'id'}:
        raise ValueError(f'Unknown request fields: {sorted(request.keys() - required - {"id"})}')
    valid_task_name(request['task_name'])
    validate_seed(request['shs_seed'], 'shs_seed', allow_none=False)
    validate_seed(request['af3_seed'], 'af3_seed', allow_none=False)
    sequence = validate_sequence(request['sequence'])
    pairs = request['pairs']
    if not isinstance(pairs, list) or any(not isinstance(p, list) or len(p) != 2 for p in pairs):
        raise ValueError('Request pairs must be a list of two-element lists')
    parameters, paired, unpaired = split_parameters(request['parameters'])
    pair_map = build_pair_map(sequence, pairs, paired, unpaired)
    return request, sequence, pair_map, parameters


def prepare_structure(sequence, options, paired, unpaired):
    """Resolve legacy CLI structure selection into the same PairMap used by requests."""
    if 'interactions' in options:
        if 'structure' in options or 'structure_predictor' in options:
            logging.warning('Using interactions instead of structure/predictor')
        raw = parse_json(options['interactions'])
    elif 'structure' in options:
        if 'structure_predictor' in options:
            logging.warning('Using provided structure instead of predictor')
        text = options['structure']
        try:
            raw = parse_json(text)
        except json.JSONDecodeError:
            raw = text
    elif 'structure_predictor' in options:
        import structure_predictor
        pdb = options.get('pdb_id')
        raw = structure_predictor.predict(options['structure_predictor'], sequence, pdb.lower()[:4] if pdb else None)
    else:
        raise ValueError('Provide --structure, --interactions or --structure_predictor')
    rates = parse_json(options['mutation_rates']) if 'mutation_rates' in options else None
    return build_pair_map(sequence, raw, paired, unpaired, rates)


def legacy_output_name(options, parameters, sequence, paired, unpaired):
    """Preserve the old CLI naming convention; explicit output-json avoids it."""
    p = parameters
    rates = 'custom' if 'mutation_rates' in options else None
    parts = [
        f'{options.get("pdb_id")}_custom_rnamsa', f'N{p.N}', f'seed{options.get("seed")}',
        f'mru_{unpaired}', f'mrp_{paired}', f'mrs_{rates}', f'pma_{p.pair_mutation_approach}',
        f'ssi_{p.stem_single_insertion_prob}', f'sli_{p.stem_long_insertion_prob}',
        f'spd_{p.stem_pair_deletion_prob}', f'lsi_{p.loop_single_insertion_prob}',
        f'lsd_{p.loop_single_deletion_prob}', f'lli_{p.loop_long_insertion_prob}',
        f'lld_{p.loop_long_deletion_prob}', f'mif_{p.max_insertion_fraction}',
        f'mdf_{p.max_deletion_fraction}', f'maxinslen_{max(int(len(sequence)*p.max_insertion_fraction),2)}',
        f'maxdellen_{max(int(len(sequence)*p.max_deletion_fraction),2)}',
        f'wp_{p.wobble_prob}', options.get('structure_predictor') or 'none',
    ]
    return '_'.join(parts)


def show_generation(msa, pair_map, options):
    if options.get('print_msa'):
        for row in msa:
            logging.info('%s', row)
    if options.get('plot'):
        import msa_plotting
        msa_plotting.plot_final_features(msa, msa[0], pair_map.pairs,
                                        options.get('pdb_id'), options.get('show_plot', False))


def prepare_af3_json(args):
    """Read one input mode, run the pure generator, and build JSON in memory."""
    import json_generator
    options = vars(args)
    if 'request_json' in options:
        request, sequence, pair_map, parameters = load_request(options['request_json'])
        msa = MsaGenerator(parameters, seed=request['shs_seed']).generate(sequence, pair_map)
        show_generation(msa, pair_map, options)
        return json_generator.build_input_json(sequence, msa, name=request['task_name'],
                                              model_seeds=[request['af3_seed']])

    if options.get('protein_seq'):
        logging.warning('Custom protein sequences are not supported; direct input is RNA-only')
    parameter_names = {f.name for f in fields(MutationParameters)} | {'mutation_rate_paired', 'mutation_rate_unpaired'}
    parameters, paired, unpaired = split_parameters({k:v for k,v in options.items() if k in parameter_names})
    seed = options.get('seed')
    validate_seed(seed)
    af3_seed = options.get('af3_seed')
    if af3_seed is not None:
        validate_seed(af3_seed, 'af3_seed', allow_none=False)
    if 'rna_seq' in options:
        data = json_generator.build_input_json(validate_sequence(options['rna_seq']))
    else:
        data = json.loads(Path(options['input_json_path']).read_text())
        if not isinstance(data, dict) or data.get('dialect') != 'alphafold3' or not isinstance(data.get('sequences'), list):
            raise ValueError('input_json_path must contain native AF3 JSON with sequences')
    if 'max_chains' in options:
        if options['max_chains'] < 1 or len(data['sequences']) > options['max_chains']:
            raise ValueError('Input exceeds --max_chains or limit is invalid')
    generator = MsaGenerator(parameters, seed=seed)
    msas = {}
    for i, chain in enumerate(data['sequences']):
        if 'rna' not in chain:
            continue
        sequence = validate_sequence(chain['rna']['sequence'])
        pair_map = prepare_structure(sequence, options, paired, unpaired)
        msas[i] = generator.generate(sequence, pair_map)
        show_generation(msas[i], pair_map, options)
    if not msas:
        raise ValueError('Input contains no RNA chains')
    if 'task_name' in options:
        name = valid_task_name(options['task_name'])
    elif 'output_json' in options:
        # Short explicit filename also supplies the name for non-request CLI mode.
        name = Path(options['output_json']).stem
    else:
        name = legacy_output_name(options, parameters, sequence, paired, unpaired)
    return json_generator.attach_rna_msas(data, msas, name=name,
                                         model_seeds=[af3_seed] if af3_seed is not None else None)


def write_output(data, args):
    """Write exactly one final AF3 JSON; no intermediate base JSON is needed."""
    options = vars(args)
    if 'output_json' in options:
        path = Path(options['output_json'])
    else:
        path = Path(options.get('output_json_dir', 'custom_msa_json_output')) / (data['name'] + '.json')
    for name in ('request_json', 'input_json_path'):
        if name in options and path.resolve() == Path(options[name]).resolve():
            raise ValueError('Output must not overwrite the input/request file')
    path.parent.mkdir(parents=True, exist_ok=True)
    # On failure, no partially written final output appears as a successful result.
    temporary = path.with_suffix(path.suffix + '.tmp')
    try:
        temporary.write_text(json.dumps(data, indent=2, allow_nan=False))
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    logging.info('JSON written to: %s', path)
    return path


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = parse_args(argv)
    try:
        data = prepare_af3_json(args)
        write_output(data, args)
    except Exception as exc:
        logging.error('Generation failed: %s', exc)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
