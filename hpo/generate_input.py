"""Small bridge run by the SHS Python interpreter; original generator stays intact.

Use --interactions with two-element lists so the original CLI decodes JSON before
PairMap.from_raw. Supplying base native AF3 JSON supports both old/new process().
process(write=False) avoids long output filenames and swallowed main() errors.
"""
import json
import sys
from pathlib import Path


def main():
    generator_dir, request_file, output_file = map(Path, sys.argv[1:])
    request = json.loads(request_file.read_text())
    sys.path.insert(0, str(generator_dir.resolve()))
    import shs_generator as shs
    base = dict(name=request['name'], modelSeeds=[request['af3_seed']],
                dialect='alphafold3', version=1,
                sequences=[{'rna': dict(id='A', sequence=request['sequence'],
                                       modifications=[], unpairedMsa='')}])
    base_path = output_file.with_suffix('.base.json')
    base_path.write_text(json.dumps(base))
    # Empty interactions must be a truthy CLI string ('[]'), not an empty Python list.
    sys.argv = ['shs_generator.py', '--input_json_path', str(base_path),
                '--interactions', json.dumps(request['pairs']), '--seed', str(request['shs_seed'])]
    args = shs.parse_args()
    for name, value in request['parameters'].items():
        if not hasattr(args, name):
            raise ValueError(f'Unknown generator parameter: {name}')
        setattr(args, name, value)
    result = shs.MsaGenerator(args).process(write=False)
    if result is None:
        raise ValueError('Generator skipped input')
    result['name'] = request['name']
    result['modelSeeds'] = [request['af3_seed']]
    rna = result['sequences'][0]['rna']
    rows = rna['unpairedMsa'].splitlines()
    seqs = [r for r in rows if not r.startswith('>')]
    if not seqs or seqs[0] != request['sequence']:
        raise ValueError('MSA query mismatch')
    if any(len(''.join(c for c in row if not c.islower())) != len(seqs[0]) for row in seqs):
        raise ValueError('MSA aligned lengths differ')
    output_file.write_text(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
