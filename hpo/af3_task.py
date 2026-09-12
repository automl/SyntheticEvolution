"""Run inside the AF3 module environment. A success marker requires a model file."""
import json
import os
import subprocess
import sys
from pathlib import Path


def main():
    tasks = json.loads(Path(sys.argv[1]).read_text())
    index = int(os.environ.get('SLURM_ARRAY_TASK_ID', '0'))
    task = tasks[index]
    output = Path(task['output'])
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(Path(os.environ['ALPHAFOLD_BIN_DIR']) / 'run_alphafold.py'),
                    '--json_path=' + task['input'], '--output_dir=' + str(output),
                    '--model_dir=' + task['model_dir'], '--db_dir=' + os.environ['ALPHAFOLD_DATABASES'],
                    '--norun_data_pipeline'], check=True)
    models = list(output.glob('*/*_model.cif'))
    if len(models) != 1 or models[0].stat().st_size == 0:
        raise RuntimeError(f'Expected one top-level model, found {models}')
    marker = output / 'success.json'
    temporary = marker.with_suffix('.tmp')
    temporary.write_text(json.dumps({'model': str(models[0].resolve())}))
    temporary.replace(marker)


if __name__ == '__main__':
    main()
