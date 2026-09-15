"""Run inside the AF3 module environment on a GPU node. 
This file executes a single af3 prediction on the given input."""
import json
import os
import subprocess
import sys
from pathlib import Path


def main():
    """Run the AF3 task selected by the current Slurm array index."""
    af3_tasks_path = Path(sys.argv[1])
    tasks = json.loads(af3_tasks_path.read_text())
    index = int(os.environ.get('SLURM_ARRAY_TASK_ID', '0'))
    task = tasks[index]

    output = Path(task['output'])
    output.mkdir(parents=True, exist_ok=True)

    # These environment variables are provided by the AF3 module.
    alphafold_script = Path(os.environ['ALPHAFOLD_BIN_DIR']) / 'run_alphafold.py'
    database_dir = os.environ['ALPHAFOLD_DATABASES']
    alphafold_command = [
        sys.executable,  # Select the active AF3 Python environment.
        str(alphafold_script),
        '--json_path=' + task['input'],
        '--output_dir=' + str(output),
        '--model_dir=' + task['model_dir'],
        '--db_dir=' + database_dir,
        '--norun_data_pipeline',
    ]

    # A non-zero AF3 exit status propagates as an exception, so the Slurm task
    # is marked failed and the controller will not score a partial result.
    subprocess.run(alphafold_command, check=True)

    # AF3 should produce exactly one non-empty model file for this task.
    # The success marker is written only after that output has been verified.
    models = list(output.glob('*/*_model.cif'))
    if len(models) != 1 or models[0].stat().st_size == 0:
        raise RuntimeError(f'Expected one top-level model, found {models}')

    # Publish the marker atomically so readers never observe a partial JSON
    # file while the controller checks whether this task succeeded.
    marker = output / 'af3_success.json'
    temporary = marker.with_suffix('.tmp')
    temporary.write_text(json.dumps({'model': str(models[0].resolve())}))
    temporary.replace(marker)


if __name__ == '__main__':
    main()
