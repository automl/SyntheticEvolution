import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

from hpo.pipeline.run import write_json


def test_write_json_handles_simultaneous_workers(tmp_path, monkeypatch):
    target = tmp_path / "run_metadata.json"
    metadata = {
        "fingerprint": "same-run",
        "config": {"run_name": "concurrent-test"},
    }

    ready_to_replace = Barrier(2)
    original_replace = Path.replace

    def synchronized_replace(source, destination):
        if Path(destination) == target:
            # Both temporary files must be written before either rename.
            ready_to_replace.wait(timeout=5)
        return original_replace(source, destination)

    monkeypatch.setattr(Path, "replace", synchronized_replace)

    with ThreadPoolExecutor(max_workers=2) as workers:
        futures = [
            workers.submit(write_json, target, metadata)
            for _ in range(2)
        ]

        # Propagate worker exceptions: neither write may fail.
        for future in futures:
            future.result(timeout=10)

    assert json.loads(target.read_text()) == metadata

    # Both writes must also clean up their temporary files.
    assert set(tmp_path.iterdir()) == {target}