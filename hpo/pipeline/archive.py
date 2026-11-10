"""Per-RNA archives: one workspace file for everything kept from one RNA.

An archive holds the complete AF3 output (flattened, without AF3's
``<name>/`` folder), ``shs_effective_parameters.json`` and, once the pairs were extracted,
``dssr.json`` (RNA) or ``contacts.json`` (protein). Its existence marks a completed AF3 prediction.
"""

import os
import shutil
import socket
import tarfile
from pathlib import Path
from typing import Optional

SUFFIX = '.tar.gz'
MODEL_SUFFIX = '_model.cif'
DATA_SUFFIX = '_data.json'


def path_for(trial_directory: Path, task_name: str) -> Path:
    return trial_directory / f'{task_name}{SUFFIX}'


def publish(content: Path, archive: Path) -> None:
    """Pack ``content`` and publish it atomically as ``archive``.

    The archive is written next to the target and renamed, so readers never
    see a partial file, also when the controller is killed while packing.
    """
    temporary = archive.with_name(f'.{archive.name}.{socket.gethostname()}.{os.getpid()}.tmp')
    try:
        with tarfile.open(temporary, 'w:gz') as tar:
            for path in sorted(content.rglob('*')):
                if path.is_file():
                    tar.add(path, arcname=path.relative_to(content).as_posix(), recursive=False)
        temporary.replace(archive)
    finally:
        temporary.unlink(missing_ok=True)


def _top_level_member(tar: tarfile.TarFile, suffix: str) -> Optional[tarfile.TarInfo]:
    members = [
        member for member in tar.getmembers()
        if member.isfile() and '/' not in member.name and member.name.endswith(suffix)
    ]
    return members[0] if len(members) == 1 else None


def model_name(archive: Path) -> Optional[str]:
    """Name of the top-level model in a readable archive, else None.

    A missing, truncated or otherwise broken archive counts as no prediction.
    """
    if not archive.is_file():
        return None
    try:
        with tarfile.open(archive, 'r:gz') as tar:
            member = _top_level_member(tar, MODEL_SUFFIX)
            return member.name if member is not None and member.size > 0 else None
    except (OSError, EOFError, tarfile.TarError):
        return None


def extract(archive: Path, destination: Path) -> None:
    """Extract an archive written by ``publish`` into ``destination``."""
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar.getmembers():
            target = (destination / member.name).resolve()
            # Only regular files below destination (Python 3.9 has no extraction filters).
            if not member.isfile() or destination not in target.parents:
                raise ValueError(f'Unexpected archive member {member.name!r} in {archive}')
            target.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as source, target.open('wb') as stream:
                shutil.copyfileobj(source, stream)


def read_data_json(archive: Path) -> bytes:
    """AF3's ``<name>_data.json`` (the AF3 input with the synthetic MSA)."""
    with tarfile.open(archive, 'r:gz') as tar:
        member = _top_level_member(tar, DATA_SUFFIX)
        if member is None:
            raise FileNotFoundError(f'No {DATA_SUFFIX} in {archive}')
        with tar.extractfile(member) as stream:
            return stream.read()


def remove_stale_temporaries(trial_directory: Path) -> None:
    """Remove leftovers of publishes interrupted by a killed controller."""
    for path in trial_directory.glob('.*.tmp'):
        path.unlink(missing_ok=True)
