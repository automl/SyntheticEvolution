"""Per-RNA archives: publishing, validation and extraction."""
import io
import tarfile

import pytest

from hpo.pipeline import archive


def content(tmp_path, model='model'):
    directory = tmp_path / 'content'
    (directory / 'seed-42_sample-0').mkdir(parents=True)
    (directory / 'rna_00000_model.cif').write_text(model)
    (directory / 'rna_00000_data.json').write_text('{"msa": "synthetic"}')
    (directory / 'seed-42_sample-0' / 'model.cif').write_text('sample')
    return directory


def test_publish_and_extract_round_trip(tmp_path):
    target = tmp_path / 'trial' / 'rna_00000.tar.gz'
    target.parent.mkdir()
    archive.publish(content(tmp_path), target)
    # Only the archive is published; the temporary file is gone.
    assert [path.name for path in target.parent.iterdir()] == ['rna_00000.tar.gz']
    assert archive.model_name(target) == 'rna_00000_model.cif'
    assert archive.read_data_json(target) == b'{"msa": "synthetic"}'
    archive.extract(target, tmp_path / 'out')
    assert (tmp_path / 'out' / 'seed-42_sample-0' / 'model.cif').read_text() == 'sample'


def test_republishing_replaces_the_archive(tmp_path):
    target = tmp_path / 'rna_00000.tar.gz'
    directory = content(tmp_path)
    archive.publish(directory, target)
    (directory / 'dssr.json').write_text('{}')
    archive.publish(directory, target)
    with tarfile.open(target) as tar:
        assert 'dssr.json' in tar.getnames()


def test_failed_publish_leaves_no_partial_archive(tmp_path, monkeypatch):
    target = tmp_path / 'trial' / 'rna_00000.tar.gz'
    target.parent.mkdir()
    def broken_add(self, *args, **kwargs):
        raise OSError('disk full')
    monkeypatch.setattr(tarfile.TarFile, 'add', broken_add)
    with pytest.raises(OSError):
        archive.publish(content(tmp_path), target)
    assert list(target.parent.iterdir()) == []


@pytest.mark.parametrize('case', ['missing', 'truncated', 'empty_model', 'no_model'])
def test_unusable_archive_has_no_model(tmp_path, case):
    target = tmp_path / 'rna_00000.tar.gz'
    if case == 'truncated':
        archive.publish(content(tmp_path), target)
        target.write_bytes(target.read_bytes()[:40])
    elif case == 'empty_model':
        archive.publish(content(tmp_path, model=''), target)
    elif case == 'no_model':
        directory = content(tmp_path)
        (directory / 'rna_00000_model.cif').unlink()
        archive.publish(directory, target)
    assert archive.model_name(target) is None


def test_sample_models_are_not_the_top_level_model(tmp_path):
    directory = content(tmp_path)
    (directory / 'rna_00000_model.cif').unlink()
    (directory / 'seed-42_sample-0' / 'rna_00000_model.cif').write_text('sample')
    target = tmp_path / 'rna_00000.tar.gz'
    archive.publish(directory, target)
    assert archive.model_name(target) is None


def test_extract_rejects_paths_outside_destination(tmp_path):
    target = tmp_path / 'evil.tar.gz'
    with tarfile.open(target, 'w:gz') as tar:
        data = b'x'
        info = tarfile.TarInfo('../outside.txt')
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with pytest.raises(ValueError, match='Unexpected archive member'):
        archive.extract(target, tmp_path / 'out')
    assert not (tmp_path / 'outside.txt').exists()


def test_remove_stale_temporaries(tmp_path):
    (tmp_path / '.rna_00000.tar.gz.node.1.tmp').write_text('partial')
    (tmp_path / '.trial_result.json.abc.tmp').write_text('partial')
    (tmp_path / 'rna_00001.tar.gz').write_text('kept')
    archive.remove_stale_temporaries(tmp_path)
    assert [path.name for path in tmp_path.iterdir()] == ['rna_00001.tar.gz']
