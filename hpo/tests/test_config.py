"""Unit tests for configuration loading and mode-specific validation."""
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import yaml
import hpo.pipeline.trial as pipeline


class TestLoadConfig(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.config_path = Path(self.temp_dir.name) / 'config.yaml'

        # Resolve test configs locally, independently of the repository layout.
        resolver = patch.object(
            pipeline, 'repo_path', side_effect=lambda value: Path(value)
        )
        resolver.start()
        self.addCleanup(resolver.stop)

        fixture_path = (
            Path(__file__).parent
            / 'fixtures'
            / 'valid_standalone_config.yaml'
        )
        self.config = yaml.safe_load(fixture_path.read_text())

    def load(self, config, mode='standalone', **kwargs):
        self.config_path.write_text(yaml.safe_dump(config))
        return pipeline.load_config(
            self.config_path, mode=mode, **kwargs
        )

    def neps_config(self):
        config = copy.deepcopy(self.config)
        config.update({
            'evaluations': 2,
            'optimizer': 'bayesian_optimization',
            'search': {
                'N': {'type': 'integer', 'lower': 2, 'upper': 128},
            },
        })
        return config

    def test_standalone_does_not_require_neps_or_controller_settings(self):
        loaded = self.load(self.config)
        self.assertEqual(loaded, self.config)

    def test_standalone_warns_about_search_and_preserves_fixed(self):
        config = self.neps_config()
        config['fixed'] = {'mutation_rate_paired': 0.2}

        with self.assertLogs(pipeline.logger, level='WARNING') as logs:
            loaded = self.load(config)

        self.assertIn("ignoring 'search'", '\n'.join(logs.output))
        self.assertEqual(loaded['fixed'], config['fixed'])

    def test_neps_requires_nonempty_search(self):
        for search in (None, {}):
            with self.subTest(search=search):
                config = self.neps_config()
                if search is None:
                    del config['search']
                else:
                    config['search'] = search

                with self.assertRaisesRegex(ValueError, 'search'):
                    self.load(config, mode='neps')

    def test_direct_neps_execution_does_not_require_submission_settings(self):
        config = self.neps_config()
        self.assertEqual(self.load(config, mode='neps'), config)

    def test_controller_submission_requires_its_settings(self):
        config = self.neps_config()

        with self.assertRaisesRegex(ValueError, 'neps_python'):
            self.load(config, mode='neps', submitting_controller=True)

        config['neps_python'] = '/example/neps/bin/python'
        with self.assertRaisesRegex(ValueError, 'controller'):
            self.load(config, mode='neps', submitting_controller=True)

        config['controller'] = {
            'count': 1,
            'partition': 'gpu-single',
            'cpus': 1,
            'memory': '20G',
            'gres': 'gpu:A100:1',
            'time': '04:00:00',
        }
        self.assertEqual(
            self.load(config, mode='neps', submitting_controller=True),
            config,
        )

    def test_controller_submission_requires_gpu_request(self):
        config = self.neps_config()
        config['neps_python'] = '/example/neps/bin/python'
        config['controller'] = {
            'count': 1,
            'partition': 'gpu-single',
            'cpus': 1,
            'memory': '20G',
            'time': '04:00:00',
        }

        with self.assertRaisesRegex(ValueError, r'controller\.gres'):
            self.load(config, mode='neps', submitting_controller=True)

        for value in ('', None, 1):
            with self.subTest(value=value):
                invalid = copy.deepcopy(config)
                invalid['controller']['gres'] = value
                with self.assertRaisesRegex(ValueError, r'controller\.gres'):
                    self.load(invalid, mode='neps', submitting_controller=True)

    def test_export_only_requires_run_location(self):
        config = {
            'run_name': 'test-run',
            'workspace_name': 'alphafold',
        }
        self.assertEqual(self.load(config, mode='export'), config)

    def test_required_value_cannot_be_missing_empty_or_null(self):
        for case in ('missing', '', None):
            with self.subTest(value=case):
                config = copy.deepcopy(self.config)
                if case == 'missing':
                    del config['dssr']
                else:
                    config['dssr'] = case

                with self.assertRaisesRegex(ValueError, 'dssr'):
                    self.load(config)

    def test_rejects_invalid_run_names(self):
        for name in ('', '../other-run', 'run/name', 'run name', None):
            with self.subTest(name=name):
                config = {**self.config, 'run_name': name}
                with self.assertRaisesRegex(ValueError, 'run_name'):
                    self.load(config)

    def test_rejects_invalid_count(self):
        for value in (0, -1, 1.5, True, '2'):
            with self.subTest(value=value):
                config = self.neps_config()
                config['controller']['count'] = value
                config['neps_python'] = '/example/neps/bin/python'
                with self.assertRaisesRegex(ValueError, r'controller\.count'):
                    self.load(config, mode='neps', submitting_controller=True)

    def test_rejects_overlapping_fixed_and_search_parameters(self):
        config = self.neps_config()
        config['fixed'] = {'N': 20}

        with self.assertRaisesRegex(ValueError, 'both fixed and search'):
            self.load(config, mode='neps')

    def test_rejects_yaml_that_is_not_a_mapping(self):
        for value in (None, ['not', 'a', 'mapping']):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'YAML mapping'):
                    self.load(value)


if __name__ == '__main__':
    unittest.main()