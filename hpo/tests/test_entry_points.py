import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENTRY_POINTS = (
    'hpo.run_standalone_trial',
    'hpo.run_neps_hpo',
    'hpo.submit_standalone_trial',
    'hpo.submit_neps_hpo',
)


class EntryPointTests(unittest.TestCase):
    def test_package_entry_points_show_help(self):
        """Package entry points can start without submitting cluster jobs."""
        for module in ENTRY_POINTS:
            with self.subTest(module=module):
                result = subprocess.run(
                    [sys.executable, '-m', module, '--help'],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)


if __name__ == '__main__':
    unittest.main()
