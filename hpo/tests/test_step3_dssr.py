import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import hpo.pipeline.trial as pipeline
from hpo.pipeline.dssr import parse_pairs


FIXTURE = json.loads(
	(Path(__file__).parent / 'fixtures' / 'dssr_example.json').read_text()
)
SEQ = 'AGUAGUAGUAGUAGA'


class DSSRTests(unittest.TestCase):
	def test_dssr_json_parser(self):
		predicted = parse_pairs(FIXTURE, SEQ)
		self.assertEqual(
			predicted,
			[(0, 14), (1, 13), (2, 12), (3, 11), (4, 10), (5, 9)],
		)

	def test_mapping_rejects_missing_and_mismatch(self):
		with self.assertRaises(ValueError):
			parse_pairs(FIXTURE, 'A' * len(SEQ))
		broken = dict(FIXTURE, nts=FIXTURE['nts'][:-1])
		with self.assertRaises(ValueError):
			parse_pairs(broken, SEQ)

	def test_resume_and_failure_no_partial_loss(self):
		with tempfile.TemporaryDirectory() as tmp:
			directory = Path(tmp)
			config = dict(fixed={}, dssr='fake', dssr_timeout_seconds=10)
			rows = [dict(
				id='one',
				sequence=SEQ,
				pairs=[[0, 14], [1, 13], [7, 12]],
			)]
			pipeline.write_json(
				directory / 'trial_metadata.json',
				dict(fingerprint='abc', parameters={}),
			)
			pipeline.write_json(
				directory / 'af3_job_id.json',
				dict(job_id='7'),
			)
			model = directory / 'rna_00000/af3/model.cif'
			model.parent.mkdir(parents=True)
			model.write_text('mock model')
			pipeline.write_json(
				model.parent / 'af3_success.json',
				dict(model=str(model)),
			)

			with (
				patch.object(pipeline, 'wait_job') as wait,
				patch.object(
					pipeline.dssr,
					'run_dssr',
					return_value=parse_pairs(FIXTURE, SEQ),
				),
			):
				loss = pipeline.run_pipeline(config, rows, directory, 'abc')
				self.assertAlmostEqual(loss, 5 / 9)
				self.assertAlmostEqual(
					pipeline.run_pipeline(config, rows, directory, 'abc'),
					loss,
				)
				self.assertEqual(wait.call_count, 1)

			with self.assertRaises(ValueError):
				pipeline.run_pipeline(config, rows, directory, 'abc', {'N': 3})

			(directory / 'trial_result.json').unlink()
			with patch.object(pipeline, 'wait_job', side_effect=RuntimeError('failed')):
				with self.assertRaises(RuntimeError):
					pipeline.run_pipeline(config, rows, directory, 'abc')
			self.assertFalse((directory / 'trial_result.json').exists())


if __name__ == '__main__':
	unittest.main()
