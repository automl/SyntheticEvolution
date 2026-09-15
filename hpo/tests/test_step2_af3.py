import unittest
from unittest.mock import patch

import hpo.pipeline.trial as pipeline


class AF3Tests(unittest.TestCase):
	def test_parse_accounting(self):
		accounting = '\n'.join([
			'7_0|COMPLETED|0:0|',
			'7_1|FAILED+|1:0|',
			'',
			'malformed',
			'7_2||0:0|',
		])
		self.assertEqual(
			pipeline.parse_accounting(accounting),
			{
				'7_0': ('COMPLETED', '0:0'),
				'7_1': ('FAILED', '1:0'),
			},
		)

	def test_array_accounting(self):
		config = dict(wait_timeout_seconds=10, poll_seconds=0)
		with patch.object(
			pipeline,
			'command',
			side_effect=['', '7_0|COMPLETED|0:0|\n7_1|COMPLETED|0:0|'],
		):
			pipeline.wait_job('7', 2, config)

		with patch.object(
			pipeline,
			'command',
			return_value='7_0|FAILED|1:0|',
		):
			with self.assertRaises(RuntimeError):
				pipeline.wait_job('7', 1, config)


if __name__ == '__main__':
	unittest.main()
