import tempfile
import unittest
from pathlib import Path

from hpo.pipeline.run import create_dataset, normalize_pairs


class InputTests(unittest.TestCase):
	def test_normalize_pairs_rejects_invalid_and_duplicates(self):
		self.assertEqual(normalize_pairs([[0, 1], [1, 0]], 2), {(0, 1)})
		for pairs in ([[0, 0]], [[-1, 2]], [[0, 4]], [[0.0, 1]], [[True, 1]]):
			with self.subTest(pairs=pairs):
				with self.assertRaises(ValueError):
					normalize_pairs(pairs, 3)

	def test_folder_inputs(self):
		with tempfile.TemporaryDirectory() as tmp:
			folder = Path(tmp)
			for index in range(2):
				(folder / f'{index}.csv').write_text(
					f'id,sequence,pairs\nrna{index},ACGU,"[[0,3]]"\n'
				)

			self.assertEqual(len(create_dataset(folder)), 2)

			(folder / 'duplicate.csv').write_text(
				'id,sequence,pairs\nrna0,ACGU,[]\n'
			)
			with self.assertRaises(ValueError):
				create_dataset(folder)


if __name__ == '__main__':
	unittest.main()
