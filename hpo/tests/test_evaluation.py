import unittest

from hpo.pipeline.scoring import aggregate, score_pairs


class EvaluationTests(unittest.TestCase):
	def test_pair_scoring_and_empty_cases(self):
		predicted = [(0, 14), (1, 13), (2, 12), (3, 11), (4, 10), (5, 9)]
		score = score_pairs({(0, 14), (1, 13), (7, 12)}, predicted)
		self.assertEqual((score['tp'], score['fp'], score['fn']), (2, 4, 1))
		self.assertAlmostEqual(score['loss'], 5 / 9)
		self.assertEqual(score_pairs([], [])['loss'], 0)
		self.assertEqual(score_pairs([(0, 1)], [])['loss'], 1)

	def test_aggregate_rejects_empty_scores(self):
		with self.assertRaises(ValueError):
			aggregate([])


if __name__ == '__main__':
	unittest.main()
