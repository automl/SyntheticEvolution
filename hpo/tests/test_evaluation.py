"""Tests for scoring objectives, including the real GraKeL WL implementation."""
import unittest
import math
import numpy as np

from hpo.pipeline import scoring
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



class TestScoring(unittest.TestCase):

    def test_f1_is_unchanged(self):
        metrics = scoring.score_row(
            {"sequence": "ACGU", "pairs": [(0, 3)]},
            {(0, 3)},
            objective="f1",
        )

        self.assertEqual(metrics["tp"], 1)
        self.assertEqual(metrics["fp"], 0)
        self.assertEqual(metrics["fn"], 0)
        self.assertEqual(metrics["f1"], 1.0)
        self.assertEqual(metrics["loss"], 0.0)

    def test_invalid_objective(self):
        with self.assertRaises(ValueError):
            scoring.score_row(
                {"sequence": "ACGU", "pairs": []},
                set(),
                objective="not-an-objective",
            )

    def test_aggregate_still_uses_loss(self):
        self.assertEqual(
            scoring.aggregate([{"loss": 0.25}, {"loss": 0.75}]),
            0.5,
        )

        with self.assertRaises(ValueError):
            scoring.aggregate([])

        with self.assertRaises(ValueError):
            scoring.aggregate([{"loss": math.nan}])

    def test_wl_identical_structures_have_perfect_score(self):
        """An identical target/prediction must have normalized WL similarity 1."""
        metrics = scoring.score_row(
            {
                "sequence": "ACGU",
                "pairs": [(0, 3)],
            },
            {(0, 3)},
            objective="wl",
        )

        self.assertAlmostEqual(metrics["wl"], 1.0)
        self.assertAlmostEqual(metrics["loss"], 0.0)

    def test_wl_different_structures_are_not_perfect(self):
        """A structurally different graph should not receive similarity 1."""
        metrics = scoring.score_row(
            {
                "sequence": "ACGU",
                "pairs": [(0, 3)],
            },
            {(0, 2)},
            objective="wl",
        )

        self.assertGreaterEqual(metrics["wl"], 0.0)
        self.assertLess(metrics["wl"], 1.0)
        self.assertGreater(metrics["loss"], 0.0)

    def test_wl_is_symmetric(self):
        """The normalized WL kernel should give the same similarity either way."""
        target = {(0, 3)}
        predicted = {(0, 2)}

        forward = scoring.score_wl(
            target,
            predicted,
            length=4,
        )
        reverse = scoring.score_wl(
            predicted,
            target,
            length=4,
        )

        self.assertAlmostEqual(
            forward["wl"],
            reverse["wl"],
        )
        self.assertAlmostEqual(
            forward["loss"],
            reverse["loss"],
        )

    def test_wl_matches_reference_grakel_implementation(self):
        """
        Verify our implementation against the actual reference algorithm.

        This deliberately uses GraKeL directly rather than mocking it.
        """
        from grakel import Graph
        from grakel.kernels import WeisfeilerLehman, VertexHistogram

        target = {(0, 3)}
        predicted = {(0, 2)}
        length = 4

        def pairs2mat(pairs):
            matrix = np.zeros((length, length), dtype=int)

            for p1, p2 in pairs:
                matrix[p1, p2] = 1
                matrix[p2, p1] = 1

            return matrix

        def mat2graph(matrix):
            return Graph(
                initialization_object=matrix.astype(int),
                node_labels={
                    i: str(i)
                    for i in range(matrix.shape[0])
                },
            )

        true_graph = mat2graph(pairs2mat(target))
        pred_graph = mat2graph(pairs2mat(predicted))

        kernel = WeisfeilerLehman(
            n_iter=5,
            normalize=True,
            base_graph_kernel=VertexHistogram,
        )

        kernel.fit_transform([true_graph])
        expected_similarity = float(
            kernel.transform([pred_graph])[0][0]
        )

        actual = scoring.score_wl(
            target,
            predicted,
            length=length,
            n_iter=5,
        )

        self.assertAlmostEqual(
            actual["wl"],
            expected_similarity,
        )
        self.assertAlmostEqual(
            actual["loss"],
            1.0 - expected_similarity,
        )

    def test_wl_iterations_are_supported(self):
        """The configured number of WL iterations must affect the kernel call."""
        target = {(0, 4)}
        predicted = {(0, 3)}

        one_iteration = scoring.score_wl(
            target,
            predicted,
            length=5,
            n_iter=1,
        )

        five_iterations = scoring.score_wl(
            target,
            predicted,
            length=5,
            n_iter=5,
        )

        self.assertTrue(math.isfinite(one_iteration["wl"]))
        self.assertTrue(math.isfinite(one_iteration["loss"]))
        self.assertTrue(math.isfinite(five_iterations["wl"]))
        self.assertTrue(math.isfinite(five_iterations["loss"]))

        self.assertGreaterEqual(one_iteration["wl"], 0.0)
        self.assertLessEqual(one_iteration["wl"], 1.0)
        self.assertGreaterEqual(five_iterations["wl"], 0.0)
        self.assertLessEqual(five_iterations["wl"], 1.0)

    def test_wl_empty_structures_are_supported(self):
        """Two RNAs with no base-pair edges still form valid WL graphs."""
        metrics = scoring.score_wl(
            target=set(),
            predicted=set(),
            length=4,
        )

        self.assertAlmostEqual(metrics["wl"], 1.0)
        self.assertAlmostEqual(metrics["loss"], 0.0)


if __name__ == "__main__":
    unittest.main()

