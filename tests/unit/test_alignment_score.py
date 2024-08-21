# Copyright 2023-present Kensho Technologies, LLC.
import unittest

from sequence_align.pairwise import alignment_score


# Try something non-default
DEFAULT_GAP = "?"


class TestAlignmentScore(unittest.TestCase):
    def test_empty(self) -> None:
        """Score of two empty sequences should always be zero."""
        self.assertEqual(alignment_score([], []), 0)

    def test_unequal(self) -> None:
        """Should fail with sequences of different length."""
        with self.assertRaises(ValueError):
            alignment_score(["A", "B", "C"], ["D", "E"])

    def test_normal(self) -> None:
        """Score of two nonempty sequences should match and in both directions."""
        # Should be one insertion, one match, and one mismatch based on the scores
        seq_a = ["A", "B", "C"]
        seq_b = [DEFAULT_GAP, "B", "D"]
        match_score = 1.234
        mismatch_score = -1.234
        indel_score = -1.012

        expected_score = match_score + mismatch_score + indel_score

        for ab_order in [False, True]:
            seq_a_proc = seq_a if ab_order else seq_b
            seq_b_proc = seq_b if ab_order else seq_a
            with self.subTest(ab_order=ab_order):
                self.assertEqual(
                    alignment_score(
                        seq_a_proc,
                        seq_b_proc,
                        match_score=match_score,
                        mismatch_score=mismatch_score,
                        indel_score=indel_score,
                        gap=DEFAULT_GAP,
                    ),
                    expected_score,
                )
