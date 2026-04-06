# Copyright 2023-present Kensho Technologies, LLC.
from typing import Any
import unittest

from sequence_align.pairwise import needleman_wunsch_with_scores


DEFAULT_GAP = "_"


def match_mismatch(a: Any, b: Any) -> float:
    """Same score as default Needleman-Wunsch."""
    return 1.0 if a == b else -1.0


class TestNeedlemanWunschWithScores(unittest.TestCase):
    def test_empty(self) -> None:
        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            [], [], score_fn=match_mismatch, gap=DEFAULT_GAP
        )
        self.assertEqual(len(aligned_seq_a), 0)
        self.assertEqual(len(aligned_seq_b), 0)

    def test_one_empty(self) -> None:
        nonempty = ["A", "B", "C"]
        nonempty_aligned = ["A", "B", "C"]
        empty_aligned = [DEFAULT_GAP, DEFAULT_GAP, DEFAULT_GAP]

        with self.subTest(msg="AB"):
            aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
                nonempty, [], score_fn=match_mismatch, gap=DEFAULT_GAP
            )
            self.assertEqual(aligned_seq_a, nonempty_aligned)
            self.assertEqual(aligned_seq_b, empty_aligned)

        with self.subTest(msg="BA"):
            aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
                [], nonempty, score_fn=match_mismatch, gap=DEFAULT_GAP
            )
            self.assertEqual(aligned_seq_a, empty_aligned)
            self.assertEqual(aligned_seq_b, nonempty_aligned)

    def test_invalid_gap(self) -> None:
        for gap in ["A", "B", "C"]:
            with self.subTest(gap=gap):
                with self.assertRaises(ValueError):
                    needleman_wunsch_with_scores(
                        ["A", "B", "D"], ["A", "C", "D"], score_fn=match_mismatch, gap=gap
                    )

    def test_identity_score_matches_standard_nw(self) -> None:
        """When score_fn returns match_score/mismatch_score, results should match standard NW."""
        seq_a = ["G", "A", "T", "T", "A", "C", "A"]
        seq_b = ["G", "C", "A", "T", "G", "C", "G"]

        indel_score = -1.0

        exp_seq_a = ["G", DEFAULT_GAP, "A", "T", "T", "A", "C", "A"]
        exp_seq_b = ["G", "C", "A", DEFAULT_GAP, "T", "G", "C", "G"]

        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            seq_a, seq_b, score_fn=match_mismatch, indel_score=indel_score, gap=DEFAULT_GAP
        )
        self.assertEqual(aligned_seq_a, exp_seq_a)
        self.assertEqual(aligned_seq_b, exp_seq_b)

    def test_custom_continuous_scores(self) -> None:
        """Test with continuous (non-binary) scores to verify the matrix-based approach."""
        # Elements are numbers-as-strings; score by numeric proximity
        seq_a = ["1", "5", "9"]
        seq_b = ["2", "6", "8"]

        def numeric_proximity(a: str, b: str) -> float:
            return -abs(int(a) - int(b))

        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            seq_a, seq_b, score_fn=numeric_proximity, indel_score=-5.0, gap=DEFAULT_GAP
        )
        # Proximity: 1-2=-1, 5-6=-1, 9-8=-1 -> total=-3 (matched)
        # vs any gap arrangement which costs -5 per gap
        # So matching 1:1 is optimal
        self.assertEqual(aligned_seq_a, ["1", "5", "9"])
        self.assertEqual(aligned_seq_b, ["2", "6", "8"])

    def test_scores_prefer_gaps_over_bad_match(self) -> None:
        """When the score function returns very negative values, gaps should be preferred."""
        seq_a = ["A", "B", "C"]
        seq_b = ["X", "B", "Y"]

        def score_fn(a: str, b: str) -> float:
            if a == b:
                return 10.0
            return -100.0  # Very bad mismatch

        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            seq_a, seq_b, score_fn=score_fn, indel_score=-1.0, gap=DEFAULT_GAP
        )
        # Should match B:B and gap the rest rather than force A:X or C:Y mismatches
        self.assertEqual(aligned_seq_a, ["A", "B", DEFAULT_GAP, "C"])
        self.assertEqual(aligned_seq_b, [DEFAULT_GAP, "B", "Y", DEFAULT_GAP])

    def test_asymmetric_scores(self) -> None:
        """Test that asymmetric score functions are handled correctly."""
        seq_a = ["A", "B"]
        seq_b = ["B", "A"]

        def asymmetric_score(a: str, b: str) -> float:
            if a == "A" and b == "B":
                return 5.0  # A aligning to B is great
            if a == "B" and b == "A":
                return -5.0  # B aligning to A is terrible
            if a == b:
                return 1.0
            return -1.0

        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            seq_a, seq_b, score_fn=asymmetric_score, indel_score=-2.0, gap=DEFAULT_GAP
        )
        # A->B scores 5.0, B->A scores -5.0
        # Best: align A:B (score 5) + gap B + gap A = 5 + (-2) + (-2) = 1
        # vs: gap A + B:B (1) + gap A = -2 + 1 + -2 = -3
        # vs: A:B (5) + B:A (-5) = 0
        # So A:B + gaps is best
        self.assertEqual(aligned_seq_a, ["A", "B", DEFAULT_GAP])
        self.assertEqual(aligned_seq_b, [DEFAULT_GAP, "B", "A"])

    def test_non_string_elements(self) -> None:
        """Test that non-string sequences work (the function is generic over T)."""
        seq_a = [1, 2, 3]
        seq_b = [2, 3, 4]

        def score_fn(a: int, b: int) -> float:
            return 1.0 if a == b else -1.0

        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            seq_a, seq_b, score_fn=score_fn, indel_score=-1.0, gap=0
        )
        # Should align 2:2 and 3:3
        self.assertEqual(aligned_seq_a, [1, 2, 3, 0])
        self.assertEqual(aligned_seq_b, [0, 2, 3, 4])

    def test_words_with_custom_score(self) -> None:
        """Test word alignment with a custom similarity function."""
        seq_a = ["hello", "world", "foo"]
        seq_b = ["hallo", "welt", "baz", "foo"]

        def char_overlap_score(a: str, b: str) -> float:
            if a == b:
                return 2.0
            shared = len(set(a) & set(b))
            total = len(set(a) | set(b))
            return (2.0 * shared / total) - 1.0 if total > 0 else -1.0

        aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
            seq_a, seq_b, score_fn=char_overlap_score, indel_score=-1.0, gap=DEFAULT_GAP
        )
        # "hello" and "hallo" share {h, l, o} out of {h, e, a, l, o} -> 6/5 - 1 = 0.2
        # "world" and "welt" share {w, l} out of {w, o, r, l, d, e, t} -> 4/7 - 1 ~ -0.43
        # "foo" and "foo" -> 2.0
        # Best alignment should pair hello:hallo, world:welt, gap:baz, foo:foo
        self.assertEqual(aligned_seq_a, ["hello", "world", DEFAULT_GAP, "foo"])
        self.assertEqual(aligned_seq_b, ["hallo", "welt", "baz", "foo"])

    def test_exhaust_sequence(self) -> None:
        large = ["A", "B", "C", "D"]
        small = ["C", "D"]

        with self.subTest(msg="AB"):
            aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
                large, small, score_fn=match_mismatch, indel_score=0.0, gap=DEFAULT_GAP
            )
            self.assertEqual(aligned_seq_a, ["A", "B", "C", "D"])
            self.assertEqual(aligned_seq_b, [DEFAULT_GAP, DEFAULT_GAP, "C", "D"])

        with self.subTest(msg="BA"):
            aligned_seq_a, aligned_seq_b = needleman_wunsch_with_scores(
                small, large, score_fn=match_mismatch, indel_score=0.0, gap=DEFAULT_GAP
            )
            self.assertEqual(aligned_seq_a, [DEFAULT_GAP, DEFAULT_GAP, "C", "D"])
            self.assertEqual(aligned_seq_b, ["A", "B", "C", "D"])
