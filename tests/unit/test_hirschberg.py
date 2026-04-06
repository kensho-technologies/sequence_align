# Copyright 2023-present Kensho Technologies, LLC.
from copy import deepcopy
import unittest

from sequence_align.pairwise import hirschberg


DEFAULT_GAP = "-1"


class TestHirschberg(unittest.TestCase):
    def test_empty(self) -> None:
        aligned_seq_a, aligned_seq_b = hirschberg([], [])
        self.assertEqual(len(aligned_seq_a), 0)
        self.assertEqual(len(aligned_seq_b), 0)

    def test_normal(self) -> None:
        # See https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm#Example
        seq_a = ["A", "G", "T", "A", "C", "G", "C", "A"]
        seq_b = ["T", "A", "T", "G", "C"]
        exp_seq_a = ["A", "G", "T", "A", "C", "G", "C", "A"]
        exp_seq_b = [DEFAULT_GAP, DEFAULT_GAP, "T", "A", "T", "G", "C", DEFAULT_GAP]

        aligned_seq_a, aligned_seq_b = hirschberg(
            seq_a,
            seq_b,
            match_score=2.0,
            mismatch_score=-1.0,
            indel_score=-2.0,
            gap=DEFAULT_GAP,
        )
        self.assertEqual(aligned_seq_a, exp_seq_a)
        self.assertEqual(aligned_seq_b, exp_seq_b)

    # Needleman-Wunsch tests that DIFFER in results due to pivot selection
    def test_words(self) -> None:
        # Test that words also work (i.e., verify that the Python logic to map into and out of
        # integer indices works).
        seq_a = ["Hello", "world", "I", "am", "Sam"]
        seq_b = ["What's", "up", "world", "I", "am", "called", "Sam", "too"]

        exp_seq_a_aligned = [
            DEFAULT_GAP,
            DEFAULT_GAP,
            "Hello",
            "world",
            "I",
            "am",
            DEFAULT_GAP,
            "Sam",
            DEFAULT_GAP,
        ]
        exp_seq_b_aligned = [
            "What's",
            "up",
            DEFAULT_GAP,
            "world",
            "I",
            "am",
            "called",
            "Sam",
            "too",
        ]
        seq_a_aligned, seq_b_aligned = hirschberg(
            seq_a,
            seq_b,
            indel_score=0.0,  # Don't punish gaps
            gap=DEFAULT_GAP,
        )
        self.assertEqual(seq_a_aligned, exp_seq_a_aligned)
        self.assertEqual(seq_b_aligned, exp_seq_b_aligned)

    def test_encourage_gaps(self) -> None:
        seq_a = ["A", "B", "D"]
        seq_b = ["A", "C", "D"]

        # Same whether we prefer matches over mismatches or vice versa if gaps are most preferable.
        # Prefer creating a gap in sequence B over A first (see Notes in lib.rs)
        exp_seq_a = ["A", DEFAULT_GAP, DEFAULT_GAP, "B", "D", DEFAULT_GAP]
        exp_seq_b = [DEFAULT_GAP, "A", "C", DEFAULT_GAP, DEFAULT_GAP, "D"]
        with self.subTest(msg="prefer matches to mismatches"):
            aligned_seq_a, aligned_seq_b = hirschberg(
                seq_a,
                seq_b,
                match_score=1.0,
                mismatch_score=-1.0,
                indel_score=100.0,
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, exp_seq_a)
            self.assertEqual(aligned_seq_b, exp_seq_b)

        with self.subTest(msg="prefer mismatches to matches"):
            aligned_seq_a, aligned_seq_b = hirschberg(
                seq_a,
                seq_b,
                match_score=-1.0,
                mismatch_score=1.0,
                indel_score=100.0,
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, exp_seq_a)
            self.assertEqual(aligned_seq_b, exp_seq_b)

    def test_encourage_mismatches(self) -> None:
        seq_a = ["A", "B", "A"]
        seq_b = ["A", "A", "A"]

        with self.subTest(msg="prefer matches to gaps"):
            exp_seq_a = ["A", "B", "A"]
            exp_seq_b = ["A", "A", "A"]
            aligned_seq_a, aligned_seq_b = hirschberg(
                seq_a,
                seq_b,
                match_score=1.0,
                mismatch_score=100.0,
                indel_score=-1.0,
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, exp_seq_a)
            self.assertEqual(aligned_seq_b, exp_seq_b)

        with self.subTest(msg="prefer gaps to matches"):
            # Prefer creating a gap in sequence B over A first (see Notes in lib.rs)
            exp_seq_a = ["A", DEFAULT_GAP, DEFAULT_GAP, "B", "A"]
            exp_seq_b = [DEFAULT_GAP, "A", "A", "A", DEFAULT_GAP]

            aligned_seq_a, aligned_seq_b = hirschberg(
                seq_a,
                seq_b,
                match_score=-10.0,
                mismatch_score=100.0,
                indel_score=-1.0,
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, exp_seq_a)
            self.assertEqual(aligned_seq_b, exp_seq_b)

    # Needleman-Wunsch tests that should NOT differ

    def test_one_empty(self) -> None:
        nonempty = ["A", "B", "C"]
        nonempty_aligned = ["A", "B", "C"]

        empty_aligned = [DEFAULT_GAP, DEFAULT_GAP, DEFAULT_GAP]

        # Should work in either configuration
        with self.subTest(msg="AB"):
            aligned_seq_a, aligned_seq_b = hirschberg(nonempty, [], gap=DEFAULT_GAP)
            self.assertEqual(aligned_seq_a, nonempty_aligned)
            self.assertEqual(aligned_seq_b, empty_aligned)

        with self.subTest(msg="BA"):
            aligned_seq_a, aligned_seq_b = hirschberg([], nonempty, gap=DEFAULT_GAP)
            self.assertEqual(aligned_seq_a, empty_aligned)
            self.assertEqual(aligned_seq_b, nonempty_aligned)

    def test_invalid_gap(self) -> None:
        # Should fail if the gap is in either or both sequences
        for gap in ["A", "B", "C"]:
            with self.subTest(gap=gap):
                with self.assertRaises(ValueError):
                    hirschberg(["A", "B", "D"], ["A", "C", "D"], gap=gap)

    def test_exhaust_sequence(self) -> None:
        # Test that exhausting one sequence early still yields correct results
        large = ["A", "B", "C", "D"]
        small = ["C", "D"]
        large_aligned = deepcopy(large)
        small_aligned = [DEFAULT_GAP, DEFAULT_GAP, "C", "D"]

        # Should work in either configuration
        with self.subTest(msg="AB"):
            aligned_seq_a, aligned_seq_b = hirschberg(
                large,
                small,
                indel_score=0.0,  # Don't punish gaps
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, large_aligned)
            self.assertEqual(aligned_seq_b, small_aligned)

        with self.subTest(msg="BA"):
            aligned_seq_a, aligned_seq_b = hirschberg(
                small,
                large,
                indel_score=0.0,  # Don't punish gaps
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, small_aligned)
            self.assertEqual(aligned_seq_b, large_aligned)

    def test_encourage_matches(self) -> None:
        seq_a = ["A", "B", "D"]
        seq_b = ["A", "C", "D"]

        with self.subTest(msg="prefer mismatches to gaps"):
            exp_seq_a = ["A", "B", "D"]
            exp_seq_b = ["A", "C", "D"]
            aligned_seq_a, aligned_seq_b = hirschberg(
                seq_a,
                seq_b,
                match_score=100.0,
                mismatch_score=-1.0,
                indel_score=-10.0,
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, exp_seq_a)
            self.assertEqual(aligned_seq_b, exp_seq_b)

        with self.subTest(msg="prefer gaps to mismatches"):
            # Prefer creating a gap in sequence B over A first (see Notes in lib.rs)
            exp_seq_a = ["A", "B", DEFAULT_GAP, "D"]
            exp_seq_b = ["A", DEFAULT_GAP, "C", "D"]

            aligned_seq_a, aligned_seq_b = hirschberg(
                seq_a,
                seq_b,
                match_score=100.0,
                mismatch_score=-10.0,
                indel_score=-1.0,
                gap=DEFAULT_GAP,
            )
            self.assertEqual(aligned_seq_a, exp_seq_a)
            self.assertEqual(aligned_seq_b, exp_seq_b)

    def test_normal_needleman_wunsch(self) -> None:
        # See https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#/media/File:Needleman-Wunsch_pairwise_sequence_alignment.png  # noqa: E501
        # Use Needleman-Wunsch default scores (match=1, mismatch=-1, indel=-1)
        seq_a = ["G", "A", "T", "T", "A", "C", "A"]
        seq_b = ["G", "C", "A", "T", "G", "C", "G"]
        exp_seq_a = ["G", DEFAULT_GAP, "A", "T", "T", "A", "C", "A"]
        exp_seq_b = ["G", "C", "A", DEFAULT_GAP, "T", "G", "C", "G"]

        aligned_seq_a, aligned_seq_b = hirschberg(
            seq_a,
            seq_b,
            match_score=1.0,
            mismatch_score=-1.0,
            indel_score=-1.0,
            gap=DEFAULT_GAP,
        )
        self.assertEqual(aligned_seq_a, exp_seq_a)
        self.assertEqual(aligned_seq_b, exp_seq_b)
