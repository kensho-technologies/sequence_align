# Copyright 2023-present Kensho Technologies, LLC.
import time
from typing import Any
import unittest

from sequence_align.pairwise import needleman_wunsch_with_scores

from .utils import create_seq_pair, get_expected_perf, max_memory_usage


INDEL_SCORE = -1.0
DEFAULT_GAP = "_"

RUNTIME_SEQ_A_LEN = 5_000
RUNTIME_TRIALS = 9

MEMORY_SEQ_A_LEN = 10_000
MEMORY_TRIALS = 5


def char_overlap_score(a: str, b: str) -> float:
    """Score based on character set overlap — a continuous similarity measure."""
    if a == b:
        return 2.0
    shared = len(set(a) & set(b))
    total = len(set(a) | set(b))
    return (2.0 * shared / total) - 1.0 if total > 0 else -1.0


class TestNeedlemanWunschWithScores(unittest.TestCase):
    # Needed for mypy to not complain
    expected_perf: dict[str, Any] = dict()

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()

        cls.expected_perf = get_expected_perf("needleman_wunsch_with_scores")

    def test_runtime(self) -> None:
        seq_a, seq_b = create_seq_pair(RUNTIME_SEQ_A_LEN)

        runtimes = list()
        for _ in range(RUNTIME_TRIALS):
            start_t = time.perf_counter()
            needleman_wunsch_with_scores(
                seq_a,
                seq_b,
                score_fn=char_overlap_score,
                indel_score=INDEL_SCORE,
                gap=DEFAULT_GAP,
            )
            end_t = time.perf_counter()
            runtimes.append(end_t - start_t)

        median_runtime = sorted(runtimes)[len(runtimes) // 2]
        exp_median = self.expected_perf["runtime"]["median"]
        diff = median_runtime - exp_median
        sign = "+" if diff > 0 else "-"
        diff_pct = abs(diff / exp_median)

        tolerance = self.expected_perf["runtime"]["tolerance"]
        self.assertLessEqual(
            abs(diff_pct),
            tolerance,
            msg=f"""Expected runtime to be within {tolerance * 100.0}% of {exp_median:.3f}s.
Got {median_runtime:.3f}s ({sign}{100.0 * diff_pct}%) instead.

Consider adjusting the median number and/or tolerance if this change in performance is expected.""",
        )

    def test_memory(self) -> None:
        seq_a, seq_b = create_seq_pair(MEMORY_SEQ_A_LEN)

        max_mems = list()
        for _ in range(MEMORY_TRIALS):
            max_mem = max_memory_usage(
                needleman_wunsch_with_scores,
                (seq_a, seq_b),
                {
                    "score_fn": char_overlap_score,
                    "indel_score": INDEL_SCORE,
                    "gap": DEFAULT_GAP,
                },
            )
            max_mems.append(max_mem)

        median_max_mem = sorted(max_mems)[len(max_mems) // 2]
        exp_median = self.expected_perf["memory"]["median"]
        diff = median_max_mem - exp_median
        sign = "+" if diff > 0 else "-"
        diff_pct = abs(diff / exp_median)

        tolerance = self.expected_perf["memory"]["tolerance"]
        self.assertLessEqual(
            abs(diff_pct),
            tolerance,
            msg=f"""Expected memory to be within {tolerance * 100.0}% of {exp_median:.3f}MiB.
Got {median_max_mem:.3f}MiB ({sign}{100.0 * diff_pct}%) instead.

Consider adjusting the median number and/or tolerance if this change in performance is expected.""",
        )
