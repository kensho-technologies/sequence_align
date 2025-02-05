# Copyright 2023-present Kensho Technologies, LLC.
import time
from typing import Any
import unittest

from sequence_align.pairwise import needleman_wunsch

from .utils import create_seq_pair, get_expected_perf, max_memory_usage


# Fix these so that we run with the same scores, even if defaults change
MATCH_SCORE = 1.0
MISMATCH_SCORE = -1.0
INDEL_SCORE = -1.0
DEFAULT_GAP = "_"

RUNTIME_SEQ_A_LEN = 5_000
RUNTIME_TRIALS = 9

MEMORY_SEQ_A_LEN = 10_000
MEMORY_TRIALS = 5


class TestNeedlemanWunsch(unittest.TestCase):
    # Needed for mypy to not complain
    expected_perf: dict[str, Any] = dict()

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()

        cls.expected_perf = get_expected_perf("needleman_wunsch")

    def test_runtime(self) -> None:
        seq_a, seq_b = create_seq_pair(RUNTIME_SEQ_A_LEN)

        runtimes = list()
        for _ in range(RUNTIME_TRIALS):
            start_t = time.perf_counter()
            needleman_wunsch(
                seq_a,
                seq_b,
                match_score=MATCH_SCORE,
                mismatch_score=MISMATCH_SCORE,
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
            msg=f"""Expected runtime to be within {tolerance * 100.}% of {exp_median:.3f}s.
Got {median_runtime:.3f}s ({sign}{100. * diff_pct}%) instead.

Consider adjusting the median number and/or tolerance if this change in performance is expected.""",
        )

    def test_memory(self) -> None:
        seq_a, seq_b = create_seq_pair(MEMORY_SEQ_A_LEN)

        max_mems = list()
        for _ in range(MEMORY_TRIALS):
            max_mem = max_memory_usage(
                needleman_wunsch,
                (seq_a, seq_b),
                {
                    "match_score": MATCH_SCORE,
                    "mismatch_score": MISMATCH_SCORE,
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
            msg=f"""Expected memory to be within {tolerance * 100.}% of {exp_median:.3f}MiB.
Got {median_max_mem:.3f}MiB ({sign}{100. * diff_pct}%) instead.

Consider adjusting the median number and/or tolerance if this change in performance is expected.""",
        )
