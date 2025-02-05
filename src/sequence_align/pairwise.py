# Copyright 2023-present Kensho Technologies, LLC.
from typing import Sequence

from sequence_align import _sequence_align  # type: ignore


_GAP_VAL = -1


def _entry2idx(
    seq_a: Sequence[str],
    seq_b: Sequence[str],
    gap: str,
    allow_gap: bool = False,
) -> tuple[dict[int, str], list[int], list[int]]:
    symbols = set(seq_a).union(set(seq_b))
    if not allow_gap and gap in symbols:
        raise ValueError(f'Gap entry "{gap}" found in seq_a and/or seq_b; must not exist in either')

    symbols_without_gap = symbols - {gap}
    idx2symbol: dict[int, str] = {
        _GAP_VAL: gap,
        **{idx: symbol for idx, symbol in enumerate(sorted(symbols_without_gap))},
    }
    symbol2idx = {symbol: idx for idx, symbol in idx2symbol.items()}

    seq_a_indices = [symbol2idx[symbol] for symbol in seq_a]
    seq_b_indices = [symbol2idx[symbol] for symbol in seq_b]

    return (idx2symbol, seq_a_indices, seq_b_indices)


def _idx2entry(
    idx2symbol: dict[int, str],
    seq_a_indices_aligned: list[int],
    seq_b_indices_aligned: list[int],
    gap: str,
) -> tuple[list[str], list[str]]:
    seq_a_aligned = [gap if idx == _GAP_VAL else idx2symbol[idx] for idx in seq_a_indices_aligned]
    seq_b_aligned = [gap if idx == _GAP_VAL else idx2symbol[idx] for idx in seq_b_indices_aligned]
    return (seq_a_aligned, seq_b_aligned)


def needleman_wunsch(
    seq_a: Sequence[str],
    seq_b: Sequence[str],
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    indel_score: float = -1.0,
    gap: str = "-",
) -> tuple[list[str], list[str]]:
    """Compute an optimal global pairwise alignment using the Needleman-Wunsch algorithm.

    Args:
        seq_a: First sequence in pair to align.
        seq_b: Second sequence in pair to align.
        match_score: Score to apply for transitions where the sequences match each other at a given
            index. Defaults to 1.
        mismatch_score: Score to apply for transitions where the sequences do _not_ match each other
            at a given index. Defaults to -1.
        indel_score: Score to apply for insertion/deletion transitions where one sequence advances
            without the other advancing (thus inserting a gap). Defaults to -1.
        gap: Value to use for marking a gap in one sequence in the final output. Cannot be present
            in `seq_a` and/or `seq_b`. Defaults to "-".

    Returns:
        Sequences A and B, respectively, aligned to each other with gaps represented by `gap`.

    Raises:
        ValueError: If `gap` is found in `seq_a` and/or `seq_b`.

    Note:
        Unlike other implementations, this only considers a **single backpointer** when backtracing
        the optimal pairwise alignment, rather than potentially two or three backpointers for each
        cell if the scores are equal. Rather, this will prioritize "up" transitions (*i.e.*, gap in
        `seq_b`) over "left" transitions (*i.e.*, gap in `seq_a`), which in turn is prioritized over
        "diagonal" transitions (*i.e.*, no gap). This is a somewhat arbitrary distinction, but is
        consistent and leads to a simpler implementation that is both faster and uses less memory.

        This takes O(mn) time and O(mn) space complexity, where m and n are the lengths of the two
        sequences, respectively.

        See https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm for more information.
    """
    # First, map the sequences to integers
    idx2symbol, seq_a_indices, seq_b_indices = _entry2idx(seq_a, seq_b, gap)

    # Now, run alignment in Rust
    seq_a_indices_aligned, seq_b_indices_aligned = _sequence_align.needleman_wunsch(
        seq_a_indices,
        seq_b_indices,
        match_score=match_score,
        mismatch_score=mismatch_score,
        indel_score=indel_score,
        gap_val=_GAP_VAL,
    )

    # Finally, map back and return
    return _idx2entry(idx2symbol, seq_a_indices_aligned, seq_b_indices_aligned, gap)


def hirschberg(
    seq_a: Sequence[str],
    seq_b: Sequence[str],
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    indel_score: float = -1.0,
    gap: str = "-",
) -> tuple[list[str], list[str]]:
    """Compute an optimal global pairwise alignment using the Hirschberg algorithm.

    Args:
        seq_a: First sequence in pair to align.
        seq_b: Second sequence in pair to align.
        match_score: Score to apply for transitions where the sequences match each other at a given
            index. Defaults to 1.
        mismatch_score: Score to apply for transitions where the sequences do _not_ match each other
            at a given index. Defaults to -1.
        indel_score: Score to apply for insertion/deletion transitions where one sequence advances
            without the other advancing (thus inserting a gap). Defaults to -1.
        gap: Value to use for marking a gap in one sequence in the final output. Cannot be present
            in `seq_a` and/or `seq_b`. Defaults to "-".

    Returns:
        Sequences A and B, respectively, aligned to each other with gaps represented by `gap`.

    Raises:
        ValueError: If `gap` is found in `seq_a` and/or `seq_b`.

    Note:
        This is a modification of Needleman-Wunsch that reduces space complexity from
        O(mn) to O(min(m, n))) and returns the corresponding aligned sequences, with any gaps
        represented by `gap`. **Please note** that this will not necessarily generate the exact same
        alignments as those returned by the `needleman_wunsch` implementation in this package!

        Unlike other implementations, this only considers a **single backpointer** when backtracing
        the optimal pairwise alignment, rather than potentially two or three backpointers for each
        cell if the scores are equal. Rather, this will prioritize "up" transitions (*i.e.*, gap in
        `seq_b`) over "left" transitions (*i.e.*, gap in `seq_a`), which in turn is prioritized over
        "diagonal" transitions (*i.e.*, no gap). This is a somewhat arbitrary distinction, but is
        consistent and leads to a simpler implementation that is both faster and uses less memory.

        This takes O(mn) time and O(min(m, n)) space complexity, where m and n are the lengths of
        the two sequences, respectively.

        See https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm for more information.
    """
    # First, map the sequences to integers
    idx2symbol, seq_a_indices, seq_b_indices = _entry2idx(seq_a, seq_b, gap)

    # Now, run alignment in Rust
    seq_a_indices_aligned, seq_b_indices_aligned = _sequence_align.hirschberg(
        seq_a_indices,
        seq_b_indices,
        match_score=match_score,
        mismatch_score=mismatch_score,
        indel_score=indel_score,
        gap_val=_GAP_VAL,
    )

    # Finally, map back and return
    return _idx2entry(idx2symbol, seq_a_indices_aligned, seq_b_indices_aligned, gap)


def alignment_score(
    aligned_seq_a: Sequence[str],
    aligned_seq_b: Sequence[str],
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    indel_score: float = -1.0,
    gap: str = "-",
) -> float:
    """Compute the alignment score for the pair of sequences.

    Args:
        aligned_seq_a: First aligned sequence.
        aligned_seq_b: Second aligned sequence.
        match_score: Score to apply for transitions where the sequences match each other at a given
            index. Defaults to 1.
        mismatch_score: Score to apply for transitions where the sequences do _not_ match each other
            at a given index. Defaults to -1.
        indel_score: Score to apply for insertion/deletion transitions where one sequence advances
            without the other advancing (thus inserting a gap). Defaults to -1.
        gap: Value to use for marking gaps in the aligned sequences. Defaults to "-".

    Returns:
        Needleman-Wunsch alignment score representing the sum of match, mismatch and
        insertion/deletion transition scores for the provided alignment.

    Note:
        See `NWScore()` function at
        https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm for more information.
    """
    # First, map the sequences to integers
    _, aligned_seq_a_indices, aligned_seq_b_indices = _entry2idx(
        aligned_seq_a, aligned_seq_b, gap, allow_gap=True
    )

    # Now, get the score
    return float(
        _sequence_align.alignment_score(
            aligned_seq_a_indices,
            aligned_seq_b_indices,
            match_score=match_score,
            mismatch_score=mismatch_score,
            indel_score=indel_score,
            gap_val=_GAP_VAL,
        )
    )
