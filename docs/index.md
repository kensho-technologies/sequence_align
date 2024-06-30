## Introduction

Efficient implementations of Needleman-Wunsch and other sequence alignment algorithms written in Rust with Python bindings via PyO3.

## Installation

```bash
pip install sequence_align
```

## Quickstart

```python
seq_a = ["G", "A", "T", "T", "A", "C", "A"]
seq_b = ["G", "C", "A", "T", "G", "C", "G"]

aligned_seq_a, aligned_seq_b = needleman_wunsch(
    seq_a,
    seq_b,
    match_score=1.0,
    mismatch_score=-1.0,
    indel_score=-1.0,
    gap="_",
)
```
