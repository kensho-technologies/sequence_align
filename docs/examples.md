## Needleman-Wunsch

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

# Expects ["G", "_", "A", "T", "T", "A", "C", "A"]
print(aligned_seq_a)

# Expects ["G", "C", "A", "_", "T", "G", "C", "G"]
print(aligned_seq_b)
```


## Hirschberg

```python
seq_a = ["A", "G", "T", "A", "C", "G", "C", "A"]
seq_b = ["T", "A", "T", "G", "C"]

aligned_seq_a, aligned_seq_b = hirschberg(
    seq_a,
    seq_b,
    match_score=2.0,
    mismatch_score=-1.0,
    indel_score=-2.0,
    gap="_",
)

# Expects ["A", "G", "T", "A", "C", "G", "C", "A"]
print(aligned_seq_a)

# Expects ["_", "_", "T", "A", "T", "G", "C", "_"]
print(aligned_seq_b)
```
