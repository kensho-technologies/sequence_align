# Copyright 2023-present Kensho Technologies, LLC.
import multiprocessing as mp
import os
import random
from typing import Any, Callable

import psutil
import yaml


DEFAULT_GAP = "_"

CHARS = ["A", "C", "G", "T"]
PERTURB_PROB = 0.1
DELETE_PROB = 0.5
MUTATE_PROB = 0.5
INSERT_PROB = 0.5

EXPECTED_PERF_YML = os.path.join(os.path.dirname(__file__), "expected_perf.yml")
MEM_CHECK_INTERVAL = 0.01  # seconds


def _create_perturbed_seq(seq_a: list[str]) -> list[str]:
    random.seed(1234)

    seq_b = list()
    for this_entry in seq_a:
        if random.random() < PERTURB_PROB:
            # Perturb this entry -- either delete it or modify it
            if random.random() < DELETE_PROB:
                continue
            else:
                # Modify -- change the entry itself and/or insert entry/entries after it
                if random.random() < MUTATE_PROB:
                    this_entry_before = this_entry
                    this_entry = random.choice(
                        [char for char in CHARS if char != this_entry_before]
                    )
                seq_b.append(this_entry)

                while random.random() < INSERT_PROB:
                    new_entry = random.choice(CHARS)
                    seq_b.append(new_entry)
        else:
            # Leave entry as-is
            seq_b.append(this_entry)

    return seq_b


def create_seq_pair(seq_a_len: int) -> tuple[list[str], list[str]]:
    """Create a pair of sequences, where the second is a perturbed version of the first."""
    seq_a = random.choices(CHARS, k=seq_a_len)
    seq_b = _create_perturbed_seq(seq_a)
    return (seq_a, seq_b)


def get_expected_perf(key: str) -> Any:
    """Load the expected performance dictionary for the provided key."""
    with open(EXPECTED_PERF_YML, "r") as fd:
        expected_perf_full = yaml.safe_load(fd)
    return expected_perf_full[key]


def max_memory_usage(
    func: Callable[..., Any], args: tuple[Any, ...], kwargs: dict[str, Any]
) -> int:
    """Run the given function in a separate process and return the maximum memory usage in MiB."""
    max_mem = 0
    mp_proc = mp.Process(target=func, args=args, kwargs=kwargs)
    mp_proc.start()
    psutil_proc = psutil.Process(pid=mp_proc.pid)
    while True:
        try:
            max_mem = max(max_mem, psutil_proc.memory_info().rss)  # kB
            mp_proc.join(timeout=MEM_CHECK_INTERVAL)
            if mp_proc.exitcode is not None:
                break
        except psutil.NoSuchProcess:
            # Process finished since last check -- exit early
            break
    return max_mem // (1024 * 1024)
