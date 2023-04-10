#!/usr/bin/env bash
# Copyright 2023-present Kensho Technologies, LLC.
python -m pytest -s --cov=src/sequence_align "$@"
