#!/usr/bin/env bash
# Copyright 2023-present Kensho Technologies, LLC.
set -euxo pipefail
python -m pytest "$@"
