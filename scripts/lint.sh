#!/usr/bin/env bash
# Copyright 2021-present Kensho Technologies, LLC.

# Treat undefined variables and non-zero exits in pipes as errors.
set -uo pipefail

# Ensure that the "**" glob operator is applied recursively.
# Make globs that do not match return null values.
shopt -s globstar nullglob

# Break on first error.
set -e

# Parse input arguments.
any_run_only_set=0
run_fast_linters=0  # copyright line check, ruff check, ruff format, cargo fmt
run_mypy=0
run_cargo_clippy=0
fix=0
for i in "$@"; do
    case $i in
        --run-only-fast-linters )
            if [ "$any_run_only_set" -eq 1 ]; then
                echo "Multiple run-only options set, this is not supported.";
                exit 1;
            fi
            any_run_only_set=1
            run_fast_linters=1
            shift;;

        --run-only-mypy )
            if [ "$any_run_only_set" -eq 1 ]; then
                echo "Multiple run-only options set, this is not supported.";
                exit 1;
            fi
            any_run_only_set=1
            run_mypy=1
            shift;;

        --run-only-cargo-clippy )
            if [ "$any_run_only_set" -eq 1 ]; then
                echo "Multiple run-only options set, this is not supported.";
                exit 1;
            fi
            any_run_only_set=1
            run_cargo_clippy=1
            shift;;

        --fix )
            fix=1
            echo "running in FIX mode"
            shift;;

        *)
            echo "Unknown option: $i";
            exit 1;;
    esac
done

if [ "$any_run_only_set" -eq 0 ]; then
    run_fast_linters=1
    run_mypy=1
    run_cargo_clippy=1
fi

# Make sure the current working directory for this script is the root directory.
cd "$(git -C "$(dirname "${0}")" rev-parse --show-toplevel )"

# Continue on error to allow ignoring certain linters.
# Errors are manually aggregated at the end.
set +e

if [ "$run_fast_linters" -eq 1 ]; then
    echo -e '*** Running copyright line check... ***\n'
    ./scripts/copyright_line_check.sh
    copyright_line_check_exit_code=$?
    echo -e "\n*** End of copyright line check run; exit: $copyright_line_check_exit_code ***\n"

    echo -e '*** Running ruff check... ***\n'
    if [ "$fix" -eq 1 ]; then
        ruff check --fix .
        ruff_check_exit_code=$?
    else
        ruff check .
        ruff_check_exit_code=$?
    fi
    echo -e "\n*** End of ruff check run; exit: $ruff_check_exit_code ***\n"

    echo -e '*** Running ruff format... ***\n'
    if [ "$fix" -eq 1 ]; then
        ruff format .
        ruff_format_exit_code=$?
    else
        ruff format --check --diff .
        ruff_format_exit_code=$?
    fi
    echo -e "\n*** End of ruff format run; exit: $ruff_format_exit_code ***\n"

    echo -e '\n*** Running cargo fmt...\n'
    if [ "$fix" -eq 1 ]; then
        cargo fmt -v --all --manifest-path=./rust/Cargo.toml
        cargo_fmt_exit_code=$?
    else
        cargo fmt -v --all --manifest-path=./rust/Cargo.toml --check
        cargo_fmt_exit_code=$?
    fi
    echo -e "\n*** End of cargo fmt run, exit: $cargo_fmt_exit_code ***\n"
fi

if [ "$run_mypy" -eq 1 ]; then
    echo -e '*** Running mypy... ***\n'
    mypy .
    mypy_exit_code=$?
    echo -e "\n*** End of mypy run, exit: $mypy_exit_code ***\n"
fi

if [ "$run_cargo_clippy" -eq 1 ]; then
    # Warn about pedantic stuff; deny all other defaults
    echo -e '\n*** Running cargo clippy...\n'
    cargo_clippy_flags="--manifest-path=./rust/Cargo.toml --release -- -W clippy::pedantic -D clippy::all"
    cargo_clippy_exit_code=0
    if [ "$fix" -eq 1 ]; then
        cargo clippy --fix --allow-dirty --allow-staged $cargo_clippy_flags
        cargo_clippy_exit_code=$?
    fi
    if [ "$cargo_clippy_exit_code" -eq 0 ]; then
        # Even if running in fix mode and the fixes "worked", may still report errors upon check
        cargo clippy --no-deps $cargo_clippy_flags
        cargo_clippy_exit_code=$?
    fi
    echo -e "\n*** End of cargo clippy run, exit: $cargo_clippy_exit_code ***\n"
fi

if  [[
        (
            ("$run_fast_linters" == 1) && (
                ("$copyright_line_check_exit_code" != "0") ||
                ("$ruff_check_exit_code" != "0") ||
                ("$ruff_format_exit_code" != "0") ||
                ("$cargo_fmt_exit_code" != "0")
            )
        ) || (
            ("$run_mypy" == 1) && ("$mypy_exit_code" != "0")
        ) || (
            ("$run_cargo_clippy" == 1) && ("$cargo_clippy_exit_code" != "0")
        )
    ]]; then
    echo -e "\n*** Lint failed. ***\n"

    if [ "$run_fast_linters" -eq 1 ]; then
        echo -e "copyright line check exit: $copyright_line_check_exit_code"
        echo -e "ruff check exit: $ruff_check_exit_code"
        echo -e "ruff format exit: $ruff_format_exit_code"
        echo -e "cargo fmt exit: $cargo_fmt_exit_code"
    fi
    if [ "$run_mypy" -eq 1 ]; then
        echo -e "mypy exit: $mypy_exit_code"
    fi
    if [ "$run_cargo_clippy" -eq 1 ]; then
        echo -e "cargo clippy exit: $cargo_clippy_exit_code"
    fi

    exit 1
fi

echo -e "\n*** Lint successful. ***\n"
