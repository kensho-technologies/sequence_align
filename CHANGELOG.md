# Changelog

All notable changes to this package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## v0.4.0 - 2026-04-06

### Added

- `needleman_wunsch_with_scores()` function supporting custom pairwise scoring functions for alignment, enabling continuous similarity measures (e.g., spatial proximity, text edit distance) instead of binary match/mismatch.
- CHANGELOG.md and link from pyproject.toml for PyPI visibility.

### Changed

- Update Python version support to 3.10-3.14 (drop 3.9, add 3.14).

## v0.3.0 - 2025-03-05

### Changed

- Update Python version support to 3.9-3.13.
- Update GitHub Actions versions to fix wheel builds.

## v0.2.0 - 2024-08-22

### Added

- `alignment_score()` function to Python API for computing Needleman-Wunsch alignment scores on pre-aligned sequences.

## v0.1.2 - 2024-05-18

### Fixed

- Broken 0.1.1 wheels and LICENSE file.

### Changed

- PEP 639 compliance with license-file.
- Update minimum Python version to 3.8.

## v0.1.1 - 2023-04-13

### Fixed

- Bug fixes ([#10](https://github.com/kensho-technologies/sequence_align/issues/10), [#2](https://github.com/kensho-technologies/sequence_align/issues/2)).

## v0.1.0 - 2023-04-05

### Added

- Initial release with Needleman-Wunsch and Hirschberg algorithm implementations.
- Rust core with Python bindings via PyO3.
- Python 3.8-3.11 support.
