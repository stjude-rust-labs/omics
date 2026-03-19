# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Added in-place `move_forward(&mut self)` and `move_backward(&mut self)`
  methods to `Coordinate` that modify the position without cloning
  ([#7](https://github.com/stjude-rust-labs/omics/pull/7)).

### Revise

- Renamed consuming `move_forward(self)` and `move_backward(self)` to
  `into_move_forward` and `into_move_backward` respectively
  ([#7](https://github.com/stjude-rust-labs/omics/pull/7)).
- Fixed unnecessary `self.contig.clone()` in consuming move methods, which
  already own `self`
  ([#7](https://github.com/stjude-rust-labs/omics/pull/7)).

## 0.2.0 - 01-03-2025

### Added

- An initial set of benchmarks for some of the foundation data elements
  ([#3](https://github.com/stjude-rust-labs/omics/pull/3)).

### Revise

- This crate was completely reworked. See the pull request for more details
  ([#3](https://github.com/stjude-rust-labs/omics/pull/3)).
- Implemented `From<String>` for `Contig`
  ([#4](https://github.com/stjude-rust-labs/omics/pull/4)).
- Added the `into_start()`, `into_end()`, and `coordinate_at_offset()` methods
  for `Interval` ([#4](https://github.com/stjude-rust-labs/omics/pull/4)).
- Added the `into_equivalent_base()` and `into_equivalent_interbase()` methods
  for interbase and base intervals respectively
  ([#4](https://github.com/stjude-rust-labs/omics/pull/4)).

## 0.1.0 - 10-09-2024

### Added

- Added the initial version of `omics-coordinate` crate.
