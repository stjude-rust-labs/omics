# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.1.0 - 11-15-2026

### Added

* Added exact-parity SIMD dispatch for byte-slice global and local affine-gap
  alignment with Apple Silicon NEON and Linux `x86_64` AVX2 wavefront kernels.
  Dispatch logs the selected backend (`debug`, once per process) and each
  call's lane width and scalar fallbacks (`trace`) through `tracing`.
* Added the initial version of the `omics-alignment` crate, including a complete
  checked CIGAR model with `OperationKind`, `Operation`, and `Cigar` types supporting
  `FromStr`, `Display`, and aggregate length validation
  ([#18](https://github.com/stjude-rust-labs/omics/pull/18)).
* Added `Alignment`, which eagerly validates that a starting reference coordinate,
  a starting query coordinate, and a `Cigar` form a traversal that stays within
  representable bounds; `Alignment::steps` yields one `Step` per CIGAR operation
  losslessly, and `Step::is_aligned` supports filtering to aligned operations
  (`M`, `=`, `X`) without discarding operation boundaries
  ([#18](https://github.com/stjude-rust-labs/omics/pull/18)).
* Added generic global and local affine-gap alignment with deterministic canonical
  CIGAR output using `=`, `X`, `I`, and `D`, plus returned half-open local input
  ranges that place each result on the original sequences
  ([#19](https://github.com/stjude-rust-labs/omics/pull/19)).

### Changed

* Raised the minimum supported Rust version to `1.87`
  ([#22](https://github.com/stjude-rust-labs/omics/pull/22)).
* Moved to the Rust 2024 edition and updated dependencies
  ([#23](https://github.com/stjude-rust-labs/omics/pull/23)).
