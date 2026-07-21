# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

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
