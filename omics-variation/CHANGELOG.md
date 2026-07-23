# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Added the `MNV`, insertion, deletion, and delins small-variant kinds
  alongside the existing SNV, backed by a shared `Alteration` core and a new
  `Kind` classification. The top-level `Variant` now parses and dispatches to
  each kind and exposes a parsimony `normalize` method. Each kind exposes a
  `try_new` taking string-like inputs and a `TryFrom<(Coordinate, Alteration)>`
  conversion for pre-built alleles
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* Added `reference_interval` and `alternate_interval` APIs for variant kinds
  whose reference and alternate alleles occupy different coordinate systems or
  spans. `SNV` and `MNV` variants use a single `interval` method because the
  reference and alternate spans are identical
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* Added Criterion benchmarks for qualified variant parsing, rejected
  coordinate-system qualifiers, display, normalization, and interval queries
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* Added a structural-variant tier under the `structural` module, modeling
  oriented `Breakend`s, single novel `Adjacency` junctions carrying an optional
  non-templated insertion, and a `StructuralVariant` event whose `Kind`
  (deletion, insertion, tandem duplication, inversion, translocation, breakend,
  or complex) is derived from breakend geometry rather than stored. The tier
  parses and serializes a compact crate-local string format and ships Criterion
  benchmarks ([#16](https://github.com/stjude-rust-labs/omics/pull/16)).
* Added a copy-number tier with strandless half-open regions, typed absolute
  `Count` observations, nonzero `Ploidy` baselines, canonical
  `contig:start-end(i):copies` serialization, baseline-relative `Change`
  classification, strict base-2 and base-10 logarithmic conversions, and
  top-level `CopyNumber*` aliases.

### Changed

* Moved the small-variant tier into a `small` module. The previous paths remain
  available through re-exports, including `omics_variation::variant`
  ([#16](https://github.com/stjude-rust-labs/omics/pull/16)).
* Gave the structural `Kind::Translocation` a `Join` payload alongside its
  `Locality`, so a co-linear fusion is distinguished from a fold-back
  (inverted) one, including across contigs
  ([#16](https://github.com/stjude-rust-labs/omics/pull/16)).
* Raised the minimum supported Rust version to `1.81`
  ([#16](https://github.com/stjude-rust-labs/omics/pull/16)).
* **Breaking:** renamed the `Variant::SingleNucleotideVariation` enum variant to
  `Variant::Snv` for consistency with the new kind variants (`Mnv`,
  `Insertion`, `Deletion`, `Delins`)
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* **Breaking:** reworked the top-level `Error`: removed `Error::ParseError` in
  favor of `InvalidFormat`, `Coordinate`, `ReferenceSequence`,
  `AlternateSequence`, `Alteration`, `Kind`, and `NormalizeOverflow` variants
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* **Breaking:** replaced `snv::Error::Relation` with `snv::Error::Alteration`
  now that SNVs are built on the shared `Alteration` core
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* **Breaking:** required serialized variant positions to include `(b)` for base
  coordinates or `(i)` for interbase coordinates. `SNV`, `MNV`, deletion, and
  delins records accept only `(b)`, while insertion records accept only `(i)`
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).

### Crate Updates

- `omics-coordinate`: bumped to v0.4.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.4.0))

## 0.3.0 - 03-19-2026

### Changed

* Used `thiserror` for better error messages and improved existing error
  messages ([#5](https://github.com/stjude-rust-labs/omics/pull/5)).

### Crate Updates

- `omics-coordinate`: bumped to v0.3.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.3.0))
- `omics-molecule`: bumped to v0.2.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-molecule-v0.2.0))

## 0.2.0 - 01-03-2025

### Crate Updates

- `omics-coordinate`: bumped to v0.2.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.2.0))

## 0.1.0 - 10-09-2024

### Added

- Added the initial version of `omics-variation` crate.
