# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Added the optional `alignment` component (`--features alignment`), which re-exports
  `omics-alignment` as `omics::alignment` and automatically enables the `coordinate`
  component ([#18](https://github.com/stjude-rust-labs/omics/pull/18)).
* Added the `omics::alignment::algorithm` API under the existing `alignment`
  feature, including the `global` and `local` affine-gap aligners
  ([#19](https://github.com/stjude-rust-labs/omics/pull/19)).

### Changed

* Raised the minimum supported Rust version to `1.87`
  ([#22](https://github.com/stjude-rust-labs/omics/pull/22)).
* Moved to the Rust 2024 edition and updated dependencies
  ([#23](https://github.com/stjude-rust-labs/omics/pull/23)).

### Crate Updates

- `omics-coordinate`: bumped to v0.5.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.5.0))
- `omics-molecule`: bumped to v0.3.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-molecule-v0.3.0))
- `omics-variation`: bumped to v0.4.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-variation-v0.4.0))

## 0.4.0 - 03-19-2026

### Crate Updates

- `omics-coordinate`: bumped to v0.4.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.4.0))

## 0.3.0 - 03-19-2026

### Crate Updates

- `omics-coordinate`: bumped to v0.3.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.3.0))
- `omics-molecule`: bumped to v0.2.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-molecule-v0.2.0))
- `omics-variation`: bumped to v0.3.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-variation-v0.3.0))

## 0.2.0 - 01-03-2025

### Revise

- Promotes the `position-u64` feature of `omics-coordinate` to the top-level
  `omics` crate ([#4](https://github.com/stjude-rust-labs/omics/pull/4)).

### Crate Updates

- `omics-coordinate`: bumped to v0.2.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.2.0))
- `omics-variation`: bumped to v0.2.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-variation-v0.2.0))

## 0.1.0 - 10-09-2024

### Crate Updates

- `omics-core`: introduced at v0.1.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-core-v0.1.0))
- `omics-coordinate`: introduced at v0.1.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-coordinate-v0.1.0))
- `omics-molecule`: introduced at v0.1.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-molecule-v0.1.0))
- `omics-variation`: introduced at v0.1.0
  ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-variation-v0.1.0))
