# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Added the `MNV`, insertion, deletion, and delins small-variant kinds
  alongside the existing SNV, backed by a shared `Alteration` core and a new
  `Kind` classification. The top-level `Variant` now parses and dispatches to
  each kind and exposes a parsimony `normalize` method.

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
