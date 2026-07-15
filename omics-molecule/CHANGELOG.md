# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Added a generic `Sequence<N>` allele type, the generic sibling of the
  concrete `dna::Molecule` / `rna::Molecule`, with `FromStr` and `TryFrom<&str>`
  parsing (a lone `.` denotes the empty allele)
  ([#15](https://github.com/stjude-rust-labs/omics/pull/15)).
* Added a `Complement` trait giving the Watson-Crick complement of the DNA and
  RNA nucleotides, and a `Sequence::reverse_complement` method built on it
  ([#16](https://github.com/stjude-rust-labs/omics/pull/16)).

### Changed

* Raised the minimum supported Rust version to `1.81`
  ([#16](https://github.com/stjude-rust-labs/omics/pull/16)).

## 0.2.0 - 03-19-2026

### Changed

* Used `thiserror` for better error messages and improved existing error
  messages ([#5](https://github.com/stjude-rust-labs/omics/pull/5)).

## 0.1.0 - 10-09-2024

### Added

* Added the initial version of `omics-molecule` crate.
