# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Scoring infrastructure: `ScoringMatrix` trait, `GapPenalty`, `NucleotideMatrix`,
  and `Blosum62` substitution matrix.
- Sequence encoding: `EncodedSequence` with conversions from `omics-molecule`
  DNA, RNA, and protein types.
- CIGAR representation: `CigarOp`, `Cigar` with `Display` and `FromStr`.
- `Alignment` result type with score, CIGAR, and coordinate ranges.
- Scalar Needleman-Wunsch (global) alignment with affine gap penalties.
- Scalar Smith-Waterman (local) alignment with affine gap penalties.
