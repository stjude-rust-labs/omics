//! Biological polymers: DNA and RNA sequences.
//!
//! Each submodule provides a concrete nucleotide type and a `Molecule` wrapper
//! that represents a sequence of those nucleotides:
//!
//! * [`dna`] — Deoxyribonucleic acid sequences (`A`, `C`, `G`, `T`).
//! * [`rna`] — Ribonucleic acid sequences (`A`, `C`, `G`, `U`).
//!
//! All `Molecule` types can be parsed from one-letter code strings via
//! [`FromStr`](std::str::FromStr) and expose their inner nucleotide vectors
//! through `inner()` and `into_inner()` methods.

pub mod dna;
pub mod rna;
