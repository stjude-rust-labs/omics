//! Foundational representations of biological molecules.
//!
//! This crate provides typed representations of the biological molecules
//! relevant to omics analysis. It is organized into two layers that mirror the
//! biological hierarchy:
//!
//! * **[Compounds](`compound`)** — the small building blocks: individual
//!   nucleotides, along with traits that describe their biochemical properties
//!   and relationships.
//! * **[Polymers](`polymer`)** — larger molecules assembled from those building
//!   blocks: DNA and RNA sequences.
//!
//! # Compounds
//!
//! Compounds are the elemental units from which polymers are built. The
//! [`compound::Nucleotide`] trait is a marker trait that requires `Debug`,
//! `Display`, `Copy`, `Eq`, and `FromStr`. Every nucleotide reports its
//! [`compound::Kind`], which classifies it as either a
//! [`Purine`](compound::Kind::Purine) or a
//! [`Pyrimidine`](compound::Kind::Pyrimidine). Concrete nucleotide types are
//! provided by the [`polymer::dna`] and [`polymer::rna`] modules.
//!
//! Cross-polymer relationships between nucleotides are expressed through
//! additional traits:
//!
//! * [`Analogous`](compound::nucleotide::Analogous) — converts a nucleotide to
//!   its analog in a different molecular context (e.g., DNA `T` ↔ RNA `U`).
//! * [`Transcribe`](compound::nucleotide::Transcribe) — transcribes a DNA
//!   nucleotide to RNA (template-strand semantics).
//! * [`ReverseTranscribe`](compound::nucleotide::ReverseTranscribe) — reverse
//!   transcribes an RNA nucleotide back to DNA.
//!
//! # Polymers
//!
//! Polymers are sequences of compounds. Each polymer module provides a concrete
//! nucleotide type and a `Molecule` struct that wraps a `Vec` of those units:
//!
//! | Module | Compound type | Variants | Molecule |
//! |--------|--------------|----------|----------|
//! | [`polymer::dna`] | [`dna::Nucleotide`](polymer::dna::Nucleotide) | `A`, `C`, `G`, `T` | [`dna::Molecule`](polymer::dna::Molecule) |
//! | [`polymer::rna`] | [`rna::Nucleotide`](polymer::rna::Nucleotide) | `A`, `C`, `G`, `U` | [`rna::Molecule`](polymer::rna::Molecule) |
//!
//! Every `Molecule` can be parsed from a string of one-letter codes via
//! [`FromStr`](std::str::FromStr) and provides `inner()` and `into_inner()`
//! accessors. DNA and RNA molecules also expose a
//! [`gc_content()`](polymer::dna::Molecule::gc_content) method.
//!
//! # Quickstart
//!
//! ```
//! use omics_molecule::polymer::dna;
//! use omics_molecule::polymer::rna;
//!
//! // Parse a DNA sequence from a string of nucleotide codes.
//! let dna = "ACGT".parse::<dna::Molecule>()?;
//! assert_eq!(dna.inner().len(), 4);
//! assert_eq!(dna.gc_content(), 0.5);
//!
//! // Parse an RNA sequence.
//! let rna = "ACGU".parse::<rna::Molecule>()?;
//! assert_eq!(rna.inner().len(), 4);
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Individual nucleotides can be converted across molecular contexts:
//!
//! ```
//! use omics_molecule::compound::nucleotide::Analogous;
//! use omics_molecule::polymer::dna;
//! use omics_molecule::polymer::rna;
//!
//! // DNA thymine is analogous to RNA uracil.
//! let t = dna::Nucleotide::T;
//! let u: rna::Nucleotide = t.analogous();
//! assert_eq!(u, rna::Nucleotide::U);
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod compound;
pub mod polymer;
