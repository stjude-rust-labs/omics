//! Small molecular building blocks: nucleotides.
//!
//! This module provides the trait definitions and classification types for the
//! elemental units from which biological polymers are assembled:
//!
//! * [`Nucleotide`] — a marker trait for nucleotide types, with [`Kind`]
//!   distinguishing purines from pyrimidines.
//!
//! Concrete types that implement these traits live in the
//! [`polymer`](crate::polymer) submodules ([`dna`](crate::polymer::dna) and
//! [`rna`](crate::polymer::rna)).

mod kind;
pub mod nucleotide;

pub use kind::Kind;
pub use nucleotide::Nucleotide;
