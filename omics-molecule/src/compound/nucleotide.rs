//! Nucleotides.

pub mod relation;

pub use relation::Relation;

use crate::compound::Kind;

/// A marker trait that denotes a type of nucleotide.
pub trait Nucleotide:
    std::fmt::Debug + std::fmt::Display + Copy + Eq + PartialEq + std::str::FromStr
{
    /// Gets the [`Kind`] type for a given [`Nucleotide`].
    fn kind(&self) -> Kind;
}

/// A trait that provides methods to convert a [`Nucleotide`] to an analogous
/// [`Nucleotide`] in a different molecular context.
pub trait Analogous<T: Nucleotide>
where
    Self: Nucleotide,
{
    /// Converts a [`Nucleotide`] to an analogous [`Nucleotide`] in a different
    /// molecular context.
    fn analogous(&self) -> T;
}

/// A trait that provides methods to transcribe a [`Nucleotide`] to the
/// respective [`Nucleotide`] in a different molecular context (generally from
/// DNA to RNA).
pub trait Transcribe<T: Nucleotide>
where
    Self: Nucleotide,
{
    /// Transcribes a [`Nucleotide`] to the respective [`Nucleotide`] in a
    /// different molecular context (generally from DNA to RNA).
    fn transcribe(&self) -> T;
}

/// A trait that provides methods to reverse transcribe a [`Nucleotide`] to the
/// respective [`Nucleotide`] in a different molecular context (generally from
/// RNA to DNA).
pub trait ReverseTranscribe<T: Nucleotide>
where
    Self: Nucleotide,
{
    /// Reverse transcribes a [`Nucleotide`] to the respective [`Nucleotide`] in
    /// a different molecular context (generally from RNA to DNA).
    fn reverse_transcribe(&self) -> T;
}
