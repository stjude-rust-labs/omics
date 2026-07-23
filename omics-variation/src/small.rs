//! The small-variant tier.
//!
//! Every small variant is a reference allele and an alternate allele anchored
//! at a coordinate. [`Alteration`] is the shared core carrying the two alleles
//! and classifying them into a [`Kind`]. The typed modules in this namespace
//! add the coordinate type and kind-specific invariants.
//!
//! [`Alteration`] stores allele sequence only. It does not store a coordinate,
//! and it does not remember whether the alleles came from a string or from
//! typed construction. This keeps classification deterministic: the same
//! `reference` and `alternate` always produce the same [`Kind`].
//!
//! ```
//! use omics_molecule::polymer::dna::Nucleotide;
//! use omics_variation::variant::Alteration;
//! use omics_variation::variant::Kind;
//!
//! let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "G".parse()?)?;
//! assert_eq!(alteration.kind(), Kind::Delins);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Use [`Alteration`] when code needs to reason about alleles independent of a
//! locus. Use [`snv::Variant`], [`mnv::Variant`], [`insertion::Variant`],
//! [`deletion::Variant`], or [`delins::Variant`] when code needs a fully
//! anchored variant with span accessors.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::position::Number;
use omics_coordinate::system::Base;
use omics_coordinate::system::Interbase;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;
use thiserror::Error;

pub mod deletion;
pub mod delins;
pub mod insertion;
pub mod mnv;
pub mod snv;

/// The classification of a small variant, derived from its alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// A single-base substitution.
    Snv,
    /// A same-length substitution of two or more bases.
    Mnv,
    /// Bases inserted where none existed (empty reference).
    Insertion,
    /// Bases removed (empty alternate).
    Deletion,
    /// A combined deletion and insertion of differing lengths.
    Delins,
}

/// Builds a base interval from a start coordinate and allele length.
pub(crate) fn base_interval(start: &Coordinate<Base>, len: usize) -> Option<Interval<Base>> {
    let span = Number::try_from(len.checked_sub(1)?).ok()?;
    let end = start.clone().into_move_forward(span)?;
    Interval::try_from((start.clone(), end)).ok()
}

/// Builds a zero-width interbase interval at a boundary coordinate.
pub(crate) fn interbase_interval(coordinate: &Coordinate<Interbase>) -> Interval<Interbase> {
    // SAFETY: equal endpoints share context and form a stationary span.
    Interval::try_from((coordinate.clone(), coordinate.clone())).unwrap()
}

/// Builds a zero-width interbase interval immediately before a base coordinate.
pub(crate) fn interbase_interval_before_base(coordinate: &Coordinate<Base>) -> Interval<Interbase> {
    // SAFETY: every valid base coordinate has a preceding interbase boundary.
    let boundary = coordinate.clone().nudge_backward().unwrap();
    interbase_interval(&boundary)
}

/// An error related to an [`Alteration`].
#[derive(Error, Debug)]
pub enum Error {
    /// Both alleles were empty, which is not a variant.
    #[error("cannot create an alteration with two empty alleles")]
    BothEmpty,

    /// The reference and alternate alleles were identical.
    #[error("reference and alternate alleles are both `{0}`")]
    Identical(String),
}

/// An error raised when constructing a typed variant.
#[derive(Error, Debug)]
pub enum KindError {
    /// The alteration classified as a different kind than the target type.
    #[error("wrong alteration kind; expected {expected:?}, found {found:?}")]
    WrongKind {
        /// The kind the target type requires.
        expected: Kind,
        /// The kind the alteration actually is.
        found: Kind,
    },

    /// The variant's span cannot be represented (allele length or end
    /// coordinate overflows the position bounds).
    #[error("variant span overflows the coordinate position bounds")]
    SpanOverflow,

    /// The raw alleles could not form a valid alteration.
    #[error(transparent)]
    Alteration(#[from] Error),

    /// The coordinate could not be parsed.
    #[error(transparent)]
    Coordinate(#[from] omics_coordinate::coordinate::Error),

    /// An allele could not be parsed.
    #[error(transparent)]
    Sequence(#[from] omics_molecule::sequence::ParseError),
}

/// The reference and alternate alleles of a small variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alteration<N: Nucleotide> {
    /// The reference allele (empty for an insertion).
    reference: Sequence<N>,

    /// The alternate allele (empty for a deletion).
    alternate: Sequence<N>,
}

impl<N: Nucleotide> Alteration<N> {
    /// Attempts to create a new [`Alteration`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    ///
    /// let alteration = Alteration::<Nucleotide>::try_new("A".parse()?, "C".parse()?)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(reference: Sequence<N>, alternate: Sequence<N>) -> Result<Self, Error> {
        if reference.is_empty() && alternate.is_empty() {
            return Err(Error::BothEmpty);
        }

        if reference == alternate {
            return Err(Error::Identical(reference.to_string()));
        }

        Ok(Self {
            reference,
            alternate,
        })
    }

    /// Gets the reference allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    ///
    /// let alteration = Alteration::<Nucleotide>::try_new("A".parse()?, "C".parse()?)?;
    /// assert_eq!(alteration.reference().to_string(), "A");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &Sequence<N> {
        &self.reference
    }

    /// Gets the alternate allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    ///
    /// let alteration = Alteration::<Nucleotide>::try_new("A".parse()?, "C".parse()?)?;
    /// assert_eq!(alteration.alternate().to_string(), "C");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &Sequence<N> {
        &self.alternate
    }

    /// Classifies this [`Alteration`] into a [`Kind`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    /// use omics_variation::variant::Kind;
    ///
    /// let alteration = Alteration::<Nucleotide>::try_new("A".parse()?, "C".parse()?)?;
    /// assert_eq!(alteration.kind(), Kind::Snv);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn kind(&self) -> Kind {
        match (self.reference.len(), self.alternate.len()) {
            (0, _) => Kind::Insertion,
            (_, 0) => Kind::Deletion,
            (1, 1) => Kind::Snv,
            (x, y) if x == y => Kind::Mnv,
            _ => Kind::Delins,
        }
    }

    /// Trims shared prefix and suffix bases, returning the trimmed prefix
    /// length and the trimmed [`Alteration`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    ///
    /// let alteration = Alteration::<Nucleotide>::try_new("AAT".parse()?, "AAG".parse()?)?;
    /// let (prefix, trimmed) = alteration.trimmed();
    /// assert_eq!(prefix, 2);
    /// assert_eq!(trimmed.reference().to_string(), "T");
    /// assert_eq!(trimmed.alternate().to_string(), "G");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn trimmed(&self) -> (usize, Alteration<N>) {
        let prefix = self.reference.shared_prefix_len(&self.alternate);
        let suffix = self.reference.shared_suffix_len(&self.alternate);

        let trim = |seq: &Sequence<N>| {
            let inner = seq.inner();
            Sequence::new(inner[prefix..inner.len() - suffix].to_vec())
        };

        let reference = trim(&self.reference);
        let alternate = trim(&self.alternate);

        // SAFETY: the original alteration was valid (not both-empty, not
        // identical). Trimming the common prefix/suffix from two non-identical
        // sequences leaves at least one differing base, so the result cannot be
        // both-empty or identical.
        (
            prefix,
            Alteration {
                reference,
                alternate,
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use omics_molecule::polymer::dna;

    use super::*;

    fn alt(reference: &str, alternate: &str) -> Alteration<dna::Nucleotide> {
        Alteration::try_new(reference.parse().unwrap(), alternate.parse().unwrap()).unwrap()
    }

    #[test]
    fn it_classifies_every_kind() {
        assert_eq!(alt("A", "C").kind(), Kind::Snv);
        assert_eq!(alt("AT", "GC").kind(), Kind::Mnv);
        assert_eq!(alt(".", "AT").kind(), Kind::Insertion);
        assert_eq!(alt("AT", ".").kind(), Kind::Deletion);
        assert_eq!(alt("AT", "G").kind(), Kind::Delins);
    }

    #[test]
    fn it_rejects_both_empty() {
        let err =
            Alteration::<dna::Nucleotide>::try_new(".".parse().unwrap(), ".".parse().unwrap())
                .unwrap_err();
        assert!(matches!(err, Error::BothEmpty));
    }

    #[test]
    fn it_rejects_identical_alleles() {
        let err =
            Alteration::<dna::Nucleotide>::try_new("AT".parse().unwrap(), "AT".parse().unwrap())
                .unwrap_err();
        assert!(matches!(err, Error::Identical(_)));
    }

    #[test]
    fn it_trims_shared_flanks() {
        // AAT>AAG  ->  prefix 2, T>G
        let (prefix, trimmed) = alt("AAT", "AAG").trimmed();
        assert_eq!(prefix, 2);
        assert_eq!(trimmed.reference().to_string(), "T");
        assert_eq!(trimmed.alternate().to_string(), "G");

        // ATG>ACG  ->  prefix 1, T>C (suffix G trimmed)
        let (prefix, trimmed) = alt("ATG", "ACG").trimmed();
        assert_eq!(prefix, 1);
        assert_eq!(trimmed.reference().to_string(), "T");
        assert_eq!(trimmed.alternate().to_string(), "C");

        // Nothing shared: unchanged.
        let (prefix, trimmed) = alt("A", "C").trimmed();
        assert_eq!(prefix, 0);
        assert_eq!(trimmed.reference().to_string(), "A");
        assert_eq!(trimmed.alternate().to_string(), "C");
    }
}
