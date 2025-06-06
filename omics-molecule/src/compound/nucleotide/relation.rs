//! Relationship between an expected nucleotide and an actual nucleotide.
//!
//! [`Relation`] is the foundational struct that represents variation (or lack
//! thereof) at a single nucleotide level. The following relations are
//! supported:
//!
//! * [`Relation::Identical`]—the reference nucleotide and the alternate
//!   nucleotide are identical.
//! * [`Relation::Substitution`]—the nucleobase within a nucleotide was
//!   substituted for another nucleobase.
//! * [`Relation::Insertion`]—a nucleotide was inserted into a position where
//!   there was previously no nucleotide.
//! * [`Relation::Deletion`]—a nucleotide that previously existed was removed.
//!
//! Importantly, both nucleotides cannot be missing. Attempting to
//! create such a [`Relation`] will result in a [`Error::Empty`].

pub mod substitution;

use omics_core::MISSING_NUCLEOTIDE;
use omics_core::VARIANT_SEPARATOR;
pub use substitution::Substitution;
use thiserror::Error;

use crate::compound::Nucleotide;

/// An error related to a [`Relation`].
#[derive(Error, Debug)]
pub enum Error<N: Nucleotide> {
    /// Attempted to create a relation with no nucleotides.
    #[error("cannot create a relation with no nucleotides")]
    Empty,

    /// A parse error.
    #[error("parse error: {0}")]
    ParseError(String),

    /// A substitution error.
    #[error(transparent)]
    Substitution(#[from] substitution::Error<N>),
}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T, N> = std::result::Result<T, Error<N>>;

/// A relation between an expected [`Nucleotide`] and the existing
/// [`Nucleotide`].
#[derive(Debug, Eq, PartialEq)]
pub enum Relation<N: Nucleotide> {
    /// Two nucleotides that are identical.
    Identical(N),

    /// The nucleobase within a nucleotide was substituted for another
    /// nucleobase.
    Substitution(Substitution<N>),

    /// A nucleotide now exists where none did previously.
    Insertion(N),

    /// A nucleotide that previously existed now does not.
    Deletion(N),
}

impl<N: Nucleotide> Relation<N> {
    /// Attempts to create a new [`Relation`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::Relation;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let relation = Relation::<Nucleotide>::try_new(Some(Nucleotide::A), Some(Nucleotide::C))?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(reference: Option<N>, alternate: Option<N>) -> Result<Self, N> {
        match (reference, alternate) {
            (None, None) => Err(Error::Empty),
            (None, Some(alternate)) => Ok(Self::Insertion(alternate)),
            (Some(reference), None) => Ok(Self::Deletion(reference)),
            (Some(reference), Some(alternate)) => {
                if reference == alternate {
                    Ok(Self::Identical(reference))
                } else {
                    Ok(Self::Substitution(
                        Substitution::try_new(reference, alternate).map_err(Error::Substitution)?,
                    ))
                }
            }
        }
    }

    /// Gets the reference nucleotide from the [`Relation`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::Relation;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let relation = Relation::<Nucleotide>::try_new(Some(Nucleotide::A), Some(Nucleotide::C))?;
    /// assert_eq!(relation.reference(), Some(&Nucleotide::A));
    ///
    /// let relation = Relation::<Nucleotide>::try_new(Some(Nucleotide::A), None)?;
    /// assert_eq!(relation.reference(), Some(&Nucleotide::A));
    ///
    /// let relation = Relation::<Nucleotide>::try_new(None, Some(Nucleotide::C))?;
    /// assert_eq!(relation.reference(), None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> Option<&N> {
        match self {
            Relation::Identical(reference) => Some(reference),
            Relation::Substitution(substitution) => Some(substitution.reference()),
            Relation::Insertion(_) => None,
            Relation::Deletion(reference) => Some(reference),
        }
    }

    /// Gets the alternate nucleotide from the [`Relation`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::Relation;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let relation = Relation::<Nucleotide>::try_new(Some(Nucleotide::A), Some(Nucleotide::C))?;
    /// assert_eq!(relation.alternate(), Some(&Nucleotide::C));
    ///
    /// let relation = Relation::<Nucleotide>::try_new(Some(Nucleotide::A), None)?;
    /// assert_eq!(relation.alternate(), None);
    ///
    /// let relation = Relation::<Nucleotide>::try_new(None, Some(Nucleotide::C))?;
    /// assert_eq!(relation.alternate(), Some(&Nucleotide::C));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> Option<&N> {
        match self {
            Relation::Identical(reference) => Some(reference),
            Relation::Substitution(substitution) => Some(substitution.alternate()),
            Relation::Insertion(alternate) => Some(alternate),
            Relation::Deletion(_) => None,
        }
    }

    /// Returns a reference to the [`Substitution`] wrapped in [`Some`] if the
    /// [`Relation`] is of kind [`Relation::Substitution`]. Else, [`None`] is
    /// returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::Relation;
    /// use omics_molecule::compound::nucleotide::relation::Substitution;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let relation = Relation::try_new(Some(Nucleotide::A), Some(Nucleotide::T))?;
    /// assert_eq!(
    ///     relation.as_substitution(),
    ///     Some(&Substitution::try_new(Nucleotide::A, Nucleotide::T)?)
    /// );
    ///
    /// let relation = Relation::try_new(None, Some(Nucleotide::T))?;
    /// assert!(relation.as_substitution().is_none());
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn as_substitution(&self) -> Option<&Substitution<N>> {
        match self {
            Relation::Substitution(substitution) => Some(substitution),
            _ => None,
        }
    }

    /// Consumes `self` and returns the [`Nucleotide(s)`] that comprise this
    /// [`Relation`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::Relation;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// // An insertion.
    /// let relation = Relation::try_new(None, Some(Nucleotide::T))?;
    /// assert_eq!(relation.into_nucleotides(), (None, Some(Nucleotide::T)));
    ///
    /// // A deletion.
    /// let relation = Relation::try_new(Some(Nucleotide::A), None)?;
    /// assert_eq!(relation.into_nucleotides(), (Some(Nucleotide::A), None));
    ///
    /// // An identical pair of nucleotides.
    /// let relation = Relation::try_new(Some(Nucleotide::A), Some(Nucleotide::A))?;
    /// assert_eq!(
    ///     relation.into_nucleotides(),
    ///     (Some(Nucleotide::A), Some(Nucleotide::A))
    /// );
    ///
    /// // A substitution.
    /// let relation = Relation::try_new(Some(Nucleotide::A), Some(Nucleotide::T))?;
    /// assert_eq!(
    ///     relation.into_nucleotides(),
    ///     (Some(Nucleotide::A), Some(Nucleotide::T))
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_nucleotides(self) -> (Option<N>, Option<N>) {
        match self {
            Relation::Identical(reference) => (Some(reference.clone()), Some(reference)),
            Relation::Substitution(substitution) => {
                let (reference, alternate) = substitution.into_parts();
                (Some(reference), Some(alternate))
            }
            Relation::Insertion(alternate) => (None, Some(alternate)),
            Relation::Deletion(reference) => (Some(reference), None),
        }
    }
}

impl<N: Nucleotide> From<Relation<N>> for (Option<N>, Option<N>) {
    fn from(value: Relation<N>) -> Self {
        let (reference, alternate) = value.into_nucleotides();
        (reference, alternate)
    }
}

impl<N: Nucleotide> std::str::FromStr for Relation<N> {
    type Err = Error<N>;

    fn from_str(s: &str) -> Result<Self, N> {
        let parts = s.split(VARIANT_SEPARATOR).collect::<Vec<_>>();

        if parts.len() != 2 {
            return Err(Error::ParseError(s.to_owned()));
        }

        let mut parts = parts.into_iter();

        // SAFETY: we just ensured above that the length will always be two.
        // Since we have not taken any items from the iterator, this item will
        // always unwrap.
        let reference_nucleotide = match parts.next().unwrap() {
            MISSING_NUCLEOTIDE => None,
            v => v
                .parse::<N>()
                .map(Some)
                .map_err(|_| Error::ParseError(s.to_owned()))?,
        };

        // SAFETY: we just ensured above that the length will always be two.
        // Since we have only taken one item from the iterator, this second item
        // will always unwrap.
        let alternate_nucleotide = match parts.next().unwrap() {
            MISSING_NUCLEOTIDE => None,
            v => v
                .parse::<N>()
                .map(Some)
                .map_err(|_| Error::ParseError(s.to_owned()))?,
        };

        Self::try_new(reference_nucleotide, alternate_nucleotide)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compound::nucleotide::relation::substitution::Kind;
    use crate::polymer::dna;
    use crate::polymer::rna;

    #[test]
    fn it_correctly_identifies_all_dna_transitions()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let relation = "A:G".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        let relation = "G:A".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        let relation = "C:T".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        let relation = "T:C".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        Ok(())
    }

    #[test]
    fn it_correctly_identifies_all_dna_transversions()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let relation = "A:C".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "C:A".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "A:T".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "T:A".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "G:C".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "C:G".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "G:T".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "T:G".parse::<Relation<dna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        Ok(())
    }

    #[test]
    fn it_correctly_identifies_all_rna_transitions()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let relation = "A:G".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        let relation = "G:A".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        let relation = "C:U".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        let relation = "U:C".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transition);

        Ok(())
    }

    #[test]
    fn it_correctly_identifies_all_rna_transversions()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let relation = "A:C".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "C:A".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "A:U".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "U:A".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "G:C".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "C:G".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "G:U".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        let relation = "U:G".parse::<Relation<rna::Nucleotide>>()?;
        let substitution = relation.as_substitution().unwrap();
        assert_eq!(substitution.kind(), Kind::Transversion);

        Ok(())
    }

    #[test]
    fn it_correctly_identifies_an_insertion() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        assert_eq!(
            ".:T".parse::<Relation<dna::Nucleotide>>()?,
            Relation::Insertion(dna::Nucleotide::T)
        );

        Ok(())
    }

    #[test]
    fn it_correctly_identifies_a_deletion() -> std::result::Result<(), Box<dyn std::error::Error>> {
        assert_eq!(
            "T:.".parse::<Relation<dna::Nucleotide>>()?,
            Relation::Deletion(dna::Nucleotide::T)
        );

        Ok(())
    }

    #[test]
    fn it_does_not_allow_an_empty_relation() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        assert_eq!(
            ".:."
                .parse::<Relation<dna::Nucleotide>>()
                .unwrap_err()
                .to_string(),
            "cannot create a relation with no nucleotides"
        );

        Ok(())
    }
}
