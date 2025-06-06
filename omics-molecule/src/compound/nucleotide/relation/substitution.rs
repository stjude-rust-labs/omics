//! Substitutions.

use crate::compound::Nucleotide;

mod kind;

pub use kind::Kind;
use thiserror::Error;

/// An error related to a [`Substitution`].
#[derive(Error, Debug)]
pub enum Error<N: Nucleotide> {
    /// Attempted to create a [`Substitution`] with identical reference and
    /// alternate nucleotides.
    #[error("identical nucleotides in substitution: {0}")]
    Identical(N),
}

/// A [`Result`](std::result::Result) with an [`Error<N>`].
type Result<T, N> = std::result::Result<T, Error<N>>;

/// The substitution of a reference nucleotide wth an alternate nucleotide.
#[derive(Debug, Eq, PartialEq)]
pub struct Substitution<N: Nucleotide> {
    /// The reference allele.
    reference: N,

    /// The alternative allele.
    alternate: N,
}

impl<N: Nucleotide> Substitution<N> {
    /// Creates a new [`Substitution`].
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::relation::Substitution;
    /// use omics_molecule::compound::nucleotide::relation::substitution::Kind;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;
    ///
    /// assert_eq!(substitution.kind(), Kind::Transversion);
    /// assert_eq!(substitution.reference(), &Nucleotide::A);
    /// assert_eq!(substitution.alternate(), &Nucleotide::T);
    /// assert_eq!(substitution.into_parts(), (Nucleotide::A, Nucleotide::T));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(reference: N, alternate: N) -> Result<Self, N> {
        if reference == alternate {
            return Err(Error::Identical(reference));
        }

        Ok(Self {
            reference,
            alternate,
        })
    }

    /// Gets the [`Kind`] for this [`Substitution`].
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::relation::Substitution;
    /// use omics_molecule::compound::nucleotide::relation::substitution::Kind;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;
    ///
    /// assert_eq!(substitution.kind(), Kind::Transversion);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn kind(&self) -> Kind {
        if self.reference.kind() == self.alternate.kind() {
            Kind::Transition
        } else {
            Kind::Transversion
        }
    }

    /// Gets the reference nucleotide from this [`Substitution`].
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::relation::Substitution;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;
    ///
    /// assert_eq!(substitution.reference(), &Nucleotide::A);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &N {
        &self.reference
    }

    /// Gets the alternate nucleotide from this [`Substitution`].
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::relation::Substitution;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;
    ///
    /// assert_eq!(substitution.alternate(), &Nucleotide::T);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &N {
        &self.alternate
    }

    /// Breaks down a [`Substitution`] into its consitutient parts.
    ///
    /// ```
    /// use omics_molecule::compound::nucleotide::relation::Substitution;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;
    ///
    /// assert_eq!(substitution.into_parts(), (Nucleotide::A, Nucleotide::T));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_parts(self) -> (N, N) {
        (self.reference, self.alternate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polymer::dna::Nucleotide;

    #[test]
    fn it_correctly_creates_a_substitution() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;

        assert_eq!(substitution.kind(), Kind::Transversion);
        assert_eq!(substitution.reference(), &Nucleotide::A);
        assert_eq!(substitution.alternate(), &Nucleotide::T);
        assert_eq!(substitution.into_parts(), (Nucleotide::A, Nucleotide::T));

        Ok(())
    }

    #[test]
    fn it_correctly_refuses_to_create_a_substitution_with_identical_nucleotides()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let substitution = Substitution::try_new(Nucleotide::A, Nucleotide::T)?;

        assert_eq!(substitution.kind(), Kind::Transversion);
        assert_eq!(substitution.reference(), &Nucleotide::A);
        assert_eq!(substitution.alternate(), &Nucleotide::T);
        assert_eq!(substitution.into_parts(), (Nucleotide::A, Nucleotide::T));

        Ok(())
    }
}
