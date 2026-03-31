//! Protein.

mod amino_acid;

pub use amino_acid::AminoAcid;
use thiserror::Error;

/// An error related to a [`Molecule`].
#[derive(Error, Debug)]
pub enum Error {
    /// An error when processing an [`AminoAcid`].
    #[error(transparent)]
    AminoAcidError(#[from] amino_acid::Error),
}

/// A molecule representing a protein, which is a polymer of amino acids.
///
/// # Examples
///
/// ```
/// use omics_molecule::polymer::protein::Molecule;
///
/// let m = "MKWVTFISLLFLFSSAYS".parse::<Molecule>()?;
/// assert_eq!(m.inner().len(), 18);
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct Molecule(Vec<AminoAcid>);

impl Molecule {
    /// Gets the inner [`Vec<AminoAcid>`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::protein::Molecule;
    ///
    /// let m = "ARNDCEQGHILKMFPSTWYV".parse::<Molecule>()?;
    /// assert_eq!(m.inner().len(), 20);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn inner(&self) -> &Vec<AminoAcid> {
        self.0.as_ref()
    }

    /// Consumes the [`Molecule`] and returns the inner [`Vec<AminoAcid>`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::protein::AminoAcid;
    /// use omics_molecule::polymer::protein::Molecule;
    ///
    /// let m = "ACDE".parse::<Molecule>()?;
    /// let amino_acids = m.into_inner();
    ///
    /// assert_eq!(
    ///     amino_acids,
    ///     vec![AminoAcid::A, AminoAcid::C, AminoAcid::D, AminoAcid::E,]
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_inner(self) -> Vec<AminoAcid> {
        self.0
    }
}

impl From<Vec<AminoAcid>> for Molecule {
    fn from(v: Vec<AminoAcid>) -> Self {
        Self(v)
    }
}

impl std::str::FromStr for Molecule {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.chars()
            .map(|c| AminoAcid::try_from(c).map_err(Error::AminoAcidError))
            .collect::<Result<Vec<_>, Error>>()
            .map(Self::from)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_creates_a_molecule_from_a_vec_of_amino_acids() {
        let amino_acids = vec![AminoAcid::M, AminoAcid::K, AminoAcid::W, AminoAcid::V];

        let protein = Molecule::from(amino_acids);
        assert_eq!(protein.inner().len(), 4);
    }

    #[test]
    fn it_parses_a_molecule_from_a_valid_string() -> Result<(), Box<dyn std::error::Error>> {
        Ok("ARNDCEQGHILKMFPSTWYV".parse::<Molecule>().map(|_| ())?)
    }

    #[test]
    fn it_fails_to_parse_a_molecule_from_an_invalid_string() {
        let err = "ARJX".parse::<Molecule>().unwrap_err();
        assert_eq!(err.to_string(), "invalid amino acid `J`");
    }
}
