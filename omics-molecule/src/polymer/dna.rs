//! Deoxyribonucleic Acid.

mod nucleotide;

pub use nucleotide::Nucleotide;

/// An error related to a [`Molecule`].
#[derive(Debug)]
pub enum Error {
    /// An error when processing a [`Nucleotide`].
    NucleotideError(nucleotide::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::NucleotideError(err) => write!(f, "nucleotide error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A molecule representing Deoxyribonucleic Acid, otherwise known as DNA.
#[derive(Debug)]
pub struct Molecule(Vec<Nucleotide>);

impl Molecule {
    /// Gets the inner [`Vec<Nucleotide>`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Molecule;
    ///
    /// let m = "ACGT".parse::<Molecule>()?;
    /// assert_eq!(m.inner().len(), 4);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn inner(&self) -> &Vec<Nucleotide> {
        self.0.as_ref()
    }

    /// Consumes the [`Molecule`] and returns the inner [`Vec<Nucleotide>`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Molecule;
    /// use omics_molecule::polymer::dna::Nucleotide;
    ///
    /// let m = "ACGT".parse::<Molecule>()?;
    /// let nucleotides = m.into_inner();
    ///
    /// assert_eq!(nucleotides, vec![
    ///     Nucleotide::A,
    ///     Nucleotide::C,
    ///     Nucleotide::G,
    ///     Nucleotide::T,
    /// ]);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_inner(self) -> Vec<Nucleotide> {
        self.0
    }

    /// Gets the GC content of this [`Molecule`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Molecule;
    ///
    /// let m = "ACGT".parse::<Molecule>()?;
    /// assert_eq!(m.gc_content(), 0.5);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn gc_content(&self) -> f32 {
        let numerator = self
            .0
            .iter()
            .filter(|n| *n == &Nucleotide::C || *n == &Nucleotide::G)
            .count();

        numerator as f32 / self.0.len() as f32
    }
}

impl From<Vec<Nucleotide>> for Molecule {
    fn from(v: Vec<Nucleotide>) -> Self {
        Self(v)
    }
}

impl std::str::FromStr for Molecule {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.chars()
            .map(|c| Nucleotide::try_from(c).map_err(Error::NucleotideError))
            .collect::<Result<Vec<_>, Error>>()
            .map(Self::from)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_creates_a_molecule_from_a_vec_of_nucleotides() {
        let nucleotides = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];

        let dna = Molecule::from(nucleotides);
        assert_eq!(dna.inner().len(), 4);
    }

    #[test]
    fn it_parses_a_molecule_from_a_valid_string() -> Result<(), Box<dyn std::error::Error>> {
        Ok("ACGT".parse::<Molecule>().map(|_| ())?)
    }

    #[test]
    fn it_fails_to_parse_a_molecule_from_an_invalid_string() {
        let err = "QQQQ".parse::<Molecule>().unwrap_err();
        assert_eq!(err.to_string(), "nucleotide error: invalid nucleotide: Q");
    }
}
