//! Generic nucleotide sequences.
//!
//! A [`Sequence`] is the generic sibling of the concrete
//! [`dna::Molecule`](crate::polymer::dna::Molecule) — an ordered run of
//! [`Nucleotide`]s used to represent a variant allele. An empty sequence is a
//! valid allele (it denotes a missing side of an insertion or deletion) and is
//! rendered/parsed using [`omics_core::MISSING_NUCLEOTIDE`].

use std::str::FromStr;

use omics_core::MISSING_NUCLEOTIDE;
use thiserror::Error;

use crate::compound::Nucleotide;

/// An error encountered while parsing a [`Sequence`].
#[derive(Error, Debug)]
pub enum ParseError {
    /// The token was empty (use `.` to denote a missing allele).
    #[error("empty sequence token (use `.` for a missing allele)")]
    EmptyToken,

    /// A nucleotide within the token could not be parsed.
    #[error("invalid nucleotide in sequence: `{0}`")]
    Nucleotide(String),
}

/// An ordered run of [`Nucleotide`]s representing a variant allele.
///
/// An empty [`Sequence`] is valid and denotes a missing allele.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence<N: Nucleotide>(Vec<N>);

impl<N: Nucleotide> Sequence<N> {
    /// Creates a new [`Sequence`] from a [`Vec`] of nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let seq = Sequence::new(vec![Nucleotide::A, Nucleotide::C]);
    /// assert_eq!(seq.len(), 2);
    /// ```
    pub fn new(inner: Vec<N>) -> Self {
        Self(inner)
    }

    /// Gets the nucleotides of this [`Sequence`] as a slice.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let seq = Sequence::new(vec![Nucleotide::A]);
    /// assert_eq!(seq.inner(), &[Nucleotide::A]);
    /// ```
    pub fn inner(&self) -> &[N] {
        &self.0
    }

    /// Consumes `self` and returns the inner [`Vec`] of nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let seq = Sequence::new(vec![Nucleotide::A]);
    /// assert_eq!(seq.into_inner(), vec![Nucleotide::A]);
    /// ```
    pub fn into_inner(self) -> Vec<N> {
        self.0
    }

    /// Gets the number of nucleotides in this [`Sequence`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let seq = Sequence::new(vec![Nucleotide::A, Nucleotide::C]);
    /// assert_eq!(seq.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether this [`Sequence`] is empty (a missing allele).
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let seq = Sequence::<Nucleotide>::new(vec![]);
    /// assert!(seq.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Counts the leading nucleotides shared with `other`.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let a = "AAT".parse::<Sequence<Nucleotide>>()?;
    /// let b = "AAC".parse::<Sequence<Nucleotide>>()?;
    /// assert_eq!(a.shared_prefix_len(&b), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn shared_prefix_len(&self, other: &Sequence<N>) -> usize {
        self.0
            .iter()
            .zip(other.0.iter())
            .take_while(|(a, b)| a == b)
            .count()
    }

    /// Counts the trailing nucleotides shared with `other`, excluding any
    /// already counted as a shared prefix (so prefix + suffix never exceeds the
    /// shorter length).
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_molecule::sequence::Sequence;
    ///
    /// let a = "ATG".parse::<Sequence<Nucleotide>>()?;
    /// let b = "ACG".parse::<Sequence<Nucleotide>>()?;
    /// assert_eq!(a.shared_suffix_len(&b), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn shared_suffix_len(&self, other: &Sequence<N>) -> usize {
        let prefix = self.shared_prefix_len(other);
        let max = std::cmp::min(self.0.len(), other.0.len()) - prefix;
        self.0
            .iter()
            .rev()
            .zip(other.0.iter().rev())
            .take(max)
            .take_while(|(a, b)| a == b)
            .count()
    }
}

impl<N: Nucleotide> FromStr for Sequence<N> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::EmptyToken);
        }

        if s == MISSING_NUCLEOTIDE {
            return Ok(Self(Vec::new()));
        }

        s.chars()
            .map(|c| {
                c.to_string()
                    .parse::<N>()
                    .map_err(|_| ParseError::Nucleotide(c.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Self)
    }
}

impl<N: Nucleotide> std::fmt::Display for Sequence<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.0.is_empty() {
            return write!(f, "{MISSING_NUCLEOTIDE}");
        }

        for nucleotide in &self.0 {
            write!(f, "{nucleotide}")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polymer::dna;

    #[test]
    fn it_parses_and_displays_a_multibase_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let seq = "ACGT".parse::<Sequence<dna::Nucleotide>>()?;
        assert_eq!(seq.len(), 4);
        assert!(!seq.is_empty());
        assert_eq!(seq.to_string(), "ACGT");
        Ok(())
    }

    #[test]
    fn it_parses_and_displays_the_missing_allele() -> Result<(), Box<dyn std::error::Error>> {
        let seq = ".".parse::<Sequence<dna::Nucleotide>>()?;
        assert!(seq.is_empty());
        assert_eq!(seq.to_string(), ".");
        Ok(())
    }

    #[test]
    fn it_rejects_an_empty_token() {
        let err = "".parse::<Sequence<dna::Nucleotide>>().unwrap_err();
        assert!(matches!(err, ParseError::EmptyToken));
    }

    #[test]
    fn it_rejects_an_invalid_nucleotide() {
        let err = "AQ".parse::<Sequence<dna::Nucleotide>>().unwrap_err();
        assert!(matches!(err, ParseError::Nucleotide(_)));
    }

    #[test]
    fn it_computes_shared_prefix_and_suffix() -> Result<(), Box<dyn std::error::Error>> {
        let a = "AATG".parse::<Sequence<dna::Nucleotide>>()?;
        let b = "AACG".parse::<Sequence<dna::Nucleotide>>()?;
        assert_eq!(a.shared_prefix_len(&b), 2); // "AA"
        assert_eq!(a.shared_suffix_len(&b), 1); // "G"
        Ok(())
    }

    #[test]
    fn shared_prefix_and_suffix_do_not_overlap_count() -> Result<(), Box<dyn std::error::Error>> {
        // Identical single base: prefix counts it, suffix must not double-count.
        let a = "A".parse::<Sequence<dna::Nucleotide>>()?;
        let b = "A".parse::<Sequence<dna::Nucleotide>>()?;
        assert_eq!(a.shared_prefix_len(&b), 1);
        assert_eq!(a.shared_suffix_len(&b), 0);
        Ok(())
    }
}
