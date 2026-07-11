//! Insertions.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::system::Interbase;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;

use crate::variant::Alteration;
use crate::variant::Kind;
use crate::variant::KindError;

/// An insertion: one or more bases inserted at an interbase boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The interbase boundary at which the alternate allele is inserted.
    coordinate: Coordinate<Interbase>,

    /// The alleles (empty reference, non-empty alternate).
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N> {
    /// Attempts to create a new insertion from an [`Alteration`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Interbase;
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    /// use omics_variation::variant::insertion::Variant;
    ///
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new(".".parse()?, "AT".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl Into<Coordinate<Interbase>>,
        alteration: Alteration<N>,
    ) -> Result<Self, KindError> {
        let found = alteration.kind();
        if found != Kind::Insertion {
            return Err(KindError::WrongKind {
                expected: Kind::Insertion,
                found,
            });
        }

        Ok(Self {
            coordinate: coordinate.into(),
            alteration,
        })
    }

    /// Gets the interbase boundary of the insertion.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Interbase;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::insertion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new(".".parse()?, "AT".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.coordinate().position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn coordinate(&self) -> &Coordinate<Interbase> {
        &self.coordinate
    }

    /// Gets the inserted alternate allele.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Interbase;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::insertion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new(".".parse()?, "AT".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &Sequence<N> {
        self.alteration.alternate()
    }

    /// Gets the zero-width interbase interval at the insertion boundary.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Interbase;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::insertion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new(".".parse()?, "AT".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// let interval = variant.interbase_interval();
    /// assert_eq!(
    ///     interval.start().position().get(),
    ///     interval.end().position().get()
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn interbase_interval(&self) -> Interval<Interbase> {
        // SAFETY: a zero-width interbase interval (start == end) is always valid.
        Interval::try_new(self.coordinate.clone(), self.coordinate.clone()).unwrap()
    }

    /// Gets the underlying [`Alteration`] carrying both alleles.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Interbase;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::insertion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new(".".parse()?, "AT".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.alteration().alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alteration(&self) -> &Alteration<N> {
        &self.alteration
    }
}

#[cfg(test)]
mod tests {
    use omics_molecule::polymer::dna;

    use super::*;

    #[test]
    fn it_builds_an_insertion() {
        let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>().unwrap();
        let alteration = Alteration::try_new(".".parse().unwrap(), "AT".parse().unwrap()).unwrap();
        let variant = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap();
        assert_eq!(variant.alternate().to_string(), "AT");
    }

    #[test]
    fn it_reports_a_zero_width_interbase_interval() {
        let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>().unwrap();
        let alteration = Alteration::try_new(".".parse().unwrap(), "AT".parse().unwrap()).unwrap();
        let variant = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap();
        let interval = variant.interbase_interval();
        assert_eq!(interval.start().position().get(), 100);
        assert_eq!(interval.end().position().get(), 100);
    }

    #[test]
    fn it_rejects_a_non_insertion() {
        let coordinate = "seq0:+:100".parse::<Coordinate<Interbase>>().unwrap();
        // A deletion, not an insertion.
        let alteration = Alteration::try_new("A".parse().unwrap(), ".".parse().unwrap()).unwrap();
        let err = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap_err();
        assert!(matches!(
            err,
            KindError::WrongKind {
                expected: Kind::Insertion,
                found: Kind::Deletion
            }
        ));
    }
}
