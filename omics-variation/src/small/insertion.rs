//! Insertions.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::coordinate;
use omics_coordinate::system::Base;
use omics_coordinate::system::Interbase;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence;
use omics_molecule::sequence::Sequence;

use crate::small::Alteration;
use crate::small::Kind;
use crate::small::KindError;
use crate::small::base_interval;
use crate::small::interbase_interval;

/// An insertion of one or more bases at an interbase boundary.
///
/// Serialized top-level variants use interbase coordinates with the `(i)`
/// qualifier. Base coordinates with `(b)` are rejected because an insertion
/// occurs at a boundary between bases rather than on an existing base.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The interbase boundary at which the alternate allele is inserted.
    coordinate: Coordinate<Interbase>,

    /// An empty reference allele paired with a non-empty alternate allele.
    ///
    /// Construction guarantees `kind()` is [`Kind::Insertion`].
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N> {
    /// Attempts to create a new insertion from an interbase coordinate and the
    /// inserted bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::insertion::Variant;
    ///
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl TryInto<Coordinate<Interbase>, Error = coordinate::Error>,
        alternate: impl TryInto<Sequence<N>, Error = sequence::ParseError>,
    ) -> Result<Self, KindError> {
        let coordinate = coordinate.try_into()?;
        let alternate = alternate.try_into()?;
        let alteration = Alteration::try_new(Sequence::new(Vec::new()), alternate)?;
        Self::try_from((coordinate, alteration))
    }

    /// Gets the interbase boundary of the insertion.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::insertion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
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
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::insertion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &Sequence<N> {
        self.alteration.alternate()
    }

    /// Gets the zero-width interval spanned by the reference allele.
    ///
    /// Use this for reference-facing overlap and annotation queries. This is
    /// the interbase boundary where the insertion occurs.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::insertion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// let interval = variant.interbase_interval();
    /// assert_eq!(
    ///     interval.start().position().get(),
    ///     interval.end().position().get()
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_interval(&self) -> Interval<Interbase> {
        interbase_interval(&self.coordinate)
    }

    /// Gets the local interval spanned by the alternate allele.
    ///
    /// The alternate interval starts at the base immediately after the
    /// insertion boundary and spans the inserted allele length. This describes
    /// the local inserted sequence, not a globally projected query coordinate.
    /// Use this for local alternate-allele span queries.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::insertion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(
    ///     variant
    ///         .alternate_interval()
    ///         .unwrap()
    ///         .start()
    ///         .position()
    ///         .get(),
    ///     101
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate_interval(&self) -> Option<Interval<Base>> {
        let start = self.coordinate.clone().nudge_forward()?;
        base_interval(&start, self.alteration.alternate().len())
    }

    /// Gets the zero-width interbase interval at the insertion boundary.
    pub fn interbase_interval(&self) -> Interval<Interbase> {
        self.reference_interval()
    }

    /// Gets the underlying [`Alteration`] carrying both alleles.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::insertion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.alteration().alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alteration(&self) -> &Alteration<N> {
        &self.alteration
    }
}

impl<N: Nucleotide> TryFrom<(Coordinate<Interbase>, Alteration<N>)> for Variant<N> {
    type Error = KindError;

    /// Builds an insertion from an interbase [`Coordinate`] and a classified
    /// insertion [`Alteration`].
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
    /// let coordinate = Coordinate::<Interbase>::try_from("seq0:+:100")?;
    /// let alteration = Alteration::<Nucleotide>::try_new(".".parse()?, "AT".parse()?)?;
    /// let variant = Variant::try_from((coordinate, alteration))?;
    /// assert_eq!(variant.alternate().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_from(
        (coordinate, alteration): (Coordinate<Interbase>, Alteration<N>),
    ) -> Result<Self, Self::Error> {
        let found = alteration.kind();
        if found != Kind::Insertion {
            return Err(KindError::WrongKind {
                expected: Kind::Insertion,
                found,
            });
        }

        Ok(Self {
            coordinate,
            alteration,
        })
    }
}

#[cfg(test)]
mod tests {
    use omics_molecule::polymer::dna;

    use super::*;

    #[test]
    fn it_builds_an_insertion_from_raw_parts() {
        let variant = Variant::<dna::Nucleotide>::try_new("seq0:+:100", "AT").unwrap();
        assert_eq!(variant.alternate().to_string(), "AT");
    }

    #[test]
    fn it_builds_an_insertion_from_an_alteration() {
        let coordinate = Coordinate::<Interbase>::try_from("seq0:+:100").unwrap();
        let alteration = Alteration::<dna::Nucleotide>::try_new(
            Sequence::try_from(".").unwrap(),
            Sequence::try_from("AT").unwrap(),
        )
        .unwrap();
        let variant = Variant::try_from((coordinate, alteration)).unwrap();
        assert_eq!(variant.alternate().to_string(), "AT");
    }

    #[test]
    fn it_reports_a_zero_width_interbase_interval() {
        let variant = Variant::<dna::Nucleotide>::try_new("seq0:+:100", "AT").unwrap();
        let interval = variant.interbase_interval();
        assert_eq!(interval.start().position().get(), 100);
        assert_eq!(interval.end().position().get(), 100);
    }

    #[test]
    fn it_rejects_a_non_insertion() {
        let coordinate = Coordinate::<Interbase>::try_from("seq0:+:100").unwrap();
        // A deletion, not an insertion.
        let alteration = Alteration::<dna::Nucleotide>::try_new(
            Sequence::try_from("A").unwrap(),
            Sequence::try_from(".").unwrap(),
        )
        .unwrap();
        let err = Variant::try_from((coordinate, alteration)).unwrap_err();
        assert!(matches!(
            err,
            KindError::WrongKind {
                expected: Kind::Insertion,
                found: Kind::Deletion
            }
        ));
    }
}
