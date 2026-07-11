//! Combined deletion-insertions.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::coordinate;
use omics_coordinate::position::Number;
use omics_coordinate::system::Base;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence;
use omics_molecule::sequence::Sequence;

use crate::variant::Alteration;
use crate::variant::Kind;
use crate::variant::KindError;
use crate::variant::base_interval;

/// A combined deletion-insertion of differing lengths.
///
/// Serialized top-level variants use base coordinates with the `(b)`
/// qualifier. Interbase coordinates with `(i)` are rejected because a delins
/// replaces existing bases over a base interval.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The coordinate of the first replaced base.
    ///
    /// Construction guarantees that advancing this forward by
    /// `reference.len() - 1` stays within the position bounds, so
    /// [`interval`](Self::interval) cannot overflow.
    coordinate: Coordinate<Base>,

    /// Both alleles non-empty and of differing lengths.
    ///
    /// Construction guarantees `kind()` is [`Kind::Delins`].
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N> {
    /// Attempts to create a new delins from a coordinate and the reference and
    /// alternate bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::delins::Variant;
    ///
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl TryInto<Coordinate<Base>, Error = coordinate::Error>,
        reference: impl TryInto<Sequence<N>, Error = sequence::ParseError>,
        alternate: impl TryInto<Sequence<N>, Error = sequence::ParseError>,
    ) -> Result<Self, KindError> {
        let coordinate = coordinate.try_into()?;
        let reference = reference.try_into()?;
        let alternate = alternate.try_into()?;
        let alteration = Alteration::try_new(reference, alternate)?;
        Self::try_from((coordinate, alteration))
    }

    /// Gets the coordinate of the first replaced base.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::delins::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(variant.coordinate().position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn coordinate(&self) -> &Coordinate<Base> {
        &self.coordinate
    }

    /// Gets the reference allele.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::delins::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &Sequence<N> {
        self.alteration.reference()
    }

    /// Gets the alternate allele.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::delins::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(variant.alternate().to_string(), "G");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &Sequence<N> {
        self.alteration.alternate()
    }

    /// Gets the interval spanned by the reference allele.
    ///
    /// Use this for reference-facing overlap and annotation queries. This is
    /// the interval replaced on the reference.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::delins::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(variant.interval().end().position().get(), 101);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_interval(&self) -> Interval<Base> {
        // SAFETY: construction validated that this span is representable.
        base_interval(&self.coordinate, self.alteration.reference().len()).unwrap()
    }

    /// Gets the local interval spanned by the alternate allele.
    ///
    /// The alternate interval starts at the same coordinate as the reference
    /// interval and spans the alternate allele length. This describes the
    /// local replacement sequence, not a globally projected query coordinate.
    /// Use this for local alternate-allele span queries.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::delins::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(
    ///     variant
    ///         .alternate_interval()
    ///         .unwrap()
    ///         .start()
    ///         .position()
    ///         .get(),
    ///     100
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate_interval(&self) -> Option<Interval<Base>> {
        base_interval(&self.coordinate, self.alteration.alternate().len())
    }

    /// Gets the interval spanned by the reference allele.
    pub fn interval(&self) -> Interval<Base> {
        self.reference_interval()
    }

    /// Gets the underlying [`Alteration`] carrying both alleles.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::delins::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT", "G")?;
    /// assert_eq!(variant.alteration().alternate().to_string(), "G");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alteration(&self) -> &Alteration<N> {
        &self.alteration
    }
}

impl<N: Nucleotide> TryFrom<(Coordinate<Base>, Alteration<N>)> for Variant<N> {
    type Error = KindError;

    /// Builds a delins from a base [`Coordinate`] and a classified delins
    /// [`Alteration`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    /// use omics_variation::variant::delins::Variant;
    ///
    /// let coordinate = Coordinate::<Base>::try_from("seq0:+:100")?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "G".parse()?)?;
    /// let variant = Variant::try_from((coordinate, alteration))?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_from(
        (coordinate, alteration): (Coordinate<Base>, Alteration<N>),
    ) -> Result<Self, Self::Error> {
        let found = alteration.kind();
        if found != Kind::Delins {
            return Err(KindError::WrongKind {
                expected: Kind::Delins,
                found,
            });
        }

        // Validate the span is representable now so `interval` cannot panic.
        let span = Number::try_from(alteration.reference().len() - 1)
            .map_err(|_| KindError::SpanOverflow)?;
        coordinate
            .clone()
            .into_move_forward(span)
            .ok_or(KindError::SpanOverflow)?;

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
    fn it_builds_a_delins_from_raw_parts() {
        let variant = Variant::<dna::Nucleotide>::try_new("seq0:+:100", "AT", "G").unwrap();
        assert_eq!(variant.reference().to_string(), "AT");
        assert_eq!(variant.alternate().to_string(), "G");
        assert_eq!(variant.interval().end().position().get(), 101);
    }

    #[test]
    fn it_builds_a_delins_from_an_alteration() {
        let coordinate = Coordinate::<Base>::try_from("seq0:+:100").unwrap();
        let alteration = Alteration::<dna::Nucleotide>::try_new(
            Sequence::try_from("AT").unwrap(),
            Sequence::try_from("G").unwrap(),
        )
        .unwrap();
        let variant = Variant::try_from((coordinate, alteration)).unwrap();
        assert_eq!(variant.reference().to_string(), "AT");
    }

    #[test]
    fn it_rejects_non_delins() {
        let coordinate = Coordinate::<Base>::try_from("seq0:+:100").unwrap();
        // An MNV, not a delins.
        let alteration = Alteration::<dna::Nucleotide>::try_new(
            Sequence::try_from("AT").unwrap(),
            Sequence::try_from("GC").unwrap(),
        )
        .unwrap();
        let err = Variant::try_from((coordinate, alteration)).unwrap_err();
        assert!(matches!(
            err,
            KindError::WrongKind {
                expected: Kind::Delins,
                found: Kind::Mnv
            }
        ));
    }

    #[test]
    fn it_rejects_a_span_that_overflows() {
        // Two-base reference delins at MAX would need MAX+1.
        let coordinate = format!("seq0:+:{}", Number::MAX);
        let err = Variant::<dna::Nucleotide>::try_new(coordinate.as_str(), "AT", "G").unwrap_err();
        assert!(matches!(err, KindError::SpanOverflow));
    }
}
