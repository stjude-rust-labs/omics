//! Deletions.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::coordinate;
use omics_coordinate::position::Number;
use omics_coordinate::system::Base;
use omics_coordinate::system::Interbase;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence;
use omics_molecule::sequence::Sequence;

use crate::variant::Alteration;
use crate::variant::Kind;
use crate::variant::KindError;
use crate::variant::base_interval;
use crate::variant::interbase_interval_before_base;

/// A deletion of one or more reference bases.
///
/// Serialized top-level variants use base coordinates with the `(b)`
/// qualifier. Interbase coordinates with `(i)` are rejected because a deletion
/// removes existing bases over a base interval.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The coordinate of the first deleted base.
    ///
    /// Construction guarantees that advancing this forward by
    /// `reference.len() - 1` stays within the position bounds, so
    /// [`interval`](Self::interval) cannot overflow.
    coordinate: Coordinate<Base>,

    /// A non-empty reference allele paired with an empty alternate allele.
    ///
    /// Construction guarantees `kind()` is [`Kind::Deletion`].
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N> {
    /// Attempts to create a new deletion from a coordinate and the deleted
    /// reference bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::deletion::Variant;
    ///
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl TryInto<Coordinate<Base>, Error = coordinate::Error>,
        reference: impl TryInto<Sequence<N>, Error = sequence::ParseError>,
    ) -> Result<Self, KindError> {
        let coordinate = coordinate.try_into()?;
        let reference = reference.try_into()?;
        let alteration = Alteration::try_new(reference, Sequence::new(Vec::new()))?;
        Self::try_from((coordinate, alteration))
    }

    /// Gets the coordinate of the first deleted base.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::deletion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.coordinate().position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn coordinate(&self) -> &Coordinate<Base> {
        &self.coordinate
    }

    /// Gets the deleted reference allele.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::deletion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &Sequence<N> {
        self.alteration.reference()
    }

    /// Gets the interval spanned by the reference allele.
    ///
    /// Use this for reference-facing overlap and annotation queries. This is
    /// the interval deleted from the reference.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::deletion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.interval().end().position().get(), 101);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_interval(&self) -> Interval<Base> {
        // SAFETY: construction validated that this span is representable.
        base_interval(&self.coordinate, self.alteration.reference().len()).unwrap()
    }

    /// Gets the zero-width interval spanned by the alternate allele.
    ///
    /// A deletion has an empty alternate allele. Its alternate interval is the
    /// interbase boundary immediately before the deleted reference interval.
    /// Use this for local alternate-allele span queries.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::deletion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(
    ///     variant.alternate_interval().start().position().get(),
    ///     variant.alternate_interval().end().position().get()
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate_interval(&self) -> Interval<Interbase> {
        interbase_interval_before_base(&self.coordinate)
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
    /// # use omics_variation::variant::deletion::Variant;
    /// let variant = Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
    /// assert_eq!(variant.alteration().reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alteration(&self) -> &Alteration<N> {
        &self.alteration
    }
}

impl<N: Nucleotide> TryFrom<(Coordinate<Base>, Alteration<N>)> for Variant<N> {
    type Error = KindError;

    /// Builds a deletion from a base [`Coordinate`] and a classified deletion
    /// [`Alteration`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    /// use omics_variation::variant::deletion::Variant;
    ///
    /// let coordinate = Coordinate::<Base>::try_from("seq0:+:100")?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, ".".parse()?)?;
    /// let variant = Variant::try_from((coordinate, alteration))?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_from(
        (coordinate, alteration): (Coordinate<Base>, Alteration<N>),
    ) -> Result<Self, Self::Error> {
        let found = alteration.kind();
        if found != Kind::Deletion {
            return Err(KindError::WrongKind {
                expected: Kind::Deletion,
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
    fn it_builds_a_deletion_from_raw_parts() {
        let variant = Variant::<dna::Nucleotide>::try_new("seq0:+:100", "AT").unwrap();
        assert_eq!(variant.reference().to_string(), "AT");
        assert_eq!(variant.interval().start().position().get(), 100);
        assert_eq!(variant.interval().end().position().get(), 101);
    }

    #[test]
    fn it_builds_a_deletion_from_an_alteration() {
        let coordinate = Coordinate::<Base>::try_from("seq0:+:100").unwrap();
        let alteration = Alteration::<dna::Nucleotide>::try_new(
            Sequence::try_from("AT").unwrap(),
            Sequence::try_from(".").unwrap(),
        )
        .unwrap();
        let variant = Variant::try_from((coordinate, alteration)).unwrap();
        assert_eq!(variant.reference().to_string(), "AT");
    }

    #[test]
    fn it_rejects_a_non_deletion_alteration() {
        let coordinate = Coordinate::<Base>::try_from("seq0:+:100").unwrap();
        // An SNV, not a deletion.
        let alteration = Alteration::<dna::Nucleotide>::try_new(
            Sequence::try_from("A").unwrap(),
            Sequence::try_from("C").unwrap(),
        )
        .unwrap();
        let err = Variant::try_from((coordinate, alteration)).unwrap_err();
        assert!(matches!(
            err,
            KindError::WrongKind {
                expected: Kind::Deletion,
                found: Kind::Snv
            }
        ));
    }

    #[test]
    fn it_rejects_an_empty_reference() {
        let err = Variant::<dna::Nucleotide>::try_new("seq0:+:100", ".").unwrap_err();
        assert!(matches!(
            err,
            KindError::Alteration(crate::variant::Error::BothEmpty)
        ));
    }

    #[test]
    fn it_rejects_a_span_that_overflows() {
        // Two-base deletion at MAX would need MAX+1.
        let coordinate = format!("seq0:+:{}", Number::MAX);
        let err = Variant::<dna::Nucleotide>::try_new(coordinate.as_str(), "AT").unwrap_err();
        assert!(matches!(err, KindError::SpanOverflow));
    }
}
