//! Deletions.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::position::Number;
use omics_coordinate::system::Base;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;

use crate::variant::Alteration;
use crate::variant::Kind;
use crate::variant::KindError;

/// A deletion: one or more reference bases removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The coordinate of the first deleted base.
    ///
    /// [`try_new`](Self::try_new) guarantees that advancing this forward by
    /// `reference.len() - 1` stays within the position bounds, so
    /// [`interval`](Self::interval) cannot overflow.
    coordinate: Coordinate<Base>,

    /// A non-empty reference allele paired with an empty alternate allele.
    ///
    /// [`try_new`](Self::try_new) guarantees `kind()` is [`Kind::Deletion`].
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N> {
    /// Attempts to create a new deletion from an [`Alteration`].
    ///
    /// The span (`start + reference.len() - 1`) is validated here so that
    /// [`interval`](Self::interval) never panics on a constructed value.
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
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, ".".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl Into<Coordinate<Base>>,
        alteration: Alteration<N>,
    ) -> Result<Self, KindError> {
        let found = alteration.kind();
        if found != Kind::Deletion {
            return Err(KindError::WrongKind {
                expected: Kind::Deletion,
                found,
            });
        }

        let coordinate = coordinate.into();

        // Validate the span is representable now so `interval()` cannot panic.
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

    /// Gets the coordinate of the first deleted base.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::deletion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, ".".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
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
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::deletion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, ".".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &Sequence<N> {
        self.alteration.reference()
    }

    /// Gets the interval spanned by the deleted bases.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::deletion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, ".".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.interval().end().position().get(), 101);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn interval(&self) -> Interval<Base> {
        let span = Number::try_from(self.alteration.reference().len() - 1)
            .expect("validated at construction");
        let end = self
            .coordinate
            .clone()
            .into_move_forward(span)
            .expect("validated at construction");
        // SAFETY: start and end share a contig/strand and end >= start.
        Interval::try_new(self.coordinate.clone(), end).unwrap()
    }

    /// Gets the underlying [`Alteration`] carrying both alleles.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::deletion::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, ".".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.alteration().reference().to_string(), "AT");
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

    fn deletion(pos: &str, reference: &str) -> Variant<dna::Nucleotide> {
        let coordinate = format!("seq0:+:{pos}").parse::<Coordinate<Base>>().unwrap();
        let alteration =
            Alteration::try_new(reference.parse().unwrap(), ".".parse().unwrap()).unwrap();
        Variant::try_new(coordinate, alteration).unwrap()
    }

    #[test]
    fn it_builds_a_deletion_and_computes_its_interval() {
        let variant = deletion("100", "AT");
        assert_eq!(variant.reference().to_string(), "AT");
        let interval = variant.interval();
        assert_eq!(interval.start().position().get(), 100);
        assert_eq!(interval.end().position().get(), 101);
    }

    #[test]
    fn it_rejects_a_non_deletion_alteration() {
        let coordinate = "seq0:+:100".parse::<Coordinate<Base>>().unwrap();
        // An SNV, not a deletion.
        let alteration = Alteration::try_new("A".parse().unwrap(), "C".parse().unwrap()).unwrap();
        let err = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap_err();
        assert!(matches!(
            err,
            KindError::WrongKind {
                expected: Kind::Deletion,
                found: Kind::Snv
            }
        ));
    }

    #[test]
    fn it_rejects_a_span_that_overflows() {
        let coordinate = format!("seq0:+:{}", Number::MAX)
            .parse::<Coordinate<Base>>()
            .unwrap();
        // Two-base deletion at MAX would need MAX+1.
        let alteration = Alteration::try_new("AT".parse().unwrap(), ".".parse().unwrap()).unwrap();
        let err = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap_err();
        assert!(matches!(err, KindError::SpanOverflow));
    }
}
