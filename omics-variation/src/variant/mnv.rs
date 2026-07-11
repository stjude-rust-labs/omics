//! Multi-nucleotide variants.

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::position::Number;
use omics_coordinate::system::Base;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;

use crate::variant::Alteration;
use crate::variant::Kind;
use crate::variant::KindError;

/// A multi-nucleotide variant: an equal-length substitution of two or more
/// bases.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The coordinate of the first substituted base.
    ///
    /// [`try_new`](Self::try_new) guarantees that advancing this forward by
    /// `reference.len() - 1` stays within the position bounds, so
    /// [`interval`](Self::interval) cannot overflow.
    coordinate: Coordinate<Base>,

    /// Equal-length, non-empty, differing reference and alternate alleles.
    ///
    /// [`try_new`](Self::try_new) guarantees `kind()` is [`Kind::Mnv`].
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N> {
    /// Attempts to create a new MNV from an [`Alteration`].
    ///
    /// The span is validated here so that [`interval`](Self::interval) never
    /// panics on a constructed value.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::variant::Alteration;
    /// use omics_variation::variant::mnv::Variant;
    ///
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "GC".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.reference().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl Into<Coordinate<Base>>,
        alteration: Alteration<N>,
    ) -> Result<Self, KindError> {
        let found = alteration.kind();
        if found != Kind::Mnv {
            return Err(KindError::WrongKind {
                expected: Kind::Mnv,
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

    /// Gets the coordinate of the first substituted base.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::mnv::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "GC".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
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
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::mnv::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "GC".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
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
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::mnv::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "GC".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.alternate().to_string(), "GC");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &Sequence<N> {
        self.alteration.alternate()
    }

    /// Gets the interval spanned by the reference bases.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_coordinate::Coordinate;
    /// # use omics_coordinate::system::Base;
    /// # use omics_molecule::polymer::dna::Nucleotide;
    /// # use omics_variation::variant::Alteration;
    /// # use omics_variation::variant::mnv::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "GC".parse()?)?;
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
    /// # use omics_variation::variant::mnv::Variant;
    /// let coordinate = "seq0:+:100".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("AT".parse()?, "GC".parse()?)?;
    /// let variant = Variant::try_new(coordinate, alteration)?;
    /// assert_eq!(variant.alteration().alternate().to_string(), "GC");
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
    fn it_builds_an_mnv() {
        let coordinate = "seq0:+:100".parse::<Coordinate<Base>>().unwrap();
        let alteration = Alteration::try_new("AT".parse().unwrap(), "GC".parse().unwrap()).unwrap();
        let variant = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap();
        assert_eq!(variant.reference().to_string(), "AT");
        assert_eq!(variant.alternate().to_string(), "GC");
        assert_eq!(variant.interval().end().position().get(), 101);
    }

    #[test]
    fn it_rejects_non_mnv() {
        let coordinate = "seq0:+:100".parse::<Coordinate<Base>>().unwrap();
        // A delins, not an MNV.
        let alteration = Alteration::try_new("AT".parse().unwrap(), "G".parse().unwrap()).unwrap();
        let err = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap_err();
        assert!(matches!(
            err,
            KindError::WrongKind {
                expected: Kind::Mnv,
                found: Kind::Delins
            }
        ));
    }

    #[test]
    fn it_rejects_a_span_that_overflows() {
        let coordinate = format!("seq0:+:{}", Number::MAX)
            .parse::<Coordinate<Base>>()
            .unwrap();
        // Two-base MNV at MAX would need MAX+1.
        let alteration = Alteration::try_new("AT".parse().unwrap(), "GC".parse().unwrap()).unwrap();
        let err = Variant::<dna::Nucleotide>::try_new(coordinate, alteration).unwrap_err();
        assert!(matches!(err, KindError::SpanOverflow));
    }
}
