//! Single nucleotide variations.

use std::str::FromStr;

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::Strand;
use omics_coordinate::coordinate;
use omics_coordinate::system::Base;
use omics_core::VARIANT_SEPARATOR;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;
use thiserror::Error;

use crate::PositionQualifier;
use crate::qualified_position;
use crate::split_qualified_position;
use crate::variant::Alteration;
use crate::variant::Kind;
use crate::variant::KindError;
use crate::variant::base_interval;

/// A parse error related to a [`Variant`].
#[derive(Error, Debug)]
pub enum ParseError<N: Nucleotide>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    /// An invalid format was encountered when parsing a [`Variant`].
    #[error("invalid format: {0}")]
    InvalidFormat(String),

    /// An issue occurred when parsing the coordinate of the [`Variant`].
    #[error(transparent)]
    CoordinateError(#[from] coordinate::Error),

    /// The coordinate system qualifier was missing or invalid.
    #[error("position `{position}` must end with `(b)`")]
    CoordinateSystemQualifier {
        /// The unqualified or incorrectly qualified position token.
        position: String,
    },

    /// An issue occurred when parsing the reference nucleotide of the
    /// [`Variant`].
    #[error("reference nucleotide error: {0}")]
    ReferenceNucleotide(<N as FromStr>::Err),

    /// An issue occurred when parsing the alternate nucleotide of the
    /// [`Variant`].
    #[error("alternate nucleotide error: {0}")]
    AlternateNucleotide(<N as FromStr>::Err),
}

/// An error related to a [`Variant`].
#[derive(Error, Debug)]
pub enum Error<N: Nucleotide>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    /// Attempted to create a [`Variant`] with identical reference and
    /// alternate nucleotides.
    #[error("identical nucleotides for snv: {0}")]
    Identical(N),

    /// Unsuccessfully attempted to parse a [`Variant`] from a string.
    #[error(transparent)]
    Parse(#[from] ParseError<N>),

    /// The alleles did not form a valid alteration.
    #[error(transparent)]
    Alteration(crate::variant::Error),
}

/// A single nucleotide variant.
///
/// Serialized `SNV` strings use base coordinates with the `(b)` qualifier.
/// Interbase coordinates with `(i)` are rejected because an `SNV` substitutes
/// one existing base.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant<N: Nucleotide> {
    /// The coordinate of the substituted base.
    coordinate: Coordinate<Base>,

    /// Reference and alternate alleles, each exactly one base.
    ///
    /// Construction guarantees `kind()` is [`Kind::Snv`].
    alteration: Alteration<N>,
}

impl<N: Nucleotide> Variant<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    /// Attempts to create a new [`Variant`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::base::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = Variant::<dna::Nucleotide>::try_new(
    ///     "seq0:+:1".parse::<Coordinate>()?,
    ///     dna::Nucleotide::A,
    ///     dna::Nucleotide::T,
    /// )?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        coordinate: impl Into<Coordinate<Base>>,
        reference_nucleotide: impl Into<N>,
        alternate_nucleotide: impl Into<N>,
    ) -> Result<Self, Error<N>> {
        let coordinate = coordinate.into();
        let reference_nucleotide = reference_nucleotide.into();
        let alternate_nucleotide = alternate_nucleotide.into();

        let alteration = Alteration::try_new(
            Sequence::new(vec![reference_nucleotide]),
            Sequence::new(vec![alternate_nucleotide]),
        )
        .map_err(|err| match err {
            crate::variant::Error::Identical(_) => Error::Identical(reference_nucleotide),
            other => Error::Alteration(other),
        })?;

        Ok(Self {
            coordinate,
            alteration,
        })
    }
}

impl<N: Nucleotide> Variant<N> {
    /// Gets the [`Coordinate`] for this [`Variant`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::base::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1(b):A:T".parse::<Variant<dna::Nucleotide>>()?;
    ///
    /// assert_eq!(variant.coordinate().contig().as_str(), "seq0");
    /// assert_eq!(variant.coordinate().strand(), Strand::Positive);
    /// assert_eq!(variant.coordinate().position().get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn coordinate(&self) -> &Coordinate<Base> {
        &self.coordinate
    }

    /// Gets the reference nucleotide as a [`Nucleotide`] from the [`Variant`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1(b):A:T".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.reference(), dna::Nucleotide::A);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> N {
        // SAFETY: a single nucleotide variant always has exactly one reference
        // base.
        self.alteration.reference().inner()[0]
    }

    /// Gets the alternate nucleotide as a [`Nucleotide`] from the [`Variant`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1(b):A:T".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.alternate(), dna::Nucleotide::T);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> N {
        // SAFETY: a single nucleotide variant always has exactly one alternate
        // base.
        self.alteration.alternate().inner()[0]
    }

    /// Gets the interval spanned by the reference allele.
    ///
    /// Use this for reference-facing overlap and annotation queries.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
    /// let interval = variant.interval();
    /// assert_eq!(interval.start().position().get(), 1);
    /// assert_eq!(interval.end().position().get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_interval(&self) -> Interval<Base> {
        // SAFETY: an SNV reference allele always has exactly one base.
        base_interval(&self.coordinate, self.alteration.reference().len()).unwrap()
    }

    /// Gets the interval spanned by the alternate allele.
    ///
    /// The alternate interval is the same base as the reference interval
    /// because an `SNV` substitutes one existing base with one alternate base.
    /// Use this for local alternate-allele span queries.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
    /// let interval = variant.alternate_interval();
    /// assert_eq!(interval.start().position().get(), 1);
    /// assert_eq!(interval.end().position().get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate_interval(&self) -> Interval<Base> {
        self.reference_interval()
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
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.alteration().reference().to_string(), "A");
    /// assert_eq!(variant.alteration().alternate().to_string(), "C");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alteration(&self) -> &Alteration<N> {
        &self.alteration
    }
}

impl<N: Nucleotide> TryFrom<(Coordinate<Base>, Alteration<N>)> for Variant<N> {
    type Error = KindError;

    /// Builds a single nucleotide [`Variant`] from a base [`Coordinate`] and a
    /// classified single-base [`Alteration`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna::Nucleotide;
    /// use omics_variation::snv::Variant;
    /// use omics_variation::variant::Alteration;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// let alteration = Alteration::<Nucleotide>::try_new("A".parse()?, "C".parse()?)?;
    /// let variant = Variant::try_from((coordinate, alteration))?;
    /// assert_eq!(variant.reference(), Nucleotide::A);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_from(
        (coordinate, alteration): (Coordinate<Base>, Alteration<N>),
    ) -> Result<Self, Self::Error> {
        let found = alteration.kind();
        if found != Kind::Snv {
            return Err(KindError::WrongKind {
                expected: Kind::Snv,
                found,
            });
        }

        Ok(Self {
            coordinate,
            alteration,
        })
    }
}

impl<N: Nucleotide> std::str::FromStr for Variant<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    type Err = Error<N>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(VARIANT_SEPARATOR).collect::<Vec<_>>();
        let num_parts = parts.len();

        if num_parts != 4 && num_parts != 5 {
            return Err(Error::Parse(ParseError::InvalidFormat(s.to_owned())));
        }

        let mut parts = parts.into_iter();

        let (contig, strand, position) = match num_parts {
            4 => {
                // SAFETY: we just ensured that the number of parts is four.
                // Since we have not taken any items from the iterator, these
                // two items will always unwrap.
                let contig = parts.next().unwrap();
                let position = parts.next().unwrap();
                (contig, Strand::Positive.to_string(), position)
            }
            5 => {
                // SAFETY: we just ensured that the number of parts is five.
                // Since we have not taken any items from the iterator, these
                // three items will always unwrap.
                let contig = parts.next().unwrap();
                let strand = parts.next().unwrap();
                let position = parts.next().unwrap();
                (contig, strand.to_string(), position)
            }
            // SAFETY: we ensured above that the number of parts must be either four or five.
            _ => unreachable!(),
        };
        let (position, qualifier) = split_qualified_position(position).ok_or_else(|| {
            Error::Parse(ParseError::CoordinateSystemQualifier {
                position: position.to_owned(),
            })
        })?;

        if qualifier != PositionQualifier::Base {
            return Err(Error::Parse(ParseError::CoordinateSystemQualifier {
                position: format!("{position}{}", qualifier.suffix()),
            }));
        }
        let coordinate = [contig, strand.as_str(), position].join(VARIANT_SEPARATOR);

        let coordinate = match coordinate.parse::<Coordinate<Base>>() {
            Ok(coordinate) => coordinate,
            Err(err) => return Err(Error::Parse(ParseError::CoordinateError(err))),
        };

        // SAFETY: in all cases above, we leave two items in the iterator. Since we have
        // not taken any items yet, this will always unwrap.
        let reference_nucleotide = parts
            .next()
            .unwrap()
            .parse::<N>()
            .map_err(|err| Error::Parse(ParseError::ReferenceNucleotide(err)))?;

        // SAFETY: in all cases above, we leave two items in the iterator. Since we have
        // only taken one item so far, this will always unwrap.
        let alternate_nucleotide = parts
            .next()
            .unwrap()
            .parse::<N>()
            .map_err(|err| Error::Parse(ParseError::AlternateNucleotide(err)))?;

        Self::try_new(coordinate, reference_nucleotide, alternate_nucleotide)
    }
}

impl<N: Nucleotide> std::fmt::Display for Variant<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let coordinate = [
            self.coordinate().contig().to_string(),
            self.coordinate().strand().to_string(),
            qualified_position(self.coordinate().position().get(), PositionQualifier::Base),
        ]
        .join(VARIANT_SEPARATOR);

        let parts = [
            coordinate.as_str(),
            &self.reference().to_string(),
            &self.alternate().to_string(),
        ];

        write!(f, "{}", parts.join(self::VARIANT_SEPARATOR))
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::Strand;
    use omics_molecule::polymer::dna;
    use omics_molecule::polymer::rna;

    use super::*;

    #[test]
    fn it_creates_a_variant_in_a_dna_context() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:+:1(b):A:C".parse::<Variant<dna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), dna::Nucleotide::A);
        assert_eq!(variant.alternate(), dna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_in_a_rna_context() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:+:1(b):U:C".parse::<Variant<rna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), rna::Nucleotide::U);
        assert_eq!(variant.alternate(), rna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_on_the_negative_strand_in_a_dna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:-:1(b):A:C".parse::<Variant<dna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Negative);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), dna::Nucleotide::A);
        assert_eq!(variant.alternate(), dna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_on_the_negative_strand_in_a_rna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:-:1(b):U:C".parse::<Variant<rna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Negative);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), rna::Nucleotide::U);
        assert_eq!(variant.alternate(), rna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_with_no_specified_strand_in_a_dna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:1(b):A:C".parse::<Variant<dna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), dna::Nucleotide::A);
        assert_eq!(variant.alternate(), dna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_with_no_specified_strand_in_a_rna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:1(b):U:C".parse::<Variant<rna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), rna::Nucleotide::U);
        assert_eq!(variant.alternate(), rna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_fails_when_creating_a_variant_with_identical_nucleotides() {
        let err = "seq0:+:1(b):A:A"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(err.to_string(), "identical nucleotides for snv: A");
    }

    #[test]
    fn it_fails_when_attempting_to_represent_an_insertion() {
        let err = "seq0:+:1(b):.:A"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "reference nucleotide error: invalid nucleotide `.`"
        );
    }

    #[test]
    fn it_fails_when_attempting_to_represent_a_deletion() {
        let err = "seq0:+:1(b):A:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "alternate nucleotide error: invalid nucleotide `.`"
        );
    }

    #[test]
    fn it_fails_when_attempting_to_represent_an_empty_pair() {
        let err = "seq0:+:1(b):.:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "reference nucleotide error: invalid nucleotide `.`"
        );
    }

    #[test]
    fn it_fails_without_a_base_coordinate_system_qualifier() {
        let err = "seq0:+:1:A:C"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(err.to_string(), "position `1` must end with `(b)`");
    }

    #[test]
    fn it_fails_with_an_interbase_coordinate_system_qualifier() {
        let err = "seq0:+:1(i):A:C"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(err.to_string(), "position `1(i)` must end with `(b)`");
    }
}
