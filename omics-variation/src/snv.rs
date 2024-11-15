//! Single nucleotide variations.

use std::str::FromStr;

use omics_coordinate::Coordinate;
use omics_coordinate::Strand;
use omics_coordinate::coordinate;
use omics_coordinate::system::Base;
use omics_core::VARIANT_SEPARATOR;
use omics_molecule::compound::Nucleotide;
use omics_molecule::compound::nucleotide::relation;
use omics_molecule::compound::nucleotide::relation::Relation;

/// A parse error related to a [`Variant`].
#[derive(Debug)]
pub enum ParseError<N: Nucleotide>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    /// An invalid format was encountered when parsing a [`Variant`].
    InvalidFormat(String),

    /// An issue occurred when parsing the coordinate of the [`Variant`].
    CoordinateError(coordinate::Error),

    /// An issue occurred when parsing the reference nucleotide of the
    /// [`Variant`].
    ReferenceNucleotide(<N as FromStr>::Err),

    /// An issue occurred when parsing the alternate nucleotide of the
    /// [`Variant`].
    AlternateNucleotide(<N as FromStr>::Err),
}

impl<N: Nucleotide> std::fmt::Display for ParseError<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InvalidFormat(value) => write!(f, "invalid format: {value}"),
            ParseError::CoordinateError(err) => write!(f, "coordinate error: {err}"),
            ParseError::ReferenceNucleotide(err) => write!(f, "reference nucleotide error: {err}"),
            ParseError::AlternateNucleotide(err) => write!(f, "alternate nucleotide error: {err}"),
        }
    }
}

impl<N: Nucleotide> std::error::Error for ParseError<N> where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display
{
}

/// An error related to a [`Variant`].
#[derive(Debug)]
pub enum Error<N: Nucleotide>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    /// Attempted to create a [`Variant`] with identical reference and
    /// alternate nucleotides.
    Identical(N),

    /// Unsuccessfully attempted to parse a [`Variant`] from a string.
    Parse(ParseError<N>),

    /// An error constructing a relation.
    Relation(relation::Error<N>),
}

impl<N: Nucleotide> std::fmt::Display for Error<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Identical(nucleotide) => {
                write!(f, "identical nucleotides for snv: {nucleotide}")
            }
            Error::Parse(err) => write!(f, "parse error: {err}"),
            Error::Relation(err) => write!(f, "relation error: {err}"),
        }
    }
}

impl<N: Nucleotide> std::error::Error for Error<N> where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display
{
}

/// A single nucleotide variant.
#[derive(Debug)]
pub struct Variant<N: Nucleotide> {
    /// The coordinate.
    coordinate: Coordinate<Base>,

    /// The relation.
    relation: Relation<N>,
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

        let relation = Relation::try_new(Some(reference_nucleotide), Some(alternate_nucleotide))
            .map_err(Error::Relation)?;

        if let Relation::Identical(nucleotide) = relation {
            return Err(Error::Identical(nucleotide));
        }

        Ok(Self {
            coordinate,
            relation,
        })
    }

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
    /// let variant = Variant::<dna::Nucleotide>::try_new(
    ///     "seq0:+:1".parse::<Coordinate>()?,
    ///     dna::Nucleotide::A,
    ///     dna::Nucleotide::T,
    /// )?;
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
    /// use omics_coordinate::base::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1:A:T".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.reference(), &dna::Nucleotide::A);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &N {
        // SAFETY: because a single nucleotide variant is guaranteed to have a
        // reference nucleotide within the inner [`Relation`], this will
        // always unwrap successfully.
        self.relation.reference().unwrap()
    }

    /// Gets the alternate nucleotide as a [`Nucleotide`] from the [`Variant`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::base::Coordinate;
    /// use omics_coordinate::system::Base;
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::snv::Variant;
    ///
    /// let variant = "seq0:+:1:A:T".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.alternate(), &dna::Nucleotide::T);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate(&self) -> &N {
        // SAFETY: because a single nucleotide variant is guaranteed to have a
        // alternate nucleotide within the inner [`Relation`], this will
        // always unwrap successfully.
        self.relation.alternate().unwrap()
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

        let coordinate = match num_parts {
            4 => {
                let positive = Strand::Positive.to_string();

                // SAFETY: we just ensured that the number of parts is four.
                // Since we have not taken any items from the iterator, these
                // two items will always unwrap.
                [
                    parts.next().unwrap(),
                    positive.as_str(),
                    parts.next().unwrap(),
                ]
                .join(VARIANT_SEPARATOR)
            }
            5 => {
                // SAFETY: we just ensured that the number of parts is five.
                // Since we have not taken any items from the iterator, these
                // three items will always unwrap.
                [
                    parts.next().unwrap(),
                    parts.next().unwrap(),
                    parts.next().unwrap(),
                ]
                .join(VARIANT_SEPARATOR)
            }
            // SAFETY: we ensured above that the number of parts must be either four or five.
            _ => unreachable!(),
        };

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

impl<N: Nucleotide> std::fmt::Display for Variant<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let coordinate = self.coordinate().to_string();

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
    use omics_molecule::polymer::dna;
    use omics_molecule::polymer::rna;

    use super::*;

    #[test]
    fn it_creates_a_variant_in_a_dna_context() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:+:1:A:C".parse::<Variant<dna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), &dna::Nucleotide::A);
        assert_eq!(variant.alternate(), &dna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_in_a_rna_context() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:+:1:U:C".parse::<Variant<rna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), &rna::Nucleotide::U);
        assert_eq!(variant.alternate(), &rna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_on_the_negative_strand_in_a_dna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:-:1:A:C".parse::<Variant<dna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Negative);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), &dna::Nucleotide::A);
        assert_eq!(variant.alternate(), &dna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_on_the_negative_strand_in_a_rna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:-:1:U:C".parse::<Variant<rna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Negative);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), &rna::Nucleotide::U);
        assert_eq!(variant.alternate(), &rna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_with_no_specified_strand_in_a_dna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:1:A:C".parse::<Variant<dna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), &dna::Nucleotide::A);
        assert_eq!(variant.alternate(), &dna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_creates_a_variant_with_no_specified_strand_in_a_rna_context()
    -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:1:U:C".parse::<Variant<rna::Nucleotide>>()?;

        assert_eq!(variant.coordinate().contig().as_str(), "seq0");
        assert_eq!(variant.coordinate().strand(), Strand::Positive);
        assert_eq!(variant.coordinate().position().get(), 1);
        assert_eq!(variant.reference(), &rna::Nucleotide::U);
        assert_eq!(variant.alternate(), &rna::Nucleotide::C);

        Ok(())
    }

    #[test]
    fn it_fails_when_creating_a_variant_with_identical_nucleotides() {
        let err = "seq0:+:1:A:A"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(err.to_string(), "identical nucleotides for snv: A");
    }

    #[test]
    fn it_fails_when_attempting_to_represent_an_insertion() {
        let err = "seq0:+:1:.:A"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: reference nucleotide error: parse error: invalid nucleotide: ."
        );
    }

    #[test]
    fn it_fails_when_attempting_to_represent_a_deletion() {
        let err = "seq0:+:1:A:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: alternate nucleotide error: parse error: invalid nucleotide: ."
        );
    }

    #[test]
    fn it_fails_when_attempting_to_represent_an_empty_pair() {
        let err = "seq0:+:1:.:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: reference nucleotide error: parse error: invalid nucleotide: ."
        );
    }
}
