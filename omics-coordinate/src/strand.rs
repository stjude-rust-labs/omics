//! Strands of a coordinate within a genome.

/// An error related to the parsing of a [`Strand`].
#[derive(Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Attempted to create an empty [`Strand`].
    Empty,

    /// A [`Strand`] was attempted to be parsed from the provided, invalid
    /// value.
    InvalidValue(String),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::Empty => write!(f, "empty strand"),
            ParseError::InvalidValue(value) => {
                write!(f, "invalid value for strand: {value}")
            }
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a [`Strand`].
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// A parse error.
    ParseError(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::ParseError(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A strand within a genome.
///
/// For a more in-depth discussion on this, please see [this section of the
/// docs](crate#strand).
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Strand {
    /// The positive strand (`+`).
    Positive,

    /// The negative strand (`-`).
    Negative,
}

impl Strand {
    /// Attempts to create a new [`Strand`].
    ///
    /// # Notes
    ///
    /// * This will fail with a [`ParseError::Empty`] if the provided [`String`]
    ///   is empty (a strand with no name is considered non-sensical).
    /// * This will fail with a [`ParseError::InvalidValue`] if the provided
    ///   [`String`] does not match a known [`Strand`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Strand;
    ///
    /// let contig = Strand::try_new("+")?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(value: &str) -> Result<Self, Error> {
        if value.is_empty() {
            return Err(Error::ParseError(ParseError::Empty));
        }

        match value {
            "+" => Ok(Self::Positive),
            "-" => Ok(Self::Negative),
            _ => Err(Error::ParseError(ParseError::InvalidValue(
                value.to_string(),
            ))),
        }
    }

    /// Complements a [`Strand`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Strand;
    ///
    /// // Positive complement
    ///
    /// let strand = Strand::Positive;
    /// assert_eq!(strand.complement(), Strand::Negative);
    ///
    /// // Negative complement
    ///
    /// let strand = Strand::Negative;
    /// assert_eq!(strand.complement(), Strand::Positive);
    /// ```
    pub fn complement(self) -> Strand {
        match self {
            Strand::Positive => Strand::Negative,
            Strand::Negative => Strand::Positive,
        }
    }
}

impl std::str::FromStr for Strand {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Strand::try_new(s)
    }
}

impl TryFrom<&str> for Strand {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Strand::try_new(value)
    }
}

impl TryFrom<String> for Strand {
    type Error = Error;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        Strand::try_new(&value)
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Positive => write!(f, "+"),
            Strand::Negative => write!(f, "-"),
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn it_deserializes_correctly() -> Result<(), Box<dyn std::error::Error>> {
        let strand: Strand = "+".parse()?;
        assert_eq!(strand, Strand::Positive);

        let strand: Strand = "-".parse()?;
        assert_eq!(strand, Strand::Negative);

        let err = "?".parse::<Strand>().unwrap_err();
        assert_eq!(err.to_string(), "parse error: invalid value for strand: ?");

        Ok(())
    }

    #[test]
    fn it_serializes_correctly() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(Strand::Positive.to_string(), "+");
        assert_eq!(Strand::Negative.to_string(), "-");

        Ok(())
    }
}
