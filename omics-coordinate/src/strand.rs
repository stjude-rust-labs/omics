//! Strands of a molecule.

use thiserror::Error;

///////////////////////////////////////////////////////////////////////////////////////
// Assertions
////////////////////////////////////////////////////////////////////////////////////////

const _: () = {
    // A strand should always fit into a single byte.
    assert!(size_of::<Strand>() == 1);

    /// A function to ensure that types are `Copy`.
    const fn is_copy<T: Copy>() {}
    is_copy::<Strand>();
};

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// A error related to parsing a strand.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// An invalid strand.
    ///
    /// Occurs when a strand cannot be parsed.
    #[error("invalid strand: {value}")]
    Invalid {
        /// The value that was attempted to be parsed.
        value: String,
    },
}

/// A [`Result`](std::result::Result) with an [`ParseError`].
pub type ParseResult<T> = std::result::Result<T, ParseError>;

/// A strand-related error.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// A parse error.
    #[error("parse error: {0}")]
    Parse(#[from] ParseError),
}

/// A [`Result`](std::result::Result) with an [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// Strand
////////////////////////////////////////////////////////////////////////////////////////

/// The strand of a double-stranded molecule.
///
/// For a more in-depth discussion on this, please see [this section of the
/// docs](crate#strand).
#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strand {
    /// The positive strand (`+`).
    ///
    /// This is also known as the _sense_ strand.
    Positive,

    /// The negative strand (`-`).
    ///
    /// This is also known as the _antisense_ strand.
    Negative,
}

impl Strand {
    /// Complements a strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Strand;
    ///
    /// assert_eq!(Strand::Positive.complement(), Strand::Negative);
    /// assert_eq!(Strand::Negative.complement(), Strand::Positive);
    /// ```
    pub fn complement(&self) -> Strand {
        match self {
            Strand::Positive => Strand::Negative,
            Strand::Negative => Strand::Positive,
        }
    }
}

impl std::str::FromStr for Strand {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => Err(Error::Parse(ParseError::Invalid {
                value: s.to_string(),
            })),
        }
    }
}

// NOTE: technically, this is just a duplication of of [`FromStr`] above. That
// being said, this is required for the [`Coordinate::try_new()`] call to work
// correctly with string slices.

impl TryFrom<&str> for Strand {
    type Error = Error;

    fn try_from(value: &str) -> std::result::Result<Self, Self::Error> {
        value.parse()
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
mod tests {
    use super::*;

    #[test]
    fn parse() {
        let s = "+".parse::<Strand>().unwrap();
        assert_eq!(s, Strand::Positive);

        let s = "-".parse::<Strand>().unwrap();
        assert_eq!(s, Strand::Negative);

        let err = "a".parse::<Strand>().unwrap_err();
        assert_eq!(err.to_string(), "parse error: invalid strand: a");
    }

    #[test]
    fn serialize() {
        assert_eq!(Strand::Positive.to_string(), "+");
        assert_eq!(Strand::Negative.to_string(), "-");
    }

    #[test]
    fn complement() {
        assert_eq!(Strand::Positive.complement(), Strand::Negative);
        assert_eq!(Strand::Negative.complement(), Strand::Positive);
    }
}
