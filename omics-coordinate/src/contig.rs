//! Named, contiguous molecules within a genome.

use std::ops::Deref;
use std::str::FromStr;

/// An error related to the parsing of a [`Contig`].
#[derive(Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A [`Contig`] was attempted to be parsed from an invalid value.
    InvalidValue(String),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InvalidValue(value) => {
                write!(f, "invalid value for contig: {value}")
            }
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a [`Contig`].
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// Attempted to create an empty [`Contig`].
    Empty,

    /// A parse error.
    ParseError(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Empty => write!(f, "empty"),
            Error::ParseError(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A named, contiguous molecule within a genome.
///
/// At present, a [`Contig`] is simply a wrapper around [`String`]. Notably, the
/// internal representation of [`Contig`] may change in the future (though the
/// interface will remain stable with respect to [semantic versioning](https://semver.org/)).
///
/// For a more in-depth discussion on this, please see [this section of the
/// docs](crate#contigs).
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Contig(String);

impl Contig {
    /// Attempts to create a new [`Contig`].
    ///
    /// # Notes
    ///
    /// * This will fail with a [`Error::Empty`] if the provided [`String`] is
    ///   empty (empty contig names are considered non-sensical by this crate).
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = Contig::try_new("seq0")?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new<S: Into<String>>(value: S) -> Result<Self, Error> {
        let value = value.into();

        if value.is_empty() {
            return Err(Error::Empty);
        }

        Ok(Self(value))
    }

    /// Gets the inner [`String`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = "seq0".parse::<Contig>()?;
    /// assert_eq!(contig.inner(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn inner(&self) -> &str {
        self.0.as_str()
    }

    /// Consumes the [`Contig`] and returns the inner [`String`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = "seq0".parse::<Contig>()?;
    /// assert_eq!(contig.into_inner(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_inner(self) -> String {
        self.0
    }
}

impl std::fmt::Display for Contig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for Contig {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_new(s)
    }
}

impl TryFrom<&str> for Contig {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Self::try_new(value)
    }
}

impl TryFrom<String> for Contig {
    type Error = Error;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        Self::try_new(value)
    }
}

impl Deref for Contig {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_fails_to_create_a_contig_from_an_empty_string() {
        let err = Contig::try_from(String::new()).unwrap_err();
        assert!(matches!(err, Error::Empty));
    }
}
