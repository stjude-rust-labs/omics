//! Contiguous molecules.

use thiserror::Error;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// An error related to a contig.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// An empty contig name was provided.
    #[error("contig name cannot be empty")]
    Empty,
}

/// A [`Result`](std::result::Result) with an [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// Contig
////////////////////////////////////////////////////////////////////////////////////////

/// A named, contiguous molecule within a genome.
///
/// At present, a contig is simply a wrapper around a non-empty string.
///
/// Notably, the internal representation of [`Contig`] may change in the future
/// (though the interface to this type will remain stable with respect to
/// [semantic versioning](https://semver.org/)).
///
/// For a more in-depth discussion on this, please see [this section of the
/// docs](crate#contigs).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Contig(String);

impl Contig {
    /// Attempts to create a new contig.
    ///
    /// Returns an error if the contig name is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = Contig::try_new("chr1")?;
    /// assert_eq!(contig.as_str(), "chr1");
    ///
    /// // Empty contig names are rejected
    /// assert!(Contig::try_new("").is_err());
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(value: impl Into<String>) -> Result<Self> {
        let s = value.into();
        if s.is_empty() {
            return Err(Error::Empty);
        }
        Ok(Self(s))
    }

    /// Creates a new contig without validating that the name is non-empty.
    ///
    /// # Safety
    ///
    /// This function does not validate that the contig name is non-empty.
    /// Creating a contig with an empty name may lead to invalid serialization
    /// and parsing errors.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = Contig::new_unchecked("chr1");
    /// assert_eq!(contig.as_str(), "chr1");
    /// ```
    pub fn new_unchecked(value: impl Into<String>) -> Self {
        Self(value.into())
    }

    // NOTE: an `inner()` method is explicitly not included as the type
    // dereferences `String`. This means that the `as_str()` method is usable
    // for this purpose.

    /// Consumes `self` and returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = Contig::new_unchecked("chr1");
    /// assert_eq!(contig.into_inner(), String::from("chr1"));
    /// ```
    pub fn into_inner(self) -> String {
        self.0
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl std::fmt::Display for Contig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::str::FromStr for Contig {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_new(s)
    }
}

impl TryFrom<&str> for Contig {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self> {
        Self::try_new(value)
    }
}

impl TryFrom<String> for Contig {
    type Error = Error;

    fn try_from(value: String) -> Result<Self> {
        Self::try_new(value)
    }
}

impl std::ops::Deref for Contig {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_new_valid() {
        let contig = Contig::try_new("chr1").expect("valid contig name");
        assert_eq!(contig.as_str(), "chr1");

        let contig = Contig::try_new("seq0").expect("valid contig name");
        assert_eq!(contig.as_str(), "seq0");

        let contig = Contig::try_new("X").expect("valid contig name");
        assert_eq!(contig.as_str(), "X");
    }

    #[test]
    fn try_new_empty() {
        let err = Contig::try_new("").expect_err("empty contig name should fail");
        assert_eq!(err, Error::Empty);
        assert_eq!(err.to_string(), "contig name cannot be empty");
    }

    #[test]
    fn new_unchecked() {
        let contig = Contig::new_unchecked("chr1");
        assert_eq!(contig.as_str(), "chr1");

        // new_unchecked allows empty strings (though not recommended)
        let contig = Contig::new_unchecked("");
        assert_eq!(contig.as_str(), "");
    }

    #[test]
    fn parse() {
        let contig = "chr1".parse::<Contig>().expect("contig to parse");
        assert_eq!(contig.as_str(), "chr1");
    }

    #[test]
    fn parse_empty() {
        let err = "".parse::<Contig>().expect_err("empty string should fail");
        assert_eq!(err, Error::Empty);
    }

    #[test]
    fn try_from_str() {
        let contig = Contig::try_from("chr1").expect("valid contig");
        assert_eq!(contig.as_str(), "chr1");

        let err = Contig::try_from("").expect_err("empty should fail");
        assert_eq!(err, Error::Empty);
    }

    #[test]
    fn try_from_string() {
        let contig = Contig::try_from(String::from("chr1")).expect("valid contig");
        assert_eq!(contig.as_str(), "chr1");

        let err = Contig::try_from(String::from("")).expect_err("empty should fail");
        assert_eq!(err, Error::Empty);
    }
}
