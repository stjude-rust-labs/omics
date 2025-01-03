//! Contiguous molecules.

use std::convert::Infallible;

////////////////////////////////////////////////////////////////////////////////////////
// Contig
////////////////////////////////////////////////////////////////////////////////////////

/// A named, contiguous molecule within a genome.
///
/// At present, a contig is simply a wrapper around a string. Empty contig names
/// are allowed though not recommended.
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
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = Contig::new("chr1");
    /// ```
    pub fn new(value: impl Into<String>) -> Self {
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
    /// let contig = Contig::new("chr1");
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
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Infallible> {
        Ok(Self::new(s))
    }
}

impl From<&str> for Contig {
    fn from(value: &str) -> Self {
        Self::new(value)
    }
}

impl From<String> for Contig {
    fn from(value: String) -> Self {
        Self::new(value)
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
    fn parse() {
        let contig = "chr1".parse::<Contig>().expect("contig to parse");
        assert_eq!(contig.as_str(), "chr1");
    }
}
