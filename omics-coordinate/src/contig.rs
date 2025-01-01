//! Contiguous molecules.

use std::convert::Infallible;

use string_interner::symbol::SymbolU32;

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
pub struct Contig(SymbolU32);

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
    pub fn new(value: &str) -> Self {
        let id = CORPUS.write().unwrap().get_or_intern(value);
        Self(id)
    }

    // NOTE: an `inner()` method is explicitly not included as the type
    // dereferences `String`. This means that the `as_str()` method is usable
    // for this purpose.

    /// Returns the value of the contig.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    ///
    /// let contig = Contig::new("chr1");
    /// assert_eq!(contig.value(), String::from("chr1"));
    /// ```
    pub fn value(&self) -> String {
        CORPUS.read().unwrap().resolve(self.0).unwrap().to_owned()
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl std::fmt::Display for Contig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse() {
        let contig = "chr1".parse::<Contig>().expect("contig to parse");
        assert_eq!(contig.value().as_str(), "chr1");
    }
}
