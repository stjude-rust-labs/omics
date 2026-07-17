//! Breakend orientation.

use std::str::FromStr;

use thiserror::Error;

/// A parse error related to an [`Orientation`].
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// An invalid orientation glyph was encountered.
    #[error("invalid orientation: `{0}`")]
    Invalid(String),
}

/// Which flank of a breakend's boundary is retained.
///
/// The flank is named by reference coordinate, not by strand. A breakend has no
/// strand, so [`LowerFlank`](Orientation::LowerFlank) and
/// [`HigherFlank`](Orientation::HigherFlank) always refer to the lower and
/// higher reference-coordinate side of the boundary. Whether the retained piece
/// is read forward or reverse-complemented in the derivative molecule is not
/// recorded here; that emerges from how two breakends pair, as described in the
/// [module documentation](crate::structural).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Orientation {
    /// The retained reference is the lower-coordinate flank.
    ///
    /// The junction extends onward from the boundary toward higher coordinates.
    LowerFlank,

    /// The retained reference is the higher-coordinate flank.
    HigherFlank,
}

impl Orientation {
    /// Gets the canonical ordering rank of this [`Orientation`].
    ///
    /// [`Orientation::LowerFlank`] ranks before [`Orientation::HigherFlank`].
    /// The rank is defined explicitly so canonical ordering does not depend on
    /// the incidental declaration order of the variants.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// assert!(Orientation::LowerFlank.rank() < Orientation::HigherFlank.rank());
    /// ```
    pub fn rank(&self) -> u8 {
        match self {
            Orientation::LowerFlank => 0,
            Orientation::HigherFlank => 1,
        }
    }

    /// Reports whether this [`Orientation`] is opposite to another.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// assert!(Orientation::LowerFlank.is_opposite(Orientation::HigherFlank));
    /// ```
    pub fn is_opposite(&self, other: Orientation) -> bool {
        *self != other
    }
}

impl FromStr for Orientation {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ">" => Ok(Orientation::LowerFlank),
            "<" => Ok(Orientation::HigherFlank),
            _ => Err(ParseError::Invalid(s.to_string())),
        }
    }
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Orientation::LowerFlank => write!(f, ">"),
            Orientation::HigherFlank => write!(f, "<"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_parses_and_serializes() {
        assert_eq!(">".parse::<Orientation>().unwrap(), Orientation::LowerFlank);
        assert_eq!(
            "<".parse::<Orientation>().unwrap(),
            Orientation::HigherFlank
        );
        assert_eq!(Orientation::LowerFlank.to_string(), ">");
        assert_eq!(Orientation::HigherFlank.to_string(), "<");
    }

    #[test]
    fn it_reports_opposite_orientations() {
        assert!(Orientation::LowerFlank.is_opposite(Orientation::HigherFlank));
        assert!(!Orientation::LowerFlank.is_opposite(Orientation::LowerFlank));
    }

    #[test]
    fn lower_flank_ranks_before_higher_flank() {
        assert!(Orientation::LowerFlank.rank() < Orientation::HigherFlank.rank());
    }
}
