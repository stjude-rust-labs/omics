//! Breakends.

use std::str::FromStr;

use omics_coordinate::Contig;
use omics_coordinate::Position;
use omics_coordinate::contig;
use omics_coordinate::position::Number;
use omics_coordinate::position::interbase::Position as InterbasePosition;
use omics_coordinate::system::Interbase;
use omics_core::VARIANT_SEPARATOR;
use thiserror::Error;

use crate::structural::orientation;
use crate::structural::orientation::Orientation;

/// A parse error related to a [`Breakend`].
#[derive(Error, Debug)]
pub enum ParseError {
    /// The breakend did not have the three colon-separated fields.
    #[error("invalid breakend format: `{0}`")]
    Format(String),

    /// The position did not end with the interbase `(i)` qualifier.
    #[error("position `{position}` must end with `(i)` for an interbase coordinate")]
    Qualifier {
        /// The offending position token.
        position: String,
    },

    /// The contig failed to parse.
    #[error(transparent)]
    Contig(#[from] contig::Error),

    /// The orientation failed to parse.
    #[error(transparent)]
    Orientation(#[from] orientation::ParseError),

    /// The position failed to parse.
    #[error(transparent)]
    Position(#[from] omics_coordinate::position::Error),
}

/// One oriented endpoint of a novel adjacency.
///
/// A breakend mirrors a coordinate's fields, except that the strand slot is
/// replaced by an [`Orientation`]. A breakend has no meaningful biological
/// strand, so it does not carry one.
///
/// `Hash` is intentionally not derived because `Position<Interbase>` does not
/// implement it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Breakend {
    /// The contig the breakend sits on.
    contig: Contig,

    /// Which flank of the boundary is retained.
    orientation: Orientation,

    /// The interbase boundary the breakend sits at.
    position: Position<Interbase>,
}

impl Breakend {
    /// Attempts to create a new [`Breakend`].
    ///
    /// The `position` is a raw [`Number`] rather than a generic `impl Into`,
    /// following the crate convention that non-trait constructors take
    /// positions concretely so call sites do not need type disambiguation.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// assert_eq!(breakend.position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        contig: impl TryInto<Contig, Error = contig::Error>,
        orientation: Orientation,
        position: Number,
    ) -> Result<Self, contig::Error> {
        Ok(Self {
            contig: contig.try_into()?,
            orientation,
            position: InterbasePosition::new(position),
        })
    }

    /// Gets the contig this breakend sits on.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// assert_eq!(breakend.contig().as_str(), "seq0");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        &self.contig
    }

    /// Gets the orientation of this breakend.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// assert_eq!(breakend.orientation(), Orientation::LowerFlank);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// Gets the interbase position of this breakend.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// assert_eq!(breakend.position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn position(&self) -> &Position<Interbase> {
        &self.position
    }

    /// Gets the canonical ordering key for this breakend.
    ///
    /// The key orders breakends first by contig, then by interbase position,
    /// and finally by orientation rank. It backs the canonical form of a
    /// paired adjacency.
    // Consumed by `Adjacency`'s canonical construction, added in the next task.
    #[allow(dead_code)]
    pub(crate) fn canonical_key(&self) -> (&str, Number, u8) {
        (
            self.contig.as_str(),
            self.position.get(),
            self.orientation.rank(),
        )
    }
}

impl FromStr for Breakend {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(VARIANT_SEPARATOR).collect::<Vec<_>>();
        let [contig, orientation, position] = parts.as_slice() else {
            return Err(ParseError::Format(s.to_string()));
        };

        let token = position
            .strip_suffix("(i)")
            .ok_or_else(|| ParseError::Qualifier {
                position: (*position).to_string(),
            })?;

        let position = token.parse::<InterbasePosition>()?;

        Breakend::try_new(*contig, orientation.parse::<Orientation>()?, position.get())
            .map_err(ParseError::Contig)
    }
}

impl std::fmt::Display for Breakend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{sep}{}{sep}{}(i)",
            self.contig,
            self.orientation,
            self.position.get(),
            sep = VARIANT_SEPARATOR,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_builds_and_accesses_a_breakend() {
        let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100).unwrap();
        assert_eq!(breakend.contig().as_str(), "seq0");
        assert_eq!(breakend.orientation(), Orientation::LowerFlank);
        assert_eq!(breakend.position().get(), 100);
    }

    #[test]
    fn it_parses_and_serializes() {
        let breakend = "seq0:>:100(i)".parse::<Breakend>().unwrap();
        assert_eq!(breakend.orientation(), Orientation::LowerFlank);
        assert_eq!(breakend.position().get(), 100);
        assert_eq!(breakend.to_string(), "seq0:>:100(i)");
    }

    #[test]
    fn it_rejects_a_missing_interbase_qualifier() {
        let err = "seq0:>:100".parse::<Breakend>().unwrap_err();
        assert!(matches!(err, ParseError::Qualifier { .. }));
    }
}
