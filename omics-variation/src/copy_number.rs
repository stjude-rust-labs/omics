//! Copy-number variants.

use std::cmp::Ordering;
use std::fmt;
use std::str::FromStr;

use omics_coordinate::Contig;
use omics_coordinate::Position;
use omics_coordinate::contig;
use omics_coordinate::position::Number;
use omics_coordinate::position::interbase::Position as InterbasePosition;
use omics_coordinate::system::Interbase;
use omics_core::VARIANT_SEPARATOR;
use thiserror::Error;

/// The copy-number change relative to a baseline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Change {
    /// The observed copy number is below baseline.
    Loss,
    /// The observed copy number matches baseline.
    Baseline,
    /// The observed copy number is above baseline.
    Gain,
}

/// An error related to a copy-number variant.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// The contig failed to parse.
    #[error(transparent)]
    Contig(#[from] contig::Error),

    /// The region has no span.
    #[error("copy-number region cannot be empty")]
    EmptyRegion,

    /// The region end precedes the start.
    #[error("copy-number region end must be greater than start")]
    ReversedRegion,
}

/// An error related to parsing a copy-number variant.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// The variant did not have the three colon-separated fields.
    #[error("invalid copy-number format: `{0}`")]
    Format(String),

    /// The position range did not end with the interbase qualifier.
    #[error("copy-number range `{region}` must end with `(i)` for an interbase coordinate")]
    Qualifier {
        /// The offending range token.
        region: String,
    },

    /// The position range was malformed.
    #[error("invalid copy-number range: `{range}`")]
    Range {
        /// The offending range token.
        range: String,
    },

    /// The contig failed to parse.
    #[error(transparent)]
    Contig(#[from] contig::Error),

    /// The start or end position failed to parse.
    #[error(transparent)]
    Position(#[from] omics_coordinate::position::Error),

    /// The copy count failed to parse.
    #[error("invalid copy-number count: `{copies}`")]
    Copies {
        /// The offending count token.
        copies: String,
    },

    /// The region has no span.
    #[error("copy-number region cannot be empty")]
    EmptyRegion,

    /// The region end precedes the start.
    #[error("copy-number region end must be greater than start")]
    ReversedRegion,

    /// Construction failed after parsing.
    #[error(transparent)]
    Construct(#[from] Error),
}

/// A strandless, half-open copy-number variant.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Variant {
    contig: Contig,
    start: Position<Interbase>,
    end: Position<Interbase>,
    copies: u32,
}

impl Variant {
    /// Attempts to create a new copy-number variant.
    pub fn try_new(
        contig: impl TryInto<Contig, Error = contig::Error>,
        start: Number,
        end: Number,
        copies: u32,
    ) -> Result<Self, Error> {
        let contig = contig.try_into()?;

        match start.cmp(&end) {
            Ordering::Equal => return Err(Error::EmptyRegion),
            Ordering::Greater => return Err(Error::ReversedRegion),
            Ordering::Less => {}
        }

        Ok(Self {
            contig,
            start: InterbasePosition::new(start),
            end: InterbasePosition::new(end),
            copies,
        })
    }

    /// Gets the contig this variant sits on.
    pub fn contig(&self) -> &Contig {
        &self.contig
    }

    /// Gets the start position.
    pub fn start(&self) -> &Position<Interbase> {
        &self.start
    }

    /// Gets the end position.
    pub fn end(&self) -> &Position<Interbase> {
        &self.end
    }

    /// Gets the number of copies.
    pub fn copies(&self) -> u32 {
        self.copies
    }

    /// Classifies the variant against a baseline copy number.
    pub fn change(&self, baseline: u32) -> Change {
        match self.copies.cmp(&baseline) {
            Ordering::Less => Change::Loss,
            Ordering::Equal => Change::Baseline,
            Ordering::Greater => Change::Gain,
        }
    }
}

impl FromStr for Variant {
    type Err = ParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let parts = value.split(VARIANT_SEPARATOR).collect::<Vec<_>>();
        let [contig, range, copies] = parts.as_slice() else {
            return Err(ParseError::Format(value.to_string()));
        };

        let range = range
            .strip_suffix("(i)")
            .ok_or_else(|| ParseError::Qualifier {
                region: (*range).to_string(),
            })?;

        let (start, end) = range.split_once('-').ok_or_else(|| ParseError::Range {
            range: range.to_string(),
        })?;

        if start.is_empty() || end.is_empty() || end.contains('-') {
            return Err(ParseError::Range {
                range: range.to_string(),
            });
        }

        let start = start.parse::<InterbasePosition>()?;
        let end = end.parse::<InterbasePosition>()?;
        let copies = copies.parse::<u32>().map_err(|_| ParseError::Copies {
            copies: (*copies).to_string(),
        })?;

        Variant::try_new(*contig, start.get(), end.get(), copies).map_err(ParseError::Construct)
    }
}

impl fmt::Display for Variant {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{contig}{sep}{start}-{end}(i){sep}{copies}",
            contig = self.contig,
            start = self.start.get(),
            end = self.end.get(),
            copies = self.copies,
            sep = VARIANT_SEPARATOR,
        )
    }
}
