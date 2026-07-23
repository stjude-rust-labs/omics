//! Copy-number variants.

use std::cmp::Ordering;
use std::fmt;
use std::num::NonZeroU32;
use std::str::FromStr;

use omics_coordinate::Contig;
use omics_coordinate::Position;
use omics_coordinate::contig;
use omics_coordinate::position::Number;
use omics_coordinate::position::interbase::Position as InterbasePosition;
use omics_coordinate::system::Interbase;
use omics_core::VARIANT_SEPARATOR;
use thiserror::Error;

/// The copy-number change relative to the reference ploidy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Change {
    /// The observed copy number is below the reference ploidy.
    Loss,
    /// The observed copy number matches the reference ploidy.
    Reference,
    /// The observed copy number is above the reference ploidy.
    Gain,
}

/// An absolute copy count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Count(u32);

impl Count {
    /// Creates a copy count from a raw integer value.
    pub const fn new(count: u32) -> Self {
        Self(count)
    }

    /// Gets the raw copy count.
    pub const fn get(self) -> u32 {
        self.0
    }

    /// Attempts to construct a count from a base-2 logarithmic ratio.
    pub fn try_from_log2(value: f64, ploidy: Ploidy) -> Result<Self, LogarithmicError> {
        Self::try_from_logarithmic(value, ploidy, LogarithmicBase::Two)
    }

    /// Attempts to construct a count from a base-10 logarithmic ratio.
    pub fn try_from_log10(value: f64, ploidy: Ploidy) -> Result<Self, LogarithmicError> {
        Self::try_from_logarithmic(value, ploidy, LogarithmicBase::Ten)
    }

    /// Gets the base-2 logarithmic ratio relative to the ploidy baseline.
    pub fn log2(self, ploidy: Ploidy) -> f64 {
        let ratio = f64::from(self.get()) / f64::from(ploidy.get());
        ratio.log2()
    }

    /// Gets the base-10 logarithmic ratio relative to the ploidy baseline.
    pub fn log10(self, ploidy: Ploidy) -> f64 {
        let ratio = f64::from(self.get()) / f64::from(ploidy.get());
        ratio.log10()
    }

    /// Converts a logarithmic ratio at the given base into an absolute count.
    fn try_from_logarithmic(
        value: f64,
        ploidy: Ploidy,
        base: LogarithmicBase,
    ) -> Result<Self, LogarithmicError> {
        if value.is_nan() {
            return Err(LogarithmicError::NotANumber);
        }

        if value == f64::NEG_INFINITY {
            return Ok(Self::new(0));
        }

        if value == f64::INFINITY {
            return Err(LogarithmicError::PositiveInfinity);
        }

        let copies = f64::from(ploidy.get()) * base.powf(value);

        if copies == 0.0 {
            return Err(LogarithmicError::Underflow);
        }

        let rounded = copies.round();
        if rounded == 0.0 {
            return Err(LogarithmicError::Underflow);
        }

        let tolerance = 8.0 * f64::EPSILON * copies.abs().max(1.0);
        if (copies - rounded).abs() > tolerance {
            return Err(LogarithmicError::NonIntegral);
        }

        if rounded > f64::from(u32::MAX) {
            return Err(LogarithmicError::Overflow);
        }

        Ok(Self::new(rounded as u32))
    }
}

/// A positive baseline ploidy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Ploidy(NonZeroU32);

impl Ploidy {
    /// The diploid baseline.
    // SAFETY: `2` is nonzero.
    pub const DIPLOID: Self = Self(unsafe { NonZeroU32::new_unchecked(2) });
    /// The haploid baseline.
    pub const HAPLOID: Self = Self(NonZeroU32::MIN);

    /// Attempts to create a ploidy from a raw integer value.
    pub const fn try_new(value: u32) -> Result<Self, PloidyError> {
        match NonZeroU32::new(value) {
            Some(value) => Ok(Self(value)),
            None => Err(PloidyError::Zero),
        }
    }

    /// Gets the raw ploidy value.
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

/// An error related to ploidy construction.
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum PloidyError {
    /// The ploidy was zero.
    #[error("ploidy cannot be zero")]
    Zero,
}

impl TryFrom<u32> for Ploidy {
    type Error = PloidyError;

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        Self::try_new(value)
    }
}

/// An error related to logarithmic conversion.
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogarithmicError {
    /// The logarithmic value was not a number.
    #[error("logarithmic copy-number value cannot be `nan`")]
    NotANumber,

    /// The logarithmic value was positive infinity.
    #[error("logarithmic copy-number value cannot be positive infinity")]
    PositiveInfinity,

    /// The logarithmic value underflowed below one copy.
    #[error("logarithmic copy-number value underflowed below one copy")]
    Underflow,

    /// The logarithmic value overflowed the supported count range.
    #[error("logarithmic copy-number value overflowed the `u32` count range")]
    Overflow,

    /// The logarithmic value did not map to an integral count.
    #[error("logarithmic copy-number value did not map to an integral count")]
    NonIntegral,
}

/// A logarithmic base used for copy-number conversion.
#[derive(Debug, Clone, Copy)]
enum LogarithmicBase {
    /// Base 2 ratios.
    Two,
    /// Base 10 ratios.
    Ten,
}

impl LogarithmicBase {
    /// Raises the selected base to the provided exponent.
    fn powf(self, value: f64) -> f64 {
        match self {
            Self::Two => 2.0_f64.powf(value),
            Self::Ten => 10.0_f64.powf(value),
        }
    }
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
    /// The variant did not match the canonical copy-number form.
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

    /// The ploidy failed to parse as a positive integer.
    #[error("invalid copy-number ploidy: `{ploidy}`")]
    Ploidy {
        /// The offending ploidy token.
        ploidy: String,
    },

    /// The ploidy was zero.
    #[error("copy-number ploidy cannot be zero")]
    ZeroPloidy,

    /// The region has no span.
    #[error("copy-number region cannot be empty")]
    EmptyRegion,

    /// The region end precedes the start.
    #[error("copy-number region end must be greater than start")]
    ReversedRegion,
}

impl From<Error> for ParseError {
    fn from(value: Error) -> Self {
        match value {
            Error::Contig(err) => Self::Contig(err),
            Error::EmptyRegion => Self::EmptyRegion,
            Error::ReversedRegion => Self::ReversedRegion,
        }
    }
}

/// A strandless, half-open copy-number variant.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Variant {
    /// The contig containing the affected interval.
    contig: Contig,
    /// The half-open start position.
    start: Position<Interbase>,
    /// The half-open end position.
    end: Position<Interbase>,
    /// The observed copy count over the interval.
    count: Count,
    /// The reference ploidy for the interval.
    ploidy: Ploidy,
}

impl Variant {
    /// Attempts to create a new copy-number variant.
    pub fn try_new(
        contig: impl TryInto<Contig, Error = contig::Error>,
        start: Number,
        end: Number,
        copies: u32,
        ploidy: Ploidy,
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
            count: Count::new(copies),
            ploidy,
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
    pub fn count(&self) -> Count {
        self.count
    }

    /// Gets the reference ploidy.
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    /// Classifies the variant against its stored reference ploidy.
    pub fn change(&self) -> Change {
        match self.count.get().cmp(&self.ploidy.get()) {
            Ordering::Less => Change::Loss,
            Ordering::Equal => Change::Reference,
            Ordering::Greater => Change::Gain,
        }
    }

    /// Gets the base-2 logarithmic ratio relative to the reference ploidy.
    pub fn log2(&self) -> f64 {
        self.count.log2(self.ploidy)
    }

    /// Gets the base-10 logarithmic ratio relative to the reference ploidy.
    pub fn log10(&self) -> f64 {
        self.count.log10(self.ploidy)
    }
}

impl FromStr for Variant {
    type Err = ParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let Some((contig_and_range, copies_and_ploidy)) = value.rsplit_once(VARIANT_SEPARATOR)
        else {
            return Err(ParseError::Format(value.to_string()));
        };
        let Some((contig, range)) = contig_and_range.rsplit_once(VARIANT_SEPARATOR) else {
            return Err(ParseError::Format(value.to_string()));
        };
        let Some((copies, ploidy)) = copies_and_ploidy.split_once('/') else {
            return Err(ParseError::Format(value.to_string()));
        };

        if copies.is_empty() || ploidy.is_empty() || ploidy.contains('/') {
            return Err(ParseError::Format(value.to_string()));
        }

        let range = range
            .strip_suffix("(i)")
            .ok_or_else(|| ParseError::Qualifier {
                region: range.to_string(),
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
            copies: copies.to_string(),
        })?;
        let ploidy = ploidy.parse::<u32>().map_err(|_| ParseError::Ploidy {
            ploidy: ploidy.to_string(),
        })?;
        let ploidy = Ploidy::try_new(ploidy).map_err(|_| ParseError::ZeroPloidy)?;

        Ok(Variant::try_new(
            contig,
            start.get(),
            end.get(),
            copies,
            ploidy,
        )?)
    }
}

impl fmt::Display for Variant {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{contig}{sep}{start}-{end}(i){sep}{copies}/{ploidy}",
            contig = self.contig,
            start = self.start.get(),
            end = self.end.get(),
            copies = self.count.get(),
            ploidy = self.ploidy.get(),
            sep = VARIANT_SEPARATOR,
        )
    }
}
