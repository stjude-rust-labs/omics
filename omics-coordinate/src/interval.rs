//! An interval comprised of a start and end coordinate within a particular
//! coordinate system.

use std::str::FromStr;

use crate::Contig;
use crate::Coordinate;
use crate::Position;
use crate::Strand;
use crate::contig;
use crate::position;
use crate::strand;
use crate::system::System;

pub mod one;
pub mod zero;

/// A error related to parse an [`Interval`].
#[derive(Debug)]
pub enum ParseError {
    /// Could not parse an interval from an invalid format.
    InvalidFormat(String),

    /// An invalid [`Contig`] was attempted to be parsed.
    InvalidContig(contig::Error),

    /// An invalid [`Strand`] was attempted to be parsed.
    InvalidStrand(strand::Error),

    /// An invalid [`Position`] was attempted to be parsed.
    InvalidPosition(position::Error),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InvalidFormat(v) => write!(f, "invalid format: {v}"),
            ParseError::InvalidContig(err) => write!(f, "invalid contig: {err}"),
            ParseError::InvalidStrand(err) => write!(f, "invalid strand: {err}"),
            ParseError::InvalidPosition(err) => write!(f, "invalid position: {err}"),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to an [`Interval`].
#[derive(Debug)]
pub enum Error {
    /// A coordinate error.
    Coordinate(crate::Error),

    /// Attempted to create an invalid [`Coordinate`].
    InvalidCoordinate(crate::Error),

    /// Could not create an interval for coordinates on different contigs.
    MismatchedContigs(Contig, Contig),

    /// Could not create an interval for coordinates on different strands.
    MismatchedStrands(Strand, Strand),

    /// The start position cannot equal the end position, which would result in
    /// a zero-size interval.
    ZeroSizedInterval,

    /// The start position is greater than the end position for a positive
    /// stranded interval, which is not allowed.
    StartGreaterThanEndForPositiveStrand,

    /// The end position is greater than the start position for a negative
    /// stranded interval, which is not allowed.
    EndGreaterThanStartForNegativeStrand,

    /// Attempted to clamp an interval, but the two intervals had mismatching
    /// contigs.
    MismatchedContigDuringClamp(Contig, Contig),

    /// Attempted to clamp an interval, but the two intervals had mismatching
    /// strands.
    MismatchedStrandDuringClamp(Strand, Strand),

    /// Occurs when you try to create a singular interval, but moving the
    /// position forward one returns [`None`].
    CannotIncrementSingularInterval,

    /// A parse error.
    ParseError(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Coordinate(err) => write!(f, "coordinate error: {err}"),
            Error::InvalidCoordinate(err) => write!(f, "invalid coordinate: {}", err),
            Error::MismatchedContigs(a, b) => {
                write!(f, "mismatched contigs: {a}, {b}")
            }
            Error::MismatchedStrands(a, b) => {
                write!(f, "mismatched strands: {a}, {b}")
            }
            Error::ZeroSizedInterval => write!(
                f,
                "start position equals end position, which is a zero-sized interval"
            ),
            Error::StartGreaterThanEndForPositiveStrand => {
                write!(
                    f,
                    "start position cannot be greater than the end position for a positive \
                     stranded interval"
                )
            }
            Error::EndGreaterThanStartForNegativeStrand => {
                write!(
                    f,
                    "end position cannot be greater than the start position for a negative \
                     stranded interval"
                )
            }
            Error::MismatchedContigDuringClamp(a, b) => {
                write!(f, "mismatched contig while clamping: {a}, {b}")
            }
            Error::MismatchedStrandDuringClamp(a, b) => {
                write!(f, "mismatched strand while clamping: {a}, {b}")
            }
            Error::CannotIncrementSingularInterval => {
                write!(f, "cannot move singular interval position forward one")
            }
            Error::ParseError(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// An interval consisting of a start and end coordinate.
///
/// # Overview
///
/// Intervals specify a range of coordinates, including ranges of length one (a
/// single coordinate). Intervals work differently based on whether the interval
/// is 0-based or 1-based.
///
/// * **0-based** intervals are half-open, meaning the specified start
///   coordinate is _included_ in the range and the specified end coordinate is
///   _not included_ in the range. This range follows the mathematical notation
///   `[start, end)`.
/// * **1-based** intervals are fully-closed, meaning the specified start
///   coordinate _and_ the specified end coordinate are _both included_ in the
///   range. This range follows the mathematical notation `[start, end]`.
///
/// # Serialization and Parsing
///
/// Unlike many other tools in the community, including the UCSC genome browser,
/// this crate does not use different notations for expressing 0-based and
/// 1-based intervals. Instead, a single, common notation is used:
///
/// - Ranges are provided in the form `<contig>:<strand>:<start>-<end>` (e.g.,
///   `seq0:+:0-1000`). Ranges can also be specified on the negative strand
///   (e.g., `seq0:-:1000-0`).
/// - Singular positions are provided in the form `<contig>:<strand>:<position>`
///   (e.g., `seq0:+:1`). Singular positions are always parsed as an interval
///   containing only the position specified.
///
/// For a more in-depth discussion on this, please see [this section of the
/// docs](crate#intervals).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Interval<S: System>(Coordinate<S>, Coordinate<S>);

impl<S: System> Interval<S> {
    /// Attempts to create a new [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 0-based, positive-stranded interval
    ///
    /// let start = Coordinate::<Zero>::try_new("seq0", Strand::Positive, 1)?;
    /// let end = Coordinate::<Zero>::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// // 0-based, negative-stranded interval
    ///
    /// let start = Coordinate::<Zero>::try_new("seq0", Strand::Negative, 1000)?;
    /// let end = Coordinate::<Zero>::try_new("seq0", Strand::Negative, 1)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// // 1-based, positive-stranded interval
    ///
    /// let start = Coordinate::<One>::try_new("seq0", Strand::Positive, 1)?;
    /// let end = Coordinate::<One>::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// // 1-based, negative-stranded interval
    ///
    /// let start = Coordinate::<One>::try_new("seq0", Strand::Negative, 1000)?;
    /// let end = Coordinate::<One>::try_new("seq0", Strand::Negative, 1)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(start: Coordinate<S>, end: Coordinate<S>) -> Result<Interval<S>, Error>
    where
        Interval<S>: r#trait::Interval<S>,
    {
        <Self as r#trait::Interval<S>>::try_new(start, end)
    }

    /// Gets the start [`Coordinate`] of the [`Interval`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// let start = Coordinate::<Zero>::try_new("seq0", Strand::Positive, 0)?;
    /// let end = Coordinate::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// let interval = Interval::try_new(start.clone(), end)?;
    /// assert_eq!(interval.start(), &start);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn start(&self) -> &Coordinate<S> {
        &self.0
    }

    /// Gets the end [`Coordinate`] of the [`Interval`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// let start = Coordinate::<Zero>::try_new("seq0", Strand::Positive, 0)?;
    /// let end = Coordinate::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// let interval = Interval::try_new(start, end.clone())?;
    /// assert_eq!(interval.end(), &end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> &Coordinate<S> {
        &self.1
    }

    /// Gets the [`Contig`] for this [`Interval`] by reference.
    ///
    /// Note that the contig for the start and end coordinate are guaranteed to
    /// be the same by the checks done at interval creation time, so this just
    /// returns the contig of the start coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
    /// let end = "seq0:+:1000".parse::<Coordinate<Zero>>()?;
    ///
    /// let interval = Interval::try_new(start, end)?;
    /// assert_eq!(interval.contig().inner(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        self.0.contig()
    }

    /// Gets the [`Strand`] for this [`Interval`] by reference.
    ///
    /// Note that the strand for the start and end coordinate are guaranteed to
    /// be the same by the checks done at interval creation time, so this
    /// just returns the strand of the start coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// let start = Coordinate::<Zero>::try_new("seq0", Strand::Positive, 0)?;
    /// let end = Coordinate::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// let interval = Interval::try_new(start, end)?;
    /// assert_eq!(interval.strand(), &Strand::Positive);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> &Strand {
        self.0.strand()
    }

    /// Indicates whether a [`Coordinate`] falls within an [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 0-based, positive-stranded interval
    ///
    /// let interval = "seq0:+:0-1000".parse::<Interval<Zero>>()?;
    ///
    /// assert!(interval.contains(&"seq0:+:0".parse::<Coordinate<Zero>>()?));
    /// assert!(interval.contains(&"seq0:+:999".parse::<Coordinate<Zero>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:+:0".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:-:50".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:+:1000".parse::<Coordinate<Zero>>()?));
    ///
    /// // 0-based, negative-stranded interval
    ///
    /// let interval = "seq0:-:1000-0".parse::<Interval<Zero>>()?;
    ///
    /// assert!(interval.contains(&"seq0:-:1000".parse::<Coordinate<Zero>>()?));
    /// assert!(interval.contains(&"seq0:-:1".parse::<Coordinate<Zero>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:-:1000".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:+:1000".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:-:0".parse::<Coordinate<Zero>>()?));
    ///
    /// // 1-based, positive-stranded interval
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<One>>()?;
    ///
    /// assert!(interval.contains(&"seq0:+:1".parse::<Coordinate<One>>()?));
    /// assert!(interval.contains(&"seq0:+:999".parse::<Coordinate<One>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:+:1".parse::<Coordinate<One>>()?));
    /// assert!(!interval.contains(&"seq0:-:50".parse::<Coordinate<One>>()?));
    /// assert!(!interval.contains(&"seq0:+:1001".parse::<Coordinate<One>>()?));
    ///
    /// // 1-based, negative-stranded interval
    ///
    /// let interval = "seq0:-:1000-1".parse::<Interval<One>>()?;
    ///
    /// assert!(interval.contains(&"seq0:-:1000".parse::<Coordinate<One>>()?));
    /// assert!(interval.contains(&"seq0:-:1".parse::<Coordinate<One>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:-:1000".parse::<Coordinate<One>>()?));
    /// assert!(!interval.contains(&"seq0:+:1000".parse::<Coordinate<One>>()?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contains(&self, coordinate: &Coordinate<S>) -> bool
    where
        Interval<S>: r#trait::Interval<S>,
    {
        <Self as r#trait::Interval<S>>::contains(self, coordinate)
    }

    /// Gets the distance of the [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 0-based, positive-stranded interval
    ///
    /// let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
    /// let end = "seq0:+:1000".parse::<Coordinate<Zero>>()?;
    ///
    /// let interval = Interval::<Zero>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// // 0-based, negative-stranded interval
    ///
    /// let start = "seq0:-:1000".parse::<Coordinate<Zero>>()?;
    /// let end = "seq0:-:0".parse::<Coordinate<Zero>>()?;
    ///
    /// let interval = Interval::<Zero>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// // 1-based, positive-stranded interval
    ///
    /// let start = "seq0:+:1".parse::<Coordinate<One>>()?;
    /// let end = "seq0:+:1000".parse::<Coordinate<One>>()?;
    ///
    /// let interval = Interval::<One>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// // 1-based, negative-stranded interval
    ///
    /// let start = "seq0:-:1000".parse::<Coordinate<One>>()?;
    /// let end = "seq0:-:1".parse::<Coordinate<One>>()?;
    ///
    /// let interval = Interval::<One>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn distance(&self) -> usize
    where
        Interval<S>: r#trait::Interval<S>,
    {
        <Self as r#trait::Interval<S>>::len(self)
    }

    /// Consumes self to clamp an [`Interval`] based on the value of another
    /// interval.
    ///
    /// This method can, at most, shorten the consumed interval to be within the
    /// `interval` argument's bounds. If either end of the `interval` argument
    /// falls within the consumed interval's bounds, then the returned interval
    /// will be shortened to match the `interval` argument's start/end position.
    /// If the `interval` argument completely encapsulates the consumed
    /// interval, then the original interval will be returned unmodified.
    ///
    /// Note that this function returns [`Result<Interval, Error>`] instead of
    /// an [`Option<Interval>`] because its expected that, if you are clamping
    /// an interval, you're probably doing so with the expectation that the
    /// intervals are on the same contig and strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::interval::Error;
    /// use omics_coordinate::system::One;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<One>>()?;
    /// let clamp = "seq0:+:250-2000".parse::<Interval<One>>()?;
    /// let result = interval.clamp(&clamp)?;
    /// assert_eq!(&result, &"seq0:+:250-1000".parse::<Interval<One>>()?);
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:-:2000-1000".parse::<Interval<One>>()?;
    /// let clamp = "seq0:-:3000-1250".parse::<Interval<One>>()?;
    /// let result = interval.clamp(&clamp)?;
    /// assert_eq!(&result, &"seq0:-:2000-1250".parse::<Interval<One>>()?);
    ///
    /// // Differing contigs
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<One>>()?;
    /// let clamp = "seq1:+:250-2000".parse::<Interval<One>>()?;
    /// let result = interval.clamp(&clamp);
    /// assert!(matches!(
    ///     result,
    ///     Err(Error::MismatchedContigDuringClamp(_, _))
    /// ));
    ///
    /// // Differing strands
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<One>>()?;
    /// let clamp = "seq0:-:2000-250".parse::<Interval<One>>()?;
    /// let result = interval.clamp(&clamp);
    /// assert!(matches!(
    ///     result,
    ///     Err(Error::MismatchedStrandDuringClamp(_, _))
    /// ));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn clamp(&self, b: &Interval<S>) -> Result<Interval<S>, Error>
    where
        Interval<S>: r#trait::Interval<S>,
        Position<S>: position::r#trait::Position<S>,
    {
        if self.contig() != b.contig() {
            return Err(Error::MismatchedContigDuringClamp(
                self.contig().clone(),
                b.contig().clone(),
            ));
        }

        if self.strand() != b.strand() {
            return Err(Error::MismatchedStrandDuringClamp(
                self.strand().clone(),
                b.strand().clone(),
            ));
        }

        let start_position = match self.strand() {
            Strand::Positive => std::cmp::max(self.start().position(), b.start().position()),
            Strand::Negative => std::cmp::min(self.start().position(), b.start().position()),
        };

        let end_position = match self.strand() {
            Strand::Positive => std::cmp::min(self.end().position(), b.end().position()),
            Strand::Negative => std::cmp::max(self.end().position(), b.end().position()),
        };

        let start = Coordinate::<S>::try_new(
            self.contig().clone(),
            self.strand().clone(),
            start_position.clone(),
        )
        .map_err(Error::InvalidCoordinate)?;

        let end = Coordinate::<S>::try_new(
            self.contig().clone(),
            self.strand().clone(),
            end_position.clone(),
        )
        .map_err(Error::InvalidCoordinate)?;

        Self::try_new(start, end)
    }

    /// Gets the offset of the provided [`Coordinate`] from the starting
    /// position of the interval.
    ///
    /// - If the provided coordinate _does not_ fall within the interval, then
    ///   [`None`] is returned.
    /// - If the provided coordinate _does_ fall within the interval, then the
    ///   distance from the starting position of the interval is returned
    ///   wrapped in [`Some<usize>`].
    ///
    /// Note that for [`Strand::Positive`] intervals, the magnitude should be
    /// interpreted as a _positive_ magnitude, whereas a [`Strand::Negative`]
    /// offset will be a _negative_ magnitude (though it will not be reflected
    /// in the size of the returned result).
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:+:0-1000".parse::<Interval<Zero>>()?;
    ///
    /// let target = "seq0:+:5".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, Some(5));
    ///
    /// let target = "seq0:+:999".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, Some(999));
    ///
    /// let target = "seq1:+:5".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, None);
    ///
    /// let target = "seq0:-:5".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, None);
    ///
    /// let target = "seq0:+:1000".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, None);
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:-:1000-0".parse::<Interval<Zero>>()?;
    ///
    /// let target = "seq0:-:995".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, Some(5));
    ///
    /// let target = "seq0:-:1".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, Some(999));
    ///
    /// let target = "seq1:-:955".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, None);
    ///
    /// let target = "seq0:+:995".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, None);
    ///
    /// let target = "seq:-:0".parse::<Coordinate<Zero>>()?;
    /// let offset = interval.offset(&target);
    /// assert_eq!(offset, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn offset(self: &Interval<S>, coordinate: &Coordinate<S>) -> Option<usize>
    where
        Interval<S>: r#trait::Interval<S>,
        Position<S>: position::r#trait::Position<S>,
    {
        if !self.contains(coordinate) {
            return None;
        }

        match self.strand() {
            Strand::Positive => coordinate
                .position()
                .distance_unchecked(self.start().position()),
            Strand::Negative => self
                .start()
                .position()
                .distance_unchecked(coordinate.position()),
        }
    }

    /// Attempts to translate an offset from the start of the interval to a
    /// coordinate.
    ///
    /// - If an error occurs, an [`Err(Error)`] is returned.
    /// - If the translated coordinate _does not_ fall within the interval, then
    ///   [`Ok(None)`] is returned.
    /// - If the translated coordinate _does_ fall within the interval, then the
    ///   coordinate is returned wrapped in [`Ok(Some(Position))`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:+:0-1000".parse::<Interval<Zero>>()?;
    ///
    /// assert_eq!(
    ///     interval.translate(5)?,
    ///     Some("seq0:+:5".parse::<Coordinate<Zero>>()?)
    /// );
    /// assert_eq!(
    ///     interval.translate(999)?,
    ///     Some("seq0:+:999".parse::<Coordinate<Zero>>()?)
    /// );
    /// assert_eq!(interval.translate(1000)?, None);
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:-:999-[".parse::<Interval<Zero>>()?;
    ///
    /// assert_eq!(
    ///     interval.translate(5)?,
    ///     Some("seq0:-:994".parse::<Coordinate<Zero>>()?)
    /// );
    /// assert_eq!(
    ///     interval.translate(999)?,
    ///     Some("seq0:-:0".parse::<Coordinate<Zero>>()?)
    /// );
    /// assert_eq!(interval.translate(1000)?, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn translate(&self, offset: usize) -> Result<Option<Coordinate<S>>, Error>
    where
        Interval<S>: r#trait::Interval<S>,
        Position<S>: position::r#trait::Position<S>,
    {
        let coordinate = self
            .start()
            .clone()
            .move_forward(offset)
            .map_err(Error::InvalidCoordinate)?;

        match coordinate {
            Some(c) => match self.contains(&c) {
                true => Ok(Some(c)),
                false => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Consumes `self` to return the inner [`Coordinate<S>`] that make up this
    /// [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::Zero;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:+:0-1000".parse::<Interval<Zero>>()?;
    /// let (start, end) = interval.into_coordinates();
    ///
    /// assert_eq!(start.contig().inner(), "seq0");
    /// assert_eq!(start.strand(), &Strand::Positive);
    /// assert_eq!(start.position().inner(), &Value::Usize(0));
    /// assert_eq!(end.contig().inner(), "seq0");
    /// assert_eq!(end.strand(), &Strand::Positive);
    /// assert_eq!(end.position().inner(), &Value::Usize(1000));
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:-:1000-0".parse::<Interval<Zero>>()?;
    /// let (start, end) = interval.into_coordinates();
    ///
    /// assert_eq!(start.contig().inner(), "seq0");
    /// assert_eq!(start.strand(), &Strand::Negative);
    /// assert_eq!(start.position().inner(), &Value::Usize(1000));
    /// assert_eq!(end.contig().inner(), "seq0");
    /// assert_eq!(end.strand(), &Strand::Negative);
    /// assert_eq!(end.position().inner(), &Value::Usize(0));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_coordinates(self) -> (Coordinate<S>, Coordinate<S>) {
        (self.0, self.1)
    }
}

/// Attempts to create a new [`Interval<S>`].
fn try_new<S: System>(start: Coordinate<S>, end: Coordinate<S>) -> Result<Interval<S>, Error>
where
    Interval<S>: r#trait::Interval<S>,
    Coordinate<S>: crate::r#trait::Coordinate<S>,
    Position<S>: position::r#trait::Position<S>,
{
    if start.contig() != end.contig() {
        return Err(Error::MismatchedContigs(
            start.contig().clone(),
            end.contig().clone(),
        ));
    }

    if start.strand() != end.strand() {
        return Err(Error::MismatchedStrands(
            start.strand().clone(),
            end.strand().clone(),
        ));
    }

    match start.strand() {
        Strand::Positive => {
            if start.position() > end.position() {
                return Err(Error::StartGreaterThanEndForPositiveStrand);
            }
        }
        Strand::Negative => {
            if end.position() > start.position() {
                return Err(Error::EndGreaterThanStartForNegativeStrand);
            }
        }
    }

    Ok(Interval(start, end))
}

impl<S: System> std::fmt::Display for Interval<S>
where
    Self: r#trait::Interval<S>,
    Coordinate<S>: crate::r#trait::Coordinate<S>,
    Position<S>: position::r#trait::Position<S>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}-{}",
            self.contig(),
            self.strand(),
            self.start().position(),
            self.end().position(),
        )
    }
}

impl<S: System> FromStr for Interval<S>
where
    Interval<S>: crate::interval::r#trait::Interval<S>,
    Coordinate<S>: crate::r#trait::Coordinate<S>,
    Position<S>: position::r#trait::Position<S>,
{
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(':').collect::<Vec<_>>();

        if parts.len() != 3 {
            return Err(Error::ParseError(ParseError::InvalidFormat(s.to_string())));
        }

        let mut parts = parts.iter();

        // SAFETY: we checked that there are three parts above. Given that we
        // haven't pulled anything from the iterator, we can always safely
        // unwrap this.
        let contig = parts
            .next()
            .unwrap()
            .parse::<Contig>()
            .map_err(|err| Error::ParseError(ParseError::InvalidContig(err)))?;

        // SAFETY: we checked that there are three parts above. Given that we
        // have only pulled one item from the iterator, we can always safely
        // unwrap this.
        let strand = parts
            .next()
            .unwrap()
            .parse::<Strand>()
            .map_err(|err| Error::ParseError(ParseError::InvalidStrand(err)))?;

        // SAFETY: we checked that there are three parts above. Given that we
        // have only pulled two items from the iterator, we can always safely
        // unwrap this.
        let position_parts = parts.next().unwrap().split('-').collect::<Vec<_>>();

        match position_parts.len() {
            1 => {
                // SAFETY: we just ensured that one part exists, so the direct
                // indexing of the slice will never fail.
                let position = position_parts[0]
                    .parse::<Position<S>>()
                    .map_err(|err| Error::ParseError(ParseError::InvalidPosition(err)))?;

                let start = Coordinate::try_new(contig, strand, position)
                    .map_err(Error::InvalidCoordinate)?;

                let end = start
                    .clone()
                    .move_forward(1)
                    .map_err(Error::InvalidCoordinate)?;

                match end {
                    Some(end) => Interval::try_new(start, end),
                    None => Err(Error::CannotIncrementSingularInterval),
                }
            }
            2 => {
                // SAFETY: we just ensured that two parts exist, so the direct
                // indexing of the slice for both index zero and one will never
                // fail.
                let start_position = position_parts[0]
                    .parse::<Position<S>>()
                    .map_err(|err| Error::ParseError(ParseError::InvalidPosition(err)))?;
                let end_position = position_parts[1]
                    .parse::<Position<S>>()
                    .map_err(|err| Error::ParseError(ParseError::InvalidPosition(err)))?;

                Interval::try_new(
                    Coordinate::try_new(contig.clone(), strand.clone(), start_position)
                        .map_err(Error::InvalidCoordinate)?,
                    Coordinate::try_new(contig, strand, end_position)
                        .map_err(Error::InvalidCoordinate)?,
                )
            }
            _ => Err(Error::ParseError(ParseError::InvalidFormat(s.to_string()))),
        }
    }
}

/// Traits related to an interval.
pub mod r#trait {
    use super::*;

    /// Requirements to be an interval.
    #[allow(clippy::len_without_is_empty)]
    pub trait Interval<S: System> {
        /// Attempts to create a new [`Interval<S>`].
        fn try_new(start: Coordinate<S>, end: Coordinate<S>) -> Result<super::Interval<S>, Error>;

        /// Returns whether the [`Coordinate<S>`] is contained within this
        /// [`Interval`].
        fn contains(&self, coordinate: &Coordinate<S>) -> bool;

        /// Gets the length of the [`Interval`].
        fn len(&self) -> usize;

        /// Complements the interval by swapping the strand and the start and
        /// end coordinates.
        fn complement(self) -> Result<Option<super::Interval<S>>, Error>;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::Zero;

    #[test]
    fn test_it_creates_a_valid_interval() -> Result<(), Box<dyn std::error::Error>> {
        let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
        let end = "seq0:+:10".parse::<Coordinate<Zero>>()?;

        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.distance(), 10);

        Ok(())
    }

    #[test]
    fn test_it_errors_when_contigs_differ() -> Result<(), Box<dyn std::error::Error>> {
        let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
        let end = "seq1:+:10".parse::<Coordinate<Zero>>()?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::MismatchedContigs(_, _)));
        assert_eq!(err.to_string(), "mismatched contigs: seq0, seq1");

        Ok(())
    }

    #[test]
    fn test_it_errors_when_strands_differ() -> Result<(), Box<dyn std::error::Error>> {
        let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
        let end = "seq0:-:10".parse::<Coordinate<Zero>>()?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::MismatchedStrands(_, _)));
        assert_eq!(err.to_string(), "mismatched strands: +, -");

        Ok(())
    }

    #[test]
    fn test_it_does_not_create_a_zero_sized_interval() -> Result<(), Box<dyn std::error::Error>> {
        let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
        let end = "seq0:+:0".parse::<Coordinate<Zero>>()?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::ZeroSizedInterval));

        Ok(())
    }

    #[test]
    fn test_it_does_not_allow_start_to_be_greater_than_end_for_positive_stranded_interval()
    -> Result<(), Box<dyn std::error::Error>> {
        let start = "seq0:+:10".parse::<Coordinate<Zero>>()?;
        let end = "seq0:+:0".parse::<Coordinate<Zero>>()?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::StartGreaterThanEndForPositiveStrand));

        Ok(())
    }

    #[test]
    fn test_it_does_not_allow_end_to_be_greater_than_start_for_negative_stranded_interval()
    -> Result<(), Box<dyn std::error::Error>> {
        let start = "seq0:-:0".parse::<Coordinate<Zero>>()?;
        let end = "seq0:-:10".parse::<Coordinate<Zero>>()?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::EndGreaterThanStartForNegativeStrand));

        Ok(())
    }

    #[test]
    fn test_it_clamps_correctly_for_positive_stranded_intervals()
    -> Result<(), Box<dyn std::error::Error>> {
        let interval = "seq0:+:1000-2000".parse::<Interval<Zero>>()?;

        assert!(matches!(
            interval
                .clone()
                .clamp(&"seq1:+:0-3000".parse::<Interval<Zero>>()?),
            Err(Error::MismatchedContigDuringClamp(_, _))
        ));

        assert!(matches!(
            interval
                .clone()
                .clamp(&"seq0:-:3000-0".parse::<Interval<Zero>>()?),
            Err(Error::MismatchedStrandDuringClamp(_, _))
        ));

        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:+:0-3000".parse::<Interval<Zero>>()?)?,
            "seq0:+:1000-2000".parse::<Interval<Zero>>()?
        );

        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:+:1250-3000".parse::<Interval<Zero>>()?)?,
            "seq0:+:1250-2000".parse::<Interval<Zero>>()?
        );

        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:+:0-1750".parse::<Interval<Zero>>()?)?,
            "seq0:+:1000-1750".parse::<Interval<Zero>>()?
        );

        assert_eq!(
            interval.clamp(&"seq0:+:1250-1750".parse::<Interval<Zero>>()?)?,
            "seq0:+:1250-1750".parse::<Interval<Zero>>()?
        );

        Ok(())
    }

    #[test]
    fn test_it_clamps_correctly_for_negative_stranded_intervals()
    -> Result<(), Box<dyn std::error::Error>> {
        let interval = "seq0:-:2000-1000".parse::<Interval<Zero>>()?;

        assert!(matches!(
            interval
                .clone()
                .clamp(&"seq1:-:3000-0".parse::<Interval<Zero>>()?),
            Err(Error::MismatchedContigDuringClamp(_, _))
        ));

        assert!(matches!(
            interval
                .clone()
                .clamp(&"seq0:+:0-3000".parse::<Interval<Zero>>()?),
            Err(Error::MismatchedStrandDuringClamp(_, _))
        ));

        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:-:3000-0".parse::<Interval<Zero>>()?)?,
            "seq0:-:2000-1000".parse::<Interval<Zero>>()?
        );

        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:-:3000-1250".parse::<Interval<Zero>>()?)?,
            "seq0:-:2000-1250".parse::<Interval<Zero>>()?
        );

        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:-:1750-0".parse::<Interval<Zero>>()?)?,
            "seq0:-:1750-1000".parse::<Interval<Zero>>()?
        );

        assert_eq!(
            interval.clamp(&"seq0:-:1750-1250".parse::<Interval<Zero>>()?)?,
            "seq0:-:1750-1250".parse::<Interval<Zero>>()?
        );

        Ok(())
    }

    #[test]
    fn test_parsing_intervals_works_for_valid_single_position()
    -> Result<(), Box<dyn std::error::Error>> {
        let interval = "seq0:+:1".parse::<Interval<Zero>>()?;
        let start = "seq0:+:1".parse::<Coordinate<Zero>>()?;
        let end = "seq0:+:2".parse::<Coordinate<Zero>>()?;

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);

        Ok(())
    }

    #[test]
    fn test_parsing_intervals_works_for_valid_position_range()
    -> Result<(), Box<dyn std::error::Error>> {
        // Testing positive stranded interval
        let interval = "seq0:+:1-1000".parse::<Interval<Zero>>()?;
        let start = "seq0:+:1".parse::<Coordinate<Zero>>()?;
        let end = "seq0:+:1000".parse::<Coordinate<Zero>>()?;

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);

        // Testing negative stranded interval
        let interval = "seq0:-:1000-1".parse::<Interval<Zero>>()?;
        let start = "seq0:-:1000".parse::<Coordinate<Zero>>()?;
        let end = "seq0:-:1".parse::<Coordinate<Zero>>()?;

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);

        // Testing running up until negative bound
        let interval = "seq0:-:1000-[".parse::<Interval<Zero>>()?;
        let start = "seq0:-:1000".parse::<Coordinate<Zero>>()?;
        let end = Coordinate::<Zero>::lower_bound("seq0");

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);

        Ok(())
    }

    #[test]
    fn test_various_invalid_intervals() -> Result<(), Box<dyn std::error::Error>> {
        let err = "1".parse::<Interval<Zero>>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "1-1000".parse::<Interval<Zero>>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "seq0:".parse::<Interval<Zero>>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "seq0:1-".parse::<Interval<Zero>>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "seq0:1-10000:".parse::<Interval<Zero>>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        Ok(())
    }

    #[test]
    fn test_interval_to_string() -> Result<(), Box<dyn std::error::Error>> {
        // Positive-stranded interval
        let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
        let end = "seq0:+:10".parse::<Coordinate<Zero>>()?;
        let interval = Interval::try_new(start, end)?;

        assert_eq!(interval.to_string(), "seq0:+:0-10");

        // Negative-stranded interval
        let start = "seq0:-:10".parse::<Coordinate<Zero>>()?;
        let end = "seq0:-:0".parse::<Coordinate<Zero>>()?;
        let interval = Interval::try_new(start, end)?;

        assert_eq!(interval.to_string(), "seq0:-:10-0");

        // Negative-stranded interval with negative bound
        let start = "seq0:-:10".parse::<Coordinate<Zero>>()?;
        let end = "seq0:-:[".parse::<Coordinate<Zero>>()?;
        let interval = Interval::try_new(start, end)?;

        assert_eq!(interval.to_string(), "seq0:-:10-[");

        Ok(())
    }
}
