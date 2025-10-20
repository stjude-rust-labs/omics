//! Intervals.

use std::cmp::max;
use std::cmp::min;

use thiserror::Error;

use crate::Contig;
use crate::Position;
use crate::Strand;
use crate::System;
use crate::coordinate;
use crate::coordinate::Coordinate;
use crate::position;
use crate::position::Number;
use crate::strand;
use crate::system::Base;

pub mod base;
pub mod interbase;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// An error that occurs during clamping.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ClampError {
    /// A mismatched contig error.
    ///
    /// This error occurs when one attempts to clamp an interval with another
    /// interval that is not located on the same contig.
    #[error("mismatched contigs: `{original}` and `{operand}`")]
    MismatchedContigs {
        /// The contig of the interval being clamped.
        original: Contig,

        /// The contig of the interval doing the clamping.
        operand: Contig,
    },

    /// A mismatched strand error.
    ///
    /// This error occurs when one attempts to clamp an interval with another
    /// interval that is not located on the same strand.
    #[error("mismatched strand: `{original}` and `{operand}`")]
    MismatchedStrand {
        /// The strand of the interval being clamped.
        original: Strand,

        /// The strand of the interval doing the clamping.
        operand: Strand,
    },
}

/// A [`Result`](std::result::Result) with a [`ClampError`].
pub type ClampResult<T> = std::result::Result<T, ClampError>;

/// An error related to the creation of a nonsensical interval.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum NonsensicalError {
    /// A mismatched contig error.
    ///
    /// This error occurs when one attempts to clamp an interval with another
    /// interval that is not located on the same contig.
    #[error("mismatched contigs for coordinates: `{start}` and `{end}`")]
    MismatchedContigs {
        /// The contig of the interval being clamped.
        start: Contig,

        /// The contig of the interval doing the clamping.
        end: Contig,
    },

    /// A mismatched strand error.
    ///
    /// This error occurs when one attempts to clamp an interval with another
    /// interval that is not located on the same strand.
    #[error("mismatched strands for coordinates: `{start}` and `{end}`")]
    MismatchedStrands {
        /// The strand of the interval being clamped.
        start: Strand,

        /// The strand of the interval doing the clamping.
        end: Strand,
    },

    /// A negative sized interval.
    ///
    /// This error occurs when the start of the interval comes _after_ the end
    /// of the interval.
    ///
    /// On positive stranded intervals, this is when the start position is
    /// _greater than_ the end position. On negative stranded intervals, this is
    /// when the start position is _less than_ the end position.
    #[error("negatively sized interval: start is `{start}`, end is `{end}`, strand is `{strand}`")]
    NegativelySized {
        /// The start position.
        start: Number,
        /// The end position.
        end: Number,
        /// The strand.
        strand: Strand,
    },
}

/// A [`Result`](std::result::Result) with a [`NonsensicalError`].
pub type NonsensicalResult<T> = std::result::Result<T, NonsensicalError>;

/// An error related to parsing an interval.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// An invalid format was encountered.
    #[error("invalid format: {value}")]
    Format {
        /// The value that was passed.
        value: String,
    },
}

/// A [`Result`](std::result::Result) with a [`ParseError`].
pub type ParseResult<T> = std::result::Result<T, ParseError>;

/// An error related to an interval.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// A clamping error.
    #[error("clamp error: {0}")]
    Clamp(#[from] ClampError),

    /// A coordinate error.
    #[error("coordinate error: {0}")]
    Coordinate(#[from] coordinate::Error),

    /// A nonsensical interval.
    #[error("nonsensical interval: {0}")]
    Nonsensical(#[from] NonsensicalError),

    /// One or more of the coordinates were out of bounds.
    #[error("one or more of the coordinates were out of bounds")]
    OutOfBounds,

    /// A parse error.
    #[error("parse error: {0}")]
    Parse(#[from] ParseError),

    /// A position error.
    #[error("position error: {0}")]
    Position(#[from] position::Error),

    /// A strand error.
    #[error("strand error: {0}")]
    Strand(#[from] strand::Error),
}

/// A [`Result`](std::result::Result) with an [`Error`](enum@Error).
pub type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// The `Coordinate` trait
////////////////////////////////////////////////////////////////////////////////////////

/// Traits related to a coordinate.
pub mod r#trait {
    use super::*;
    use crate::system::Base;

    /// Requirements to be an interval.
    #[allow(clippy::len_without_is_empty)]
    pub trait Interval<S: System> {
        /// Returns whether or not the entity at the in-base coordinate is
        /// contained within this interval.
        fn contains_entity(&self, coordinate: &Coordinate<Base>) -> bool;

        /// Gets the number of member contained within the interval.
        fn count_entities(&self) -> Number;
    }
}

/// An interval.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Interval<S: System> {
    /// The start coordinate.
    start: Coordinate<S>,

    /// The end coordinate.
    end: Coordinate<S>,
}

impl<S: System> Interval<S>
where
    Interval<S>: r#trait::Interval<S>,
    Position<S>: position::r#trait::Position<S>,
{
    /// Creates a new interval if the following invariants are upheld.
    ///
    /// * The contigs of the two coordinates must match.
    ///   * If this does not hold, a [`NonsensicalError::MismatchedContigs`]
    ///     will be returned.
    /// * The strands of the two coordinates must match.
    ///   * If this does not hold, a [`NonsensicalError::MismatchedStrands`]
    ///     will be returned.
    /// * The start must come _before or be equal to_ the end in that (a) on
    ///   positive strand, `start <= end`, or, (b) on the negative strand, `end
    ///   <= start`. This ensures that the interval is always oriented from
    ///   start to end of the molecule.
    ///   * If this does not hold, a [`NonsensicalError::NegativelySized`] will
    ///     be returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// // Positive strand.
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// // Negative strand.
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// // Positive strand.
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// // Negative strand.
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(start: Coordinate<S>, end: Coordinate<S>) -> Result<super::Interval<S>> {
        if start.contig() != end.contig() {
            return Err(Error::Nonsensical(NonsensicalError::MismatchedContigs {
                start: start.contig().clone(),
                end: end.contig().clone(),
            }));
        }

        if start.strand() != end.strand() {
            return Err(Error::Nonsensical(NonsensicalError::MismatchedStrands {
                start: start.strand(),
                end: end.strand(),
            }));
        }

        match start.strand() {
            Strand::Positive => {
                if start.position() > end.position() {
                    return Err(Error::Nonsensical(NonsensicalError::NegativelySized {
                        start: start.position().get(),
                        end: end.position().get(),
                        strand: start.strand(),
                    }));
                }
            }
            Strand::Negative => {
                if end.position() > start.position() {
                    return Err(Error::Nonsensical(NonsensicalError::NegativelySized {
                        start: start.position().get(),
                        end: end.position().get(),
                        strand: start.strand(),
                    }));
                }
            }
        }

        Ok(Interval { start, end })
    }

    /// Gets a reference to the start coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start.clone(), end)?;
    ///
    /// assert_eq!(interval.start(), &start);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start.clone(), end)?;
    ///
    /// assert_eq!(interval.start(), &start);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn start(&self) -> &Coordinate<S> {
        &self.start
    }

    /// Consumes `self` and returns the start coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start.clone(), end)?;
    ///
    /// assert_eq!(interval.into_start(), start);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start.clone(), end)?;
    ///
    /// assert_eq!(interval.into_start(), start);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_start(self) -> Coordinate<S> {
        self.start
    }

    /// Gets a reference to the end coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end.clone())?;
    ///
    /// assert_eq!(interval.end(), &end);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end.clone())?;
    ///
    /// assert_eq!(interval.end(), &end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> &Coordinate<S> {
        &self.end
    }

    /// Consumes `self` and returns the end coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end.clone())?;
    ///
    /// assert_eq!(interval.into_end(), end);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end.clone())?;
    ///
    /// assert_eq!(interval.into_end(), end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_end(self) -> Coordinate<S> {
        self.end
    }

    /// Consumes `self` and returns the start and end coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start.clone(), end.clone())?;
    /// let parts = interval.into_coordinates();
    ///
    /// assert_eq!(parts.0, start);
    /// assert_eq!(parts.1, end);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start.clone(), end.clone())?;
    /// let parts = interval.into_coordinates();
    ///
    /// assert_eq!(parts.0, start);
    /// assert_eq!(parts.1, end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_coordinates(self) -> (Coordinate<S>, Coordinate<S>) {
        (self.start, self.end)
    }

    /// Returns a reference to the contig.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.contig().as_str(), "seq0");
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.contig().as_str(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        self.start().contig()
    }

    /// Returns the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.strand(), Strand::Positive);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.strand(), Strand::Negative);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> Strand {
        self.start().strand()
    }

    /// Returns whether or not a coordinate is contained within this interval.
    /// Notably, when checked whether coordinates are included in the interval,
    /// both the start and end positions are considered inclusive.
    ///
    /// # Caution
    ///
    /// **This is not the method you want to use when checking if a nucleotide
    /// or amino acid at a particular position is included in the interval. This
    /// checks the coordinates themselves and, in-so-doing, considers both the
    /// start and the end positions of the interval to be inclusive.
    ///
    /// This method checks containment using coordinates in the interval's
    /// native coordinate system (the generic type `S`). If you'd like to
    /// check whether a particular nucleotide, amino acid, or other entity
    /// is contained within the interval (using in-base coordinates), use
    /// the [`contains_entity()`](Interval::contains_entity) method.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 0)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// // Coordinates on the same contig, strand, and within the interval's range
    /// // are contained within the interval.
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 0)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 5)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 10)?));
    ///
    /// // Coordinates on different contigs, strands, or outside the range are
    /// // not contained within the interval.
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq1", "+", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "-", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 11)?));
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 1)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// // Coordinates on the same contig, strand, and within the interval's range
    /// // are contained within the interval.
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 1)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 5)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 10)?));
    ///
    /// // Coordinates on different contigs, strands, or outside the range are
    /// // not contained within the interval.
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq1", "+", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "-", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 11)?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contains_coordinate(&self, coordinate: &crate::Coordinate<S>) -> bool {
        if self.contig() != coordinate.contig() {
            return false;
        }

        if self.strand() != coordinate.strand() {
            return false;
        }

        match self.strand() {
            Strand::Positive => {
                self.start().position().get() <= coordinate.position().get()
                    && self.end().position().get() >= coordinate.position().get()
            }
            Strand::Negative => {
                self.start().position().get() >= coordinate.position().get()
                    && self.end().position().get() <= coordinate.position().get()
            }
        }
    }

    /// Returns whether or not the entity at the in-base coordinate is
    /// contained within this interval.
    ///
    /// This method always works with in-base coordinates (which directly point
    /// to entities like nucleotides or amino acids), regardless of the
    /// interval's coordinate system. Use
    /// [`contains_coordinate()`](Self::contains_coordinate) if you need to
    /// check containment using the interval's native coordinate system.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 0)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// // Coordinates on the same contig, strand, and within the interval's range
    /// // are contained within the interval.
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 0)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 5)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 10)?));
    ///
    /// // Coordinates on different contigs, strands, or outside the range are
    /// // not contained within the interval.
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq1", "+", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "-", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 11)?));
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 1)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// // Coordinates on the same contig, strand, and within the interval's range
    /// // are contained within the interval.
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 1)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 5)?));
    /// assert!(interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 10)?));
    ///
    /// // Coordinates on different contigs, strands, or outside the range are
    /// // not contained within the interval.
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq1", "+", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "-", 5)?));
    /// assert!(!interval.contains_coordinate(&Coordinate::try_new("seq0", "+", 11)?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contains_entity(&self, coordinate: &Coordinate<Base>) -> bool {
        <Self as r#trait::Interval<S>>::contains_entity(self, coordinate)
    }

    /// Counts the number of entities in the interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// // Positive strand.
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.count_entities(), 10);
    ///
    /// // Negative strand.
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.count_entities(), 10);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// // Positive strand.
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.count_entities(), 11);
    ///
    /// // Negative strand.
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// assert_eq!(interval.count_entities(), 11);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn count_entities(&self) -> Number {
        <Self as r#trait::Interval<S>>::count_entities(self)
    }

    /// Consumes `self` and clamps an interval by another interval.
    ///
    /// Clamping is an operation whereby the ends of an interval are restricted
    /// to the range of the argument passed in with a tendency to restrict
    /// towards the middle of the interval.
    ///
    /// # Summary
    ///
    /// * If the interval being operated on is completely contained within the
    ///   argument interval, the interval being operated on is returned.
    ///
    /// ```text
    /// ╔═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════ →
    /// 10    11    12    13    14    15    16    17    18    19    20       |
    ///                   ●───────────────────────● [13, 17]                 | Original Interval
    ///       ●───────────────────────────────────────────────● [11, 19]     | Argument Interval
    /// ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄
    ///                   ●───────────────────────● [13, 17]                 | Resulting Interval
    ///
    ///
    /// Here, no modifications were made to the original interval, as neither
    /// the start nor the end of the interval would be restricted by the
    /// argument interval.
    /// ```
    ///
    /// * If the argument interval is completely within the interval being
    ///   operated on, the argument interval will clamp both sides of the
    ///   original interval, and the argument interval will be returned.
    /// ```text
    /// ╔═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════ →
    /// 10    11    12    13    14    15    16    17    18    19    20       |
    ///       ●───────────────────────────────────────────────● [11, 19]     | Original Interval
    ///                   ●───────────────────────● [13, 17]                 | Argument Interval
    /// ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄
    ///                   ●───────────────────────● [13, 17]                 | Resulting Interval
    ///
    ///
    /// Here, both the start and the end position of the original interval were
    /// restricted by the start and end of the argument interval respectively.
    /// ```
    ///
    /// * If the argument interval would restrict the length of one side of the
    ///   subject interval on either end, that end is restricted to the argument
    ///   interval's value, whereas the non-restricted end is the original
    ///   interval's value.
    ///
    /// ```text
    /// ╔═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════ →
    /// 10    11    12    13    14    15    16    17    18    19    20       |
    ///       ●───────────────────────────────────● [11, 17]                 | Original Interval
    ///                   ●───────────────────────● [13, 17]                 | Argument Interval
    /// ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄
    ///                   ●───────────────────────● [13, 17]                 | Resulting Interval
    ///
    ///
    /// Here, the start of the original interval is clamped by the argument
    /// interval's start position. However, the end position of the original
    /// interval is not restricted by the argument interval's end position,
    /// so it remains the same. This results in the latter half of the interval
    /// being clamped.
    ///
    ///
    /// ╔═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════╦═════ →
    /// 10    11    12    13    14    15    16    17    18    19    20       |
    ///                   ●───────────────────────────────────● [13, 19]     | Original Interval
    ///                   ●───────────────────────● [13, 17]                 | Argument Interval
    /// ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄
    ///                   ●───────────────────────● [13, 17]                 | Resulting Interval
    ///
    ///
    /// Here, the start position of the original interval would not be
    /// restricted by the argument interval's start position, so it remains
    /// the same. However, the end position is clamped by the end position
    /// of the argument interval, so the resulting end position is that of the
    /// argument interval's end position. This results in the first half of
    /// interval being clamped.
    /// ```
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let interval = "seq0:+:10-20".parse::<Interval<Interbase>>()?;
    /// let clamped = interval.clamp("seq0:+:5-15".parse::<Interval<Interbase>>()?)?;
    /// assert_eq!(clamped, "seq0:+:10-15".parse::<Interval<Interbase>>()?);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let interval = "seq0:-:20-10".parse::<Interval<Base>>()?;
    /// let clamped = interval.clamp("seq0:-:25-15".parse::<Interval<Base>>()?)?;
    /// assert_eq!(clamped, "seq0:-:20-15".parse::<Interval<Base>>()?);
    ///
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new interval"]
    pub fn clamp(self, interval: Interval<S>) -> Result<Interval<S>> {
        let (start, end) = self.into_coordinates();
        let (operand_start, operand_end) = interval.into_coordinates();

        let (start_contig, start_strand, start) = start.into_parts();
        let (end_contig, end_strand, end) = end.into_parts();

        let (operand_contig, operand_strand, operand_start) = operand_start.into_parts();
        let (_, _, operand_end) = operand_end.into_parts();

        if start_contig != operand_contig {
            return Err(Error::Clamp(ClampError::MismatchedContigs {
                original: start_contig,
                operand: operand_contig,
            }));
        }

        if start_strand != operand_strand {
            return Err(Error::Clamp(ClampError::MismatchedStrand {
                original: start_strand,
                operand: operand_strand,
            }));
        }

        let (new_start, new_end) = match start_strand {
            Strand::Positive => (max(start, operand_start), min(end, operand_end)),
            Strand::Negative => (min(start, operand_start), max(end, operand_end)),
        };

        let start = Coordinate::<S>::new(start_contig, start_strand, new_start);
        let end = Coordinate::<S>::new(end_contig, end_strand, new_end);

        // SAFETY: both the start _and_ the end positions were originally on
        // intervals that were valid. Since we are not breaking any rules that
        // would make the intervals invalid in this method, this should always
        // unwrap.
        Ok(Self::try_new(start, end).unwrap())
    }

    /// Gets the offset of a coordinate from the start of the interval.
    ///
    /// If the coordinate is not contained within the interval, `None` is
    /// returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// let query = Coordinate::<Interbase>::try_new("seq0", "+", 15)?;
    /// assert_eq!(interval.coordinate_offset(&query).unwrap(), 5);
    ///
    /// let query = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// assert_eq!(interval.coordinate_offset(&query).unwrap(), 10);
    ///
    /// let query = Coordinate::<Interbase>::try_new("seq0", "+", 21)?;
    /// assert!(interval.coordinate_offset(&query).is_none());
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_new(start, end)?;
    ///
    /// let query = Coordinate::<Base>::try_new("seq0", "-", 15)?;
    /// assert_eq!(interval.coordinate_offset(&query).unwrap(), 5);
    ///
    /// let query = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// assert_eq!(interval.coordinate_offset(&query).unwrap(), 10);
    ///
    /// let query = Coordinate::<Base>::try_new("seq0", "-", 9)?;
    /// assert!(interval.coordinate_offset(&query).is_none());
    ///
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn coordinate_offset(&self, coordinate: &Coordinate<S>) -> Option<Number> {
        if !self.contains_coordinate(coordinate) {
            return None;
        }

        Some(
            coordinate
                .position()
                .distance_unchecked(self.start().position()),
        )
    }

    /// Returns the coordinate at the offset within the interval.
    ///
    /// This method only returns the coordinate if the coordinate falls within
    /// the interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// // Positive strand.
    ///
    /// let interval = "seq0:+:0-1000".parse::<Interval<Interbase>>()?;
    ///
    /// let expected = "seq0:+:5".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(interval.coordinate_at_offset(5).unwrap(), expected);
    ///
    /// let expected = "seq0:+:1000".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(interval.coordinate_at_offset(1000).unwrap(), expected);
    ///
    /// assert!(interval.coordinate_at_offset(1001).is_none());
    ///
    /// // Negative strand.
    ///
    /// let interval = "seq0:-:1000-0".parse::<Interval<Interbase>>()?;
    ///
    /// let expected = "seq0:-:995".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(interval.coordinate_at_offset(5).unwrap(), expected);
    ///
    /// let expected = "seq0:-:0".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(interval.coordinate_at_offset(1000).unwrap(), expected);
    ///
    /// assert_eq!(interval.coordinate_at_offset(1001), None);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// // Positive strand.
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<Base>>()?;
    ///
    /// let expected = "seq0:+:6".parse::<Coordinate<Base>>()?;
    /// assert_eq!(interval.coordinate_at_offset(5).unwrap(), expected);
    ///
    /// let expected = "seq0:+:1000".parse::<Coordinate<Base>>()?;
    /// assert_eq!(interval.coordinate_at_offset(999).unwrap(), expected);
    ///
    /// assert!(interval.coordinate_at_offset(1000).is_none());
    ///
    /// // Negative strand.
    ///
    /// let interval = "seq0:-:1000-1".parse::<Interval<Base>>()?;
    ///
    /// let expected = "seq0:-:995".parse::<Coordinate<Base>>()?;
    /// assert_eq!(interval.coordinate_at_offset(5).unwrap(), expected);
    ///
    /// let expected = "seq0:-:1".parse::<Coordinate<Base>>()?;
    /// assert_eq!(interval.coordinate_at_offset(999).unwrap(), expected);
    ///
    /// assert_eq!(interval.coordinate_at_offset(1000), None);
    ///
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn coordinate_at_offset(&self, offset: Number) -> Option<Coordinate<S>> {
        let coordinate = self.start().clone().move_forward(offset)?;

        match self.contains_coordinate(&coordinate) {
            true => Some(coordinate),
            false => None,
        }
    }

    /// Reverse complements the interval, meaning that:
    ///
    /// * the start and end positions are swapped, and
    /// * the strand is swapped.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    /// let original = Interval::try_new(start, end)?;
    ///
    /// let complemented = original.clone().reverse_complement();
    /// assert_eq!(complemented, "seq0:-:20-10".parse::<Interval<Interbase>>()?);
    ///
    /// let recomplemented = complemented.reverse_complement();
    /// assert_eq!(recomplemented, original);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let original = Interval::try_new(start, end)?;
    ///
    /// let complemented = original.clone().reverse_complement();
    /// assert_eq!(complemented, "seq0:-:20-10".parse::<Interval<Base>>()?);
    ///
    /// let recomplemented = complemented.reverse_complement();
    /// assert_eq!(recomplemented, original);
    ///
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new interval"]
    pub fn reverse_complement(self) -> super::Interval<S> {
        let (start, end) = self.into_coordinates();
        // SAFETY: because (a) intervals are inclusive of both of their start
        // and end coordinates, (b) all positions can be represented on the
        // opposite strand, and (c) swapping the start and end while also
        // swapping strand will always create the correct directionality, this will
        // always unwrap.
        Interval::try_new(end.swap_strand(), start.swap_strand()).unwrap()
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl<S: System> std::fmt::Display for Interval<S>
where
    Interval<S>: r#trait::Interval<S>,
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

impl<S: System> std::str::FromStr for Interval<S>
where
    Interval<S>: r#trait::Interval<S>,
    Position<S>: position::r#trait::Position<S>,
{
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts = s.split(':').collect::<Vec<_>>();

        if parts.len() != 3 {
            return Err(Error::Parse(ParseError::Format {
                value: s.to_string(),
            }));
        }

        let mut parts = parts.iter();

        // SAFETY: we checked that there are three parts above. Given that we
        // haven't pulled anything from the iterator, we can always safely
        // unwrap this.
        let contig = parts.next().unwrap().parse::<Contig>().map_err(|_| {
            Error::Parse(ParseError::Format {
                value: s.to_string(),
            })
        })?;

        // SAFETY: we checked that there are three parts above. Given that we
        // have only pulled one item from the iterator, we can always safely
        // unwrap this.
        let strand = parts
            .next()
            .unwrap()
            .parse::<Strand>()
            .map_err(Error::Strand)?;

        // SAFETY: we checked that there are three parts above. Given that we
        // have only pulled two items from the iterator, we can always safely
        // unwrap this.
        let positions = parts.next().unwrap().split('-').collect::<Vec<_>>();

        if positions.len() != 2 {
            return Err(Error::Parse(ParseError::Format {
                value: s.to_string(),
            }));
        }

        // SAFETY: we just ensured that two parts exist, so the direct
        // indexing of the slice for both index zero and one will never
        // fail.
        let start = positions[0]
            .parse::<Position<S>>()
            .map_err(Error::Position)?;
        let end = positions[1]
            .parse::<Position<S>>()
            .map_err(Error::Position)?;

        Interval::try_new(
            Coordinate::new(contig.clone(), strand, start),
            Coordinate::new(contig, strand, end),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::position::Error as PositionError;
    use crate::position::Number;
    use crate::position::ParseError as PositionParseError;
    use crate::strand::Error as StrandError;
    use crate::strand::ParseError as StrandParseError;
    use crate::system::Interbase;

    #[test]
    fn valid() {
        let start = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:+:9".parse::<Coordinate<Interbase>>().unwrap();

        let interval = Interval::try_new(start, end).unwrap();
        assert_eq!(interval.count_entities(), 9);
    }

    #[test]
    fn nonsensical_mismatched_contigs() {
        let start = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq1:+:10".parse::<Coordinate<Interbase>>().unwrap();

        let err = Interval::try_new(start, end).unwrap_err();
        assert_eq!(
            err,
            Error::Nonsensical(NonsensicalError::MismatchedContigs {
                start: Contig::new_unchecked("seq0"),
                end: Contig::new_unchecked("seq1")
            })
        );

        assert_eq!(
            err.to_string(),
            "nonsensical interval: mismatched contigs for coordinates: `seq0` and `seq1`"
        );
    }

    #[test]
    fn nonsensical_mismatched_strands() {
        let start = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:-:10".parse::<Coordinate<Interbase>>().unwrap();

        let err = Interval::try_new(start, end).unwrap_err();
        assert_eq!(
            err,
            Error::Nonsensical(NonsensicalError::MismatchedStrands {
                start: Strand::Positive,
                end: Strand::Negative
            })
        );

        assert_eq!(
            err.to_string(),
            "nonsensical interval: mismatched strands for coordinates: `+` and `-`"
        );
    }

    #[test]
    fn nonsensical_start_greater_than_end() {
        //===================//
        // Positive stranded //
        //===================//

        let start = "seq0:+:10".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();

        let err = Interval::try_new(start, end).unwrap_err();

        assert_eq!(
            err,
            Error::Nonsensical(NonsensicalError::NegativelySized {
                start: 10,
                end: 0,
                strand: Strand::Positive
            })
        );

        assert_eq!(
            err.to_string(),
            "nonsensical interval: negatively sized interval: start is `10`, end is `0`, strand \
             is `+`"
        );

        //===================//
        // Negative stranded //
        //===================//

        let start = "seq0:-:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:-:10".parse::<Coordinate<Interbase>>().unwrap();

        let err = Interval::try_new(start, end).unwrap_err();

        assert_eq!(
            err,
            Error::Nonsensical(NonsensicalError::NegativelySized {
                start: 0,
                end: 10,
                strand: Strand::Negative
            })
        );

        assert_eq!(
            err.to_string(),
            "nonsensical interval: negatively sized interval: start is `0`, end is `10`, strand \
             is `-`"
        );
    }

    #[test]
    fn zero_sized() {
        let start = "seq0:+:10".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:+:10".parse::<Coordinate<Interbase>>().unwrap();

        let interval = Interval::try_new(start.clone(), end.clone()).unwrap();
        assert!(interval.end().position().get() - interval.start().position().get() == 0);
        assert!(interval.contains_coordinate(&start));
        assert!(interval.contains_coordinate(&end));
        assert!(
            !interval.contains_coordinate(&"seq0:+:9".parse::<Coordinate<Interbase>>().unwrap())
        );
        assert!(
            !interval.contains_coordinate(&"seq0:+:11".parse::<Coordinate<Interbase>>().unwrap())
        );
    }

    #[test]
    fn positive_strand_clamp() {
        let interval = "seq0:+:1000-2000".parse::<Interval<Interbase>>().unwrap();

        assert_eq!(
            interval
                .clone()
                .clamp("seq1:+:0-3000".parse::<Interval<Interbase>>().unwrap()),
            Err(Error::Clamp(ClampError::MismatchedContigs {
                original: Contig::new_unchecked("seq0"),
                operand: Contig::new_unchecked("seq1")
            }))
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:-:3000-0".parse::<Interval<Interbase>>().unwrap()),
            Err(Error::Clamp(ClampError::MismatchedStrand {
                original: Strand::Positive,
                operand: Strand::Negative
            }))
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:+:0-3000".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:+:1000-2000".parse::<Interval<Interbase>>().unwrap()
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:+:1250-3000".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:+:1250-2000".parse::<Interval<Interbase>>().unwrap()
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:+:0-1750".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:+:1000-1750".parse::<Interval<Interbase>>().unwrap()
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:+:1250-1750".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:+:1250-1750".parse::<Interval<Interbase>>().unwrap()
        );
    }

    #[test]
    fn negative_strand_clamp() {
        let interval = "seq0:-:2000-1000".parse::<Interval<Interbase>>().unwrap();

        assert_eq!(
            interval
                .clone()
                .clamp("seq1:-:3000-0".parse::<Interval<Interbase>>().unwrap()),
            Err(Error::Clamp(ClampError::MismatchedContigs {
                original: Contig::new_unchecked("seq0"),
                operand: Contig::new_unchecked("seq1")
            }))
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:+:0-3000".parse::<Interval<Interbase>>().unwrap()),
            Err(Error::Clamp(ClampError::MismatchedStrand {
                original: Strand::Negative,
                operand: Strand::Positive
            }))
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:-:3000-0".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:-:2000-1000".parse::<Interval<Interbase>>().unwrap()
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:-:3000-1250".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:-:2000-1250".parse::<Interval<Interbase>>().unwrap()
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:-:1750-0".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:-:1750-1000".parse::<Interval<Interbase>>().unwrap()
        );

        assert_eq!(
            interval
                .clone()
                .clamp("seq0:-:1750-1250".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:-:1750-1250".parse::<Interval<Interbase>>().unwrap()
        );
    }

    #[test]
    fn positive_strand_offset() {
        let interval = "seq0:+:1000-2000".parse::<Interval<Interbase>>().unwrap();

        // Mismatched contigs means the interval does not contain the coordinate.
        let coordinate = "seq1:+:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        // Mismatched strands means the interval does not contain the coordinate.
        let coordinate = "seq0:-:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        // Contained within.
        let coordinate = "seq0:+:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert_eq!(interval.coordinate_offset(&coordinate).unwrap(), 0);

        let coordinate = "seq0:+:2000".parse::<Coordinate<Interbase>>().unwrap();
        assert_eq!(interval.coordinate_offset(&coordinate).unwrap(), 1000);

        // Just outside of range.
        let coordinate = "seq0:+:999".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        let coordinate = "seq0:+:2001".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());
    }

    #[test]
    fn negative_strand_offset() {
        let interval = "seq0:-:2000-1000".parse::<Interval<Interbase>>().unwrap();

        // Mismatched contigs means the interval does not contain the coordinate.
        let coordinate = "seq1:-:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        // Mismatched strands means the interval does not contain the coordinate.
        let coordinate = "seq0:+:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        // Contained within.
        let coordinate = "seq0:-:2000".parse::<Coordinate<Interbase>>().unwrap();
        assert_eq!(interval.coordinate_offset(&coordinate).unwrap(), 0);

        let coordinate = "seq0:-:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert_eq!(interval.coordinate_offset(&coordinate).unwrap(), 1000);

        // Just outside of range.
        let coordinate = "seq0:-:999".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        let coordinate = "seq0:-:2001".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());
    }

    #[test]
    fn len() {
        assert_eq!(
            "seq0:+:0-1000"
                .parse::<Interval<Interbase>>()
                .unwrap()
                .count_entities(),
            1000
        );

        assert_eq!(
            "seq0:-:1000-0"
                .parse::<Interval<Interbase>>()
                .unwrap()
                .count_entities(),
            1000
        );
        let interval = "seq0:-:2000-1000".parse::<Interval<Interbase>>().unwrap();

        // Mismatched contigs means the interval does not contain the coordinate.
        let coordinate = "seq1:-:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        // Mismatched strands means the interval does not contain the coordinate.
        let coordinate = "seq0:+:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        // Contained within.
        let coordinate = "seq0:-:2000".parse::<Coordinate<Interbase>>().unwrap();
        assert_eq!(interval.coordinate_offset(&coordinate).unwrap(), 0);

        let coordinate = "seq0:-:1000".parse::<Coordinate<Interbase>>().unwrap();
        assert_eq!(interval.coordinate_offset(&coordinate).unwrap(), 1000);

        // Just outside of range.
        let coordinate = "seq0:-:999".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());

        let coordinate = "seq0:-:2001".parse::<Coordinate<Interbase>>().unwrap();
        assert!(interval.coordinate_offset(&coordinate).is_none());
    }

    #[test]
    fn parse() {
        let value = format!("seq0:+:0-{}", Number::MAX);
        let interval = value.parse::<Interval<Interbase>>().unwrap();
        assert_eq!(interval.contig().as_str(), "seq0");
        assert_eq!(interval.strand(), Strand::Positive);
        assert_eq!(interval.start().position().get(), 0);
        assert_eq!(interval.end().position().get(), Number::MAX);

        let value = format!("seq0:-:{}-0", Number::MAX);
        let interval = value.parse::<Interval<Interbase>>().unwrap();
        assert_eq!(interval.contig().as_str(), "seq0");
        assert_eq!(interval.strand(), Strand::Negative);
        assert_eq!(interval.start().position().get(), Number::MAX);
        assert_eq!(interval.end().position().get(), 0);
    }

    #[test]
    fn parse_error() {
        let err = "1".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Parse(ParseError::Format {
                value: String::from("1")
            })
        );

        let err = "1-1000".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Parse(ParseError::Format {
                value: String::from("1-1000")
            })
        );

        let err = "seq0:".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Parse(ParseError::Format {
                value: String::from("seq0:")
            })
        );

        let err = "seq0:0-".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Parse(ParseError::Format {
                value: String::from("seq0:0-")
            })
        );

        let err = "seq0:0-10000:".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Strand(StrandError::Parse(StrandParseError::Invalid {
                value: String::from("0-10000")
            }))
        );

        let err = "seq0:+".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Parse(ParseError::Format {
                value: String::from("seq0:+")
            })
        );

        let err = "seq0:+:0".parse::<Interval<Interbase>>().unwrap_err();
        assert_eq!(
            err,
            Error::Parse(ParseError::Format {
                value: String::from("seq0:+:0")
            })
        );

        let err = "seq0:+:0-".parse::<Interval<Interbase>>().unwrap_err();
        assert!(matches!(
            err,
            Error::Position(PositionError::Parse(PositionParseError::Int { .. }))
        ));

        let err = "seq0:+:0-$".parse::<Interval<Interbase>>().unwrap_err();
        assert!(matches!(
            err,
            Error::Position(PositionError::Parse(PositionParseError::Int { .. }))
        ));
    }

    #[test]
    fn to_string() {
        // Positive-stranded interval
        let start = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:+:10".parse::<Coordinate<Interbase>>().unwrap();
        let interval = Interval::try_new(start, end).unwrap();

        assert_eq!(interval.to_string(), "seq0:+:0-10");

        // Negative-stranded interval
        let start = "seq0:-:10".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:-:0".parse::<Coordinate<Interbase>>().unwrap();
        let interval = Interval::try_new(start, end).unwrap();

        assert_eq!(interval.to_string(), "seq0:-:10-0");
    }
}
