//! Intervals.

use omics_core::VARIANT_SEPARATOR;
use thiserror::Error;

use crate::Contig;
use crate::Position;
use crate::Span;
use crate::Strand;
use crate::System;
use crate::contig;
use crate::coordinate;
use crate::coordinate::Coordinate;
use crate::coordinate::CoordinateRef;
use crate::position;
use crate::position::Number;
use crate::span;
use crate::span::Direction;
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

    /// A disjoint interval error.
    ///
    /// This error occurs when one attempts to clamp an interval with another
    /// interval that does not overlap it.
    #[error(
        "disjoint intervals: `{original_start}-{original_end}` and \
         `{operand_start}-{operand_end}` on strand `{strand}`"
    )]
    Disjoint {
        /// The start position of the interval being clamped.
        original_start: Number,

        /// The end position of the interval being clamped.
        original_end: Number,

        /// The start position of the interval doing the clamping.
        operand_start: Number,

        /// The end position of the interval doing the clamping.
        operand_end: Number,

        /// The strand shared by both intervals.
        strand: Strand,
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

    /// A contig error.
    #[error("contig error: {0}")]
    Contig(#[from] contig::Error),

    /// A coordinate error.
    #[error("coordinate error: {0}")]
    Coordinate(#[from] coordinate::Error),

    /// The span direction does not agree with the strand.
    #[error("span direction `{direction}` is incompatible with strand `{strand}`")]
    DirectionMismatch {
        /// The interval strand.
        strand: Strand,

        /// The span direction.
        direction: Direction,
    },

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
    /// The contig.
    contig: Contig,

    /// The strand.
    strand: Strand,

    /// The directed span.
    span: Span<S>,
}

impl<S: System> Interval<S>
where
    Interval<S>: r#trait::Interval<S>,
    Position<S>: position::r#trait::Position<S>,
{
    /// Creates an interval from genomic context and a directed span.
    ///
    /// Positive strands accept ascending and stationary spans. Negative
    /// strands accept descending and stationary spans.
    ///
    /// # Errors
    ///
    /// Returns [`Error::DirectionMismatch`] when the span direction does not
    /// agree with the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Span;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let span = Span::<Interbase>::try_new(10, 20)?;
    /// let interval = Interval::try_new("seq0", "+", span)?;
    /// assert_eq!(interval.to_string(), "seq0:+:10-20");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        contig: impl TryInto<Contig, Error = contig::Error>,
        strand: impl TryInto<Strand, Error = strand::Error>,
        span: Span<S>,
    ) -> Result<Self> {
        let contig = contig.try_into()?;
        let strand = strand.try_into()?;
        let compatible = matches!(
            (strand, span.direction()),
            (
                Strand::Positive,
                Direction::Ascending | Direction::Stationary
            ) | (
                Strand::Negative,
                Direction::Descending | Direction::Stationary
            )
        );

        if !compatible {
            return Err(Error::DirectionMismatch {
                strand,
                direction: span.direction(),
            });
        }

        Ok(Self {
            contig,
            strand,
            span,
        })
    }

    /// Gets a borrowed start coordinate.
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
    /// let interval = Interval::try_from((start.clone(), end))?;
    ///
    /// assert_eq!(interval.start().into_owned(), start);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_from((start.clone(), end))?;
    ///
    /// assert_eq!(interval.start().into_owned(), start);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn start(&self) -> CoordinateRef<'_, S> {
        CoordinateRef::new(&self.contig, self.strand, self.span.start())
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
    /// let interval = Interval::try_from((start.clone(), end))?;
    ///
    /// assert_eq!(interval.into_start(), start);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_from((start.clone(), end))?;
    ///
    /// assert_eq!(interval.into_start(), start);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_start(self) -> Coordinate<S> {
        let (start, _) = self.span.into_positions();
        Coordinate::new(self.contig, self.strand, start)
    }

    /// Gets a borrowed end coordinate.
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
    /// let interval = Interval::try_from((start, end.clone()))?;
    ///
    /// assert_eq!(interval.end().into_owned(), end);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_from((start, end.clone()))?;
    ///
    /// assert_eq!(interval.end().into_owned(), end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> CoordinateRef<'_, S> {
        CoordinateRef::new(&self.contig, self.strand, self.span.end())
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
    /// let interval = Interval::try_from((start, end.clone()))?;
    ///
    /// assert_eq!(interval.into_end(), end);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_from((start, end.clone()))?;
    ///
    /// assert_eq!(interval.into_end(), end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_end(self) -> Coordinate<S> {
        let (_, end) = self.span.into_positions();
        Coordinate::new(self.contig, self.strand, end)
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
    /// let interval = Interval::try_from((start.clone(), end.clone()))?;
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
    /// let interval = Interval::try_from((start.clone(), end.clone()))?;
    /// let parts = interval.into_coordinates();
    ///
    /// assert_eq!(parts.0, start);
    /// assert_eq!(parts.1, end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_coordinates(self) -> (Coordinate<S>, Coordinate<S>) {
        let (start, end) = self.span.into_positions();
        (
            Coordinate::new(self.contig.clone(), self.strand, start),
            Coordinate::new(self.contig, self.strand, end),
        )
    }

    /// Returns the directed span.
    pub fn span(&self) -> &Span<S> {
        &self.span
    }

    /// Consumes the interval and returns its directed span.
    pub fn into_span(self) -> Span<S> {
        self.span
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
    /// let interval = Interval::try_from((start, end))?;
    ///
    /// assert_eq!(interval.contig().as_str(), "seq0");
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "+", 10)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "+", 20)?;
    /// let interval = Interval::try_from((start, end))?;
    ///
    /// assert_eq!(interval.contig().as_str(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        &self.contig
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
    /// let interval = Interval::try_from((start, end))?;
    ///
    /// assert_eq!(interval.strand(), Strand::Positive);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_from((start, end))?;
    ///
    /// assert_eq!(interval.strand(), Strand::Negative);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> Strand {
        self.strand
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
    /// let interval = Interval::try_from((start, end))?;
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
    /// let interval = Interval::try_from((start, end))?;
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
        self.contig() == coordinate.contig()
            && self.strand() == coordinate.strand()
            && self.span.contains_position(coordinate.position())
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
    /// let interval = Interval::try_from((start, end))?;
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
    /// let interval = Interval::try_from((start, end))?;
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
    /// let interval = Interval::try_from((start, end))?;
    ///
    /// assert_eq!(interval.count_entities(), 10);
    ///
    /// // Negative strand.
    ///
    /// let start = Coordinate::<Interbase>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Interbase>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_from((start, end))?;
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
    /// let interval = Interval::try_from((start, end))?;
    ///
    /// assert_eq!(interval.count_entities(), 11);
    ///
    /// // Negative strand.
    ///
    /// let start = Coordinate::<Base>::try_new("seq0", "-", 20)?;
    /// let end = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let interval = Interval::try_from((start, end))?;
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
        let Self {
            contig,
            strand,
            span,
        } = self;
        let Self {
            contig: operand_contig,
            strand: operand_strand,
            span: operand_span,
        } = interval;

        if contig != operand_contig {
            return Err(Error::Clamp(ClampError::MismatchedContigs {
                original: contig,
                operand: operand_contig,
            }));
        }

        if strand != operand_strand {
            return Err(Error::Clamp(ClampError::MismatchedStrand {
                original: strand,
                operand: operand_strand,
            }));
        }

        let span = span.clamp(operand_span).map_err(|error| match error {
            span::ClampError::Disjoint {
                original_start,
                original_end,
                operand_start,
                operand_end,
            } => Error::Clamp(ClampError::Disjoint {
                original_start,
                original_end,
                operand_start,
                operand_end,
                strand,
            }),
            span::ClampError::DirectionMismatch { .. } => {
                // SAFETY: intervals on one strand always have compatible span directions.
                unreachable!()
            }
        })?;

        Ok(Self {
            contig,
            strand,
            span,
        })
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
    /// let interval = Interval::try_from((start, end))?;
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
    /// let interval = Interval::try_from((start, end))?;
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
        (self.contig() == coordinate.contig() && self.strand() == coordinate.strand())
            .then(|| self.span.position_offset(coordinate.position()))
            .flatten()
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
        self.span
            .position_at_offset(offset)
            .map(|position| Coordinate::new(self.contig.clone(), self.strand, position))
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
    /// let original = Interval::try_from((start, end))?;
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
    /// let original = Interval::try_from((start, end))?;
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
        Self {
            contig: self.contig,
            strand: self.strand.complement(),
            span: self.span.reversed(),
        }
    }
}

impl<S: System> TryFrom<(Coordinate<S>, Coordinate<S>)> for Interval<S>
where
    Interval<S>: r#trait::Interval<S>,
    Position<S>: position::r#trait::Position<S>,
{
    type Error = Error;

    fn try_from((start, end): (Coordinate<S>, Coordinate<S>)) -> Result<Self> {
        let (start_contig, start_strand, start_position) = start.into_parts();
        let (end_contig, end_strand, end_position) = end.into_parts();

        if start_contig != end_contig {
            return Err(Error::Nonsensical(NonsensicalError::MismatchedContigs {
                start: start_contig,
                end: end_contig,
            }));
        }

        if start_strand != end_strand {
            return Err(Error::Nonsensical(NonsensicalError::MismatchedStrands {
                start: start_strand,
                end: end_strand,
            }));
        }

        let strand = match start_strand {
            Strand::Positive => "+",
            Strand::Negative => "-",
        };
        let span = Span::from((start_position, end_position));

        Self::try_new(start_contig.as_str(), strand, span)
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
        write!(f, "{}:{}:{}", self.contig(), self.strand(), self.span())
    }
}

impl<S: System> std::str::FromStr for Interval<S>
where
    Interval<S>: r#trait::Interval<S>,
    Position<S>: position::r#trait::Position<S>,
{
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let (prefix, positions) =
            s.rsplit_once(VARIANT_SEPARATOR)
                .ok_or_else(|| ParseError::Format {
                    value: s.to_owned(),
                })?;
        let (contig, strand) =
            prefix
                .rsplit_once(VARIANT_SEPARATOR)
                .ok_or_else(|| ParseError::Format {
                    value: s.to_owned(),
                })?;

        let contig = contig.parse::<Contig>().map_err(|_| {
            Error::Parse(ParseError::Format {
                value: s.to_string(),
            })
        })?;
        let strand = strand.parse::<Strand>().map_err(Error::Strand)?;
        let strand = match strand {
            Strand::Positive => "+",
            Strand::Negative => "-",
        };

        let span = positions.parse::<Span<S>>().map_err(|error| match error {
            span::ParseError::Format(_) => Error::Parse(ParseError::Format {
                value: s.to_owned(),
            }),
            span::ParseError::Start(error) | span::ParseError::End(error) => Error::Position(error),
        })?;

        Interval::try_new(contig.as_str(), strand, span)
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
    use crate::system::Base;
    use crate::system::Interbase;

    #[test]
    fn valid() {
        let start = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:+:9".parse::<Coordinate<Interbase>>().unwrap();

        // SAFETY: the endpoints share context and follow the positive strand.
        let interval = Interval::try_from((start, end)).unwrap();
        assert_eq!(interval.count_entities(), 9);
    }

    #[test]
    fn nonsensical_mismatched_contigs() {
        let start = "seq0:+:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq1:+:10".parse::<Coordinate<Interbase>>().unwrap();

        let err = Interval::try_from((start, end)).unwrap_err();
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

        let err = Interval::try_from((start, end)).unwrap_err();
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

        let err = Interval::try_from((start, end)).unwrap_err();

        assert_eq!(
            err,
            Error::DirectionMismatch {
                strand: Strand::Positive,
                direction: Direction::Descending,
            }
        );

        assert_eq!(
            err.to_string(),
            "span direction `descending` is incompatible with strand `+`"
        );

        //===================//
        // Negative stranded //
        //===================//

        let start = "seq0:-:0".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:-:10".parse::<Coordinate<Interbase>>().unwrap();

        let err = Interval::try_from((start, end)).unwrap_err();

        assert_eq!(
            err,
            Error::DirectionMismatch {
                strand: Strand::Negative,
                direction: Direction::Ascending,
            }
        );

        assert_eq!(
            err.to_string(),
            "span direction `ascending` is incompatible with strand `-`"
        );
    }

    #[test]
    fn zero_sized() {
        let start = "seq0:+:10".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:+:10".parse::<Coordinate<Interbase>>().unwrap();

        // SAFETY: stationary spans are valid on the positive strand.
        let interval = Interval::try_from((start.clone(), end.clone())).unwrap();
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
    fn interval_storage_is_smaller_than_two_coordinates() {
        assert!(
            std::mem::size_of::<Interval<Interbase>>()
                < 2 * std::mem::size_of::<Coordinate<Interbase>>()
        );
    }

    #[test]
    fn endpoint_views_borrow_interval_data() -> Result<()> {
        let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
        let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
        let interval = Interval::try_from((start.clone(), end.clone()))?;

        let start_ref: CoordinateRef<'_, Interbase> = interval.start();
        let end_ref: CoordinateRef<'_, Interbase> = interval.end();
        let start_position = interval.start().position();

        assert!(std::ptr::eq(start_ref.contig(), interval.contig()));
        assert!(std::ptr::eq(end_ref.contig(), interval.contig()));
        assert_eq!(start_ref, interval.start());
        assert_eq!(end_ref, interval.end());
        assert_eq!(start_ref.into_owned(), start);
        assert_eq!(end_ref.into_owned(), end);
        assert_eq!(start_ref.cmp(&end_ref), start.cmp(&end));
        assert_eq!(start_position.get(), 10);
        assert_eq!(start_ref.to_string(), "seq0:+:10");
        assert_eq!(end_ref.to_string(), "seq0:+:20");

        Ok(())
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
    fn clamp_touching_intervals() {
        let positive = "seq0:+:10-20".parse::<Interval<Interbase>>().unwrap();
        assert_eq!(
            positive
                .clamp("seq0:+:20-30".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:+:20-20".parse::<Interval<Interbase>>().unwrap()
        );

        let negative = "seq0:-:30-20".parse::<Interval<Interbase>>().unwrap();
        assert_eq!(
            negative
                .clamp("seq0:-:20-10".parse::<Interval<Interbase>>().unwrap())
                .unwrap(),
            "seq0:-:20-20".parse::<Interval<Interbase>>().unwrap()
        );
    }

    #[test]
    fn clamp_disjoint_intervals_returns_error() {
        let positive = "seq0:+:10-20".parse::<Interval<Interbase>>().unwrap();
        assert_eq!(
            positive
                .clamp("seq0:+:21-30".parse::<Interval<Interbase>>().unwrap())
                .unwrap_err(),
            Error::Clamp(ClampError::Disjoint {
                original_start: 10,
                original_end: 20,
                operand_start: 21,
                operand_end: 30,
                strand: Strand::Positive,
            })
        );

        let negative = "seq0:-:30-21".parse::<Interval<Interbase>>().unwrap();
        assert_eq!(
            negative
                .clamp("seq0:-:20-10".parse::<Interval<Interbase>>().unwrap())
                .unwrap_err(),
            Error::Clamp(ClampError::Disjoint {
                original_start: 30,
                original_end: 21,
                operand_start: 20,
                operand_end: 10,
                strand: Strand::Negative,
            })
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
    fn count_entities() {
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

        assert_eq!(
            "seq0:+:10-10"
                .parse::<Interval<Interbase>>()
                .unwrap()
                .count_entities(),
            0
        );

        assert_eq!(
            "seq0:+:10-20"
                .parse::<Interval<Base>>()
                .unwrap()
                .count_entities(),
            11
        );

        assert_eq!(
            "seq0:-:20-10"
                .parse::<Interval<Base>>()
                .unwrap()
                .count_entities(),
            11
        );

        assert_eq!(
            "seq0:+:5-5"
                .parse::<Interval<Base>>()
                .unwrap()
                .count_entities(),
            1
        );
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
    fn contig_with_colons_round_trips() -> Result<()> {
        let start = Coordinate::<Interbase>::try_new("assembly:chr:1", "+", 10)?;
        let end = Coordinate::<Interbase>::try_new("assembly:chr:1", "+", 20)?;
        let interval = Interval::try_from((start, end))?;
        let parsed = interval.to_string().parse::<Interval<Interbase>>()?;

        assert_eq!(parsed, interval);

        Ok(())
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
        // SAFETY: the endpoints share context and follow the positive strand.
        let interval = Interval::try_from((start, end)).unwrap();

        assert_eq!(interval.to_string(), "seq0:+:0-10");

        // Negative-stranded interval
        let start = "seq0:-:10".parse::<Coordinate<Interbase>>().unwrap();
        let end = "seq0:-:0".parse::<Coordinate<Interbase>>().unwrap();
        // SAFETY: the endpoints share context and follow the negative strand.
        let interval = Interval::try_from((start, end)).unwrap();

        assert_eq!(interval.to_string(), "seq0:-:10-0");
    }
}
