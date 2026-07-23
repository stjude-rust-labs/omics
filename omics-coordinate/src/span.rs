//! Span values.

use std::cmp::Ordering;
use std::cmp::max;
use std::cmp::min;
use std::fmt;
use std::str::FromStr;

use thiserror::Error;

use crate::Position;
use crate::System;
use crate::position;
use crate::position::Number;
use crate::system::Base;
use crate::system::Interbase;

pub mod base;
pub mod interbase;

/// The separator between a span start and end position.
const SPAN_SEPARATOR: &str = "-";

/// A direction along a coordinate system.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Direction {
    /// The end lies after the start.
    Ascending,

    /// The endpoints are equal.
    Stationary,

    /// The end lies before the start.
    Descending,
}

impl Direction {
    /// Returns whether two directions can be clamped together.
    pub(crate) fn is_compatible(self, other: Self) -> bool {
        self == other || self == Self::Stationary || other == Self::Stationary
    }
}

impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value = match self {
            Self::Ascending => "ascending",
            Self::Stationary => "stationary",
            Self::Descending => "descending",
        };

        f.write_str(value)
    }
}

/// An error from clamping one span by another.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ClampError {
    /// The spans move in incompatible directions.
    #[error(
        "direction mismatch between spans `{original_start}-{original_end}` and \
         `{operand_start}-{operand_end}`; directions are `{original_direction}` and \
         `{operand_direction}`"
    )]
    DirectionMismatch {
        /// The start position of the original span.
        original_start: Number,

        /// The end position of the original span.
        original_end: Number,

        /// The direction of the original span.
        original_direction: Direction,

        /// The start position of the operand span.
        operand_start: Number,

        /// The end position of the operand span.
        operand_end: Number,

        /// The direction of the operand span.
        operand_direction: Direction,
    },

    /// The spans do not overlap.
    #[error("disjoint spans `{original_start}-{original_end}` and `{operand_start}-{operand_end}`")]
    Disjoint {
        /// The start position of the original span.
        original_start: Number,

        /// The end position of the original span.
        original_end: Number,

        /// The start position of the operand span.
        operand_start: Number,

        /// The end position of the operand span.
        operand_end: Number,
    },
}

/// A result with a span clamping error.
pub type ClampResult<T> = std::result::Result<T, ClampError>;

/// An error from building a span from raw numbers.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// The start position is invalid.
    #[error("invalid span start position")]
    Start(#[source] position::Error),

    /// The end position is invalid.
    #[error("invalid span end position")]
    End(#[source] position::Error),
}

/// A result with a span construction error.
pub type Result<T> = std::result::Result<T, Error>;

/// An error from parsing a span string.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// The string does not have the expected shape.
    #[error("invalid span format `{0}`")]
    Format(String),

    /// The start position is invalid.
    #[error("invalid span start position")]
    Start(#[source] position::Error),

    /// The end position is invalid.
    #[error("invalid span end position")]
    End(#[source] position::Error),
}

/// A result with a span parsing error.
pub type ParseResult<T> = std::result::Result<T, ParseError>;

/// Internal trait definitions for span operations.
mod r#trait {
    use super::*;

    /// Trait for creating positions from numeric values within coordinate
    /// systems.
    pub trait PositionForSpan<S: System> {
        /// Creates a new position from a numeric value.
        fn try_new_position(value: Number) -> std::result::Result<Position<S>, position::Error>;
    }

    /// Trait for coordinate span operations.
    pub trait Span<S: System> {
        /// Tests whether this span contains an entity at the given position.
        fn contains_entity(&self, position: &Position<Base>) -> bool;

        /// Returns the count of entities in this span.
        fn count_entities(&self) -> Number;

        /// Tests whether this span contains no entities.
        fn is_empty(&self) -> bool;
    }
}

impl r#trait::PositionForSpan<Base> for Position<Base> {
    fn try_new_position(value: Number) -> std::result::Result<Position<Base>, position::Error> {
        Position::<Base>::try_new(value)
    }
}

impl r#trait::PositionForSpan<Interbase> for Position<Interbase> {
    fn try_new_position(
        value: Number,
    ) -> std::result::Result<Position<Interbase>, position::Error> {
        Ok(Position::<Interbase>::new(value))
    }
}

/// A span with a typed coordinate system.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Span<S: System> {
    /// The start position of this span.
    start: Position<S>,
    /// The end position of this span.
    end: Position<S>,
}

impl<S: System> Span<S>
where
    Position<S>: r#trait::PositionForSpan<S>,
{
    /// Creates a span from start and end numbers.
    pub fn try_new(start: Number, end: Number) -> Result<Self> {
        let start = <Position<S> as r#trait::PositionForSpan<S>>::try_new_position(start)
            .map_err(Error::Start)?;
        let end = <Position<S> as r#trait::PositionForSpan<S>>::try_new_position(end)
            .map_err(Error::End)?;

        Ok(Self { start, end })
    }
}

impl<S: System> Span<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    /// Returns the start position.
    pub fn start(&self) -> &Position<S> {
        &self.start
    }

    /// Returns the end position.
    pub fn end(&self) -> &Position<S> {
        &self.end
    }

    /// Consumes the span into both positions.
    pub fn into_positions(self) -> (Position<S>, Position<S>) {
        (self.start, self.end)
    }

    /// Returns the direction from start to end.
    pub fn direction(&self) -> Direction {
        match self.start.cmp(&self.end) {
            Ordering::Less => Direction::Ascending,
            Ordering::Equal => Direction::Stationary,
            Ordering::Greater => Direction::Descending,
        }
    }

    /// Returns a span with swapped endpoints.
    #[must_use = "this method returns a new span"]
    pub fn reversed(self) -> Self {
        let (start, end) = self.into_positions();

        Self {
            start: end,
            end: start,
        }
    }

    /// Returns whether the span contains the given position.
    pub fn contains_position(&self, position: &Position<S>) -> bool {
        match self.direction() {
            Direction::Ascending => self.start <= *position && *position <= self.end,
            Direction::Stationary => self.start == *position,
            Direction::Descending => self.end <= *position && *position <= self.start,
        }
    }

    /// Returns the offset of a position from the start of the span.
    pub fn position_offset(&self, position: &Position<S>) -> Option<Number> {
        self.contains_position(position)
            .then(|| self.start.distance_unchecked(position))
    }

    /// Returns the position at an offset from the start of the span.
    pub fn position_at_offset(&self, offset: Number) -> Option<Position<S>> {
        let position = match self.direction() {
            Direction::Ascending => self.start.checked_add(offset)?,
            Direction::Stationary if offset == 0 => self.start.clone(),
            Direction::Stationary => return None,
            Direction::Descending => self.start.checked_sub(offset)?,
        };

        self.contains_position(&position).then_some(position)
    }

    /// Clamps this span by another span.
    #[must_use = "this method returns a new span"]
    pub fn clamp(self, operand: Self) -> ClampResult<Self> {
        let original_direction = self.direction();
        let operand_direction = operand.direction();
        let original_start = self.start.get();
        let original_end = self.end.get();
        let operand_start = operand.start.get();
        let operand_end = operand.end.get();

        if !original_direction.is_compatible(operand_direction) {
            return Err(ClampError::DirectionMismatch {
                original_start,
                original_end,
                original_direction,
                operand_start,
                operand_end,
                operand_direction,
            });
        }

        let lower = max(
            min(self.start.clone(), self.end.clone()),
            min(operand.start.clone(), operand.end.clone()),
        );
        let upper = min(max(self.start, self.end), max(operand.start, operand.end));

        if lower > upper {
            return Err(ClampError::Disjoint {
                original_start,
                original_end,
                operand_start,
                operand_end,
            });
        }

        Ok(match original_direction {
            Direction::Ascending => Self::from((lower, upper)),
            Direction::Stationary => Self::from((lower, upper)),
            Direction::Descending => Self::from((upper, lower)),
        })
    }
}

impl<S: System> Span<S>
where
    Span<S>: r#trait::Span<S>,
    Position<S>: position::r#trait::Position<S>,
{
    /// Returns whether the span contains the given entity.
    pub fn contains_entity(&self, position: &Position<Base>) -> bool {
        <Self as r#trait::Span<S>>::contains_entity(self, position)
    }

    /// Counts the entities contained in the span.
    pub fn count_entities(&self) -> Number {
        <Self as r#trait::Span<S>>::count_entities(self)
    }

    /// Returns whether the span contains no entities.
    pub fn is_empty(&self) -> bool {
        <Self as r#trait::Span<S>>::is_empty(self)
    }
}

impl<S: System> From<(Position<S>, Position<S>)> for Span<S> {
    fn from((start, end): (Position<S>, Position<S>)) -> Self {
        Self { start, end }
    }
}

impl<S: System> fmt::Display for Span<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

impl<S: System> FromStr for Span<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    type Err = ParseError;

    fn from_str(s: &str) -> ParseResult<Self> {
        let (start, end) = s
            .split_once(SPAN_SEPARATOR)
            .filter(|(_, end)| !end.contains(SPAN_SEPARATOR))
            .ok_or_else(|| ParseError::Format(s.to_owned()))?;

        let start = start.parse::<Position<S>>().map_err(ParseError::Start)?;
        let end = end.parse::<Position<S>>().map_err(ParseError::End)?;

        Ok(Self { start, end })
    }
}
