//! Span values.

use std::cmp::Ordering;
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

mod r#trait {
    use super::*;

    pub trait PositionForSpan<S: System> {
        fn try_new_position(value: Number) -> std::result::Result<Position<S>, position::Error>;
    }

    pub trait Span<S: System> {
        fn contains_entity(&self, position: &Position<Base>) -> bool;

        fn count_entities(&self) -> Number;

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
    start: Position<S>,
    end: Position<S>,
}

impl<S: System> Span<S>
where
    Position<S>: position::r#trait::Position<S> + r#trait::PositionForSpan<S>,
{
    /// Creates a span from start and end numbers.
    pub fn try_new(start: Number, end: Number) -> Result<Self> {
        let start = <Position<S> as r#trait::PositionForSpan<S>>::try_new_position(start)
            .map_err(Error::Start)?;
        let end = <Position<S> as r#trait::PositionForSpan<S>>::try_new_position(end)
            .map_err(Error::End)?;

        Ok(Self { start, end })
    }

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
