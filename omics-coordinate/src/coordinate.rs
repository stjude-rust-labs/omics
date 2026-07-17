//! Coordinates.

use omics_core::VARIANT_SEPARATOR;
use thiserror::Error;

use crate::Contig;
use crate::Position;
use crate::Strand;
use crate::System;
use crate::contig;
use crate::position;
use crate::position::Number;
use crate::strand;

pub mod base;
pub mod interbase;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// A parsing error related to a coordinate.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// An invalid format was encountered.
    ///
    /// This generally occurs when an incorrect number of colons (`:`) is used.
    #[error("invalid coordinate format: {value}")]
    Format {
        /// The value that was passed.
        value: String,
    },
}

/// A [`Result`](std::result::Result) with a [`ParseError`].
pub type ParseResult<T> = std::result::Result<T, ParseError>;

/// An error related to a coordinate.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// A contig error.
    #[error("contig error: {0}")]
    Contig(#[from] contig::Error),

    /// A strand error.
    #[error("strand error: {0}")]
    Strand(#[from] strand::Error),

    /// A parse error.
    #[error("parse error: {0}")]
    Parse(#[from] ParseError),

    /// A position error.
    #[error("position error: {0}")]
    Position(#[from] position::Error),
}

/// A [`Result`](std::result::Result) with an [`Error`](enum@Error).
pub type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// The `Coordinate` trait
////////////////////////////////////////////////////////////////////////////////////////

/// Traits related to a coordinate.
pub mod r#trait {
    use super::*;

    /// Requirements to be a coordinate.
    pub trait Coordinate<S: System>:
        std::fmt::Display
        + std::fmt::Debug
        + PartialEq
        + Eq
        + PartialOrd
        + Ord
        + std::str::FromStr<Err = Error>
    where
        Self: Sized,
    {
        /// Attempts to create a new coordinate.
        fn try_new<C, T>(contig: C, strand: T, position: Number) -> Result<Self>
        where
            C: TryInto<Contig>,
            C::Error: Into<contig::Error>,
            T: TryInto<Strand>,
            T::Error: Into<strand::Error>;
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Coordinate
////////////////////////////////////////////////////////////////////////////////////////

/// A borrowed coordinate.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CoordinateRef<'a, S: System> {
    /// The contig.
    contig: &'a Contig,

    /// The strand.
    strand: Strand,

    /// The position.
    position: &'a Position<S>,
}

impl<'a, S: System> CoordinateRef<'a, S> {
    /// Creates a borrowed coordinate.
    pub(crate) fn new(contig: &'a Contig, strand: Strand, position: &'a Position<S>) -> Self {
        Self {
            contig,
            strand,
            position,
        }
    }

    /// Gets the contig.
    pub fn contig(&self) -> &'a Contig {
        self.contig
    }

    /// Gets the strand.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Gets the position.
    pub fn position(&self) -> &'a Position<S> {
        self.position
    }

    /// Returns an owned coordinate.
    pub fn into_owned(self) -> Coordinate<S> {
        self.into()
    }
}

impl<S: System> From<CoordinateRef<'_, S>> for Coordinate<S> {
    fn from(value: CoordinateRef<'_, S>) -> Self {
        Self {
            system: Default::default(),
            contig: value.contig.clone(),
            strand: value.strand,
            position: value.position.clone(),
        }
    }
}

impl<'a, S: System> From<&'a Coordinate<S>> for CoordinateRef<'a, S> {
    fn from(value: &'a Coordinate<S>) -> Self {
        Self::new(&value.contig, value.strand, &value.position)
    }
}

impl<'a, 'b, S: System> From<&'a &'b Coordinate<S>> for CoordinateRef<'b, S> {
    fn from(value: &'a &'b Coordinate<S>) -> Self {
        (*value).into()
    }
}

impl<'a, S: System> From<&'a mut Coordinate<S>> for CoordinateRef<'a, S> {
    fn from(value: &'a mut Coordinate<S>) -> Self {
        Self::new(&value.contig, value.strand, &value.position)
    }
}

impl<S: System> PartialEq<Coordinate<S>> for CoordinateRef<'_, S> {
    fn eq(&self, other: &Coordinate<S>) -> bool {
        self.contig == &other.contig
            && self.strand == other.strand
            && self.position == &other.position
    }
}

impl<S: System> PartialEq<CoordinateRef<'_, S>> for Coordinate<S> {
    fn eq(&self, other: &CoordinateRef<'_, S>) -> bool {
        other == self
    }
}

impl<S: System> std::fmt::Display for CoordinateRef<'_, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}:{}", self.contig, self.strand, self.position)
    }
}

/// A coordinate.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Coordinate<S: System> {
    /// The coordinate system.
    system: S,

    /// The contig.
    contig: Contig,

    /// The strand.
    strand: Strand,

    /// The position.
    position: Position<S>,
}

impl<S: System> Coordinate<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    /// Creates a new coordinate;
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::interbase::Position;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let contig = Contig::new_unchecked("chr1");
    /// let position = Position::new(0);
    /// let strand = Strand::Positive;
    ///
    /// let coordinate = Coordinate::new(contig, strand, position);
    /// ```
    pub fn new(
        contig: impl Into<Contig>,
        strand: impl Into<Strand>,
        position: impl Into<Position<S>>,
    ) -> Self {
        let contig = contig.into();
        let strand = strand.into();
        let position = position.into();

        Self {
            system: Default::default(),
            contig,
            strand,
            position,
        }
    }

    /// Attempts to create a new coordinate.
    pub fn try_new<C, T>(contig: C, strand: T, position: Number) -> Result<Self>
    where
        Self: r#trait::Coordinate<S>,
        C: TryInto<Contig>,
        C::Error: Into<contig::Error>,
        T: TryInto<Strand>,
        T::Error: Into<strand::Error>,
    {
        <Self as r#trait::Coordinate<S>>::try_new(contig, strand, position)
    }

    /// Gets the contig for this coordinate by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        &self.contig
    }

    /// Consumes `self` and returns the inner contig from this
    /// coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(coordinate.into_contig().to_string(), String::from("seq0"));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_contig(self) -> Contig {
        self.contig
    }

    /// Gets the strand for this coordinate by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Gets the position for this coordinate by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(coordinate.position().get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn position(&self) -> &Position<S> {
        &self.position
    }

    /// Consumes `self` and returns the inner position from this
    /// coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Interbase>>()?;
    /// assert_eq!(coordinate.into_position().get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_position(self) -> Position<S> {
        self.position
    }

    /// Consumes `self` to return the parts that comprise this coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Interbase>>()?;
    ///
    /// let (contig, strand, position) = coordinate.into_parts();
    /// assert_eq!(contig.to_string(), String::from("seq0"));
    /// assert_eq!(strand, Strand::Positive);
    /// assert_eq!(position.get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_parts(self) -> (Contig, Strand, Position<S>) {
        (self.contig, self.strand, self.position)
    }

    /// Attempts to move the position forward by `magnitude` in place.
    ///
    /// This method is dependent on the strand of the coordinate:
    ///
    /// - a coordinate on the [`Strand::Positive`] moves positively, and
    /// - a coordinate on the [`Strand::Negative`] moves negatively.
    ///
    /// Returns `true` if the move succeeded, `false` if it would overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Number;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// // Interbase.
    ///
    /// let mut coordinate = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// assert!(coordinate.move_forward(10));
    /// assert_eq!(coordinate.position().get(), 10);
    ///
    /// let mut coordinate = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// assert!(coordinate.move_forward(Number::MAX));
    /// assert_eq!(coordinate.position().get(), Number::MAX);
    /// assert!(!coordinate.move_forward(1));
    ///
    /// // Base.
    ///
    /// let mut coordinate = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// assert!(coordinate.move_forward(10));
    /// assert_eq!(coordinate.position().get(), 11);
    ///
    /// let mut coordinate = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// assert!(coordinate.move_forward(Number::MAX - 1));
    /// assert_eq!(coordinate.position().get(), Number::MAX);
    /// assert!(!coordinate.move_forward(1));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_forward(&mut self, magnitude: Number) -> bool {
        if magnitude == 0 {
            return true;
        }

        let result = match self.strand {
            Strand::Positive => self.position.checked_add(magnitude),
            Strand::Negative => self.position.checked_sub(magnitude),
        };

        match result {
            Some(position) => {
                self.position = position;
                true
            }
            None => false,
        }
    }

    /// Consumes `self` and attempts to move the position forward by
    /// `magnitude`.
    ///
    /// This method is dependent on the strand of the coordinate:
    ///
    /// - a coordinate on the [`Strand::Positive`] moves positively, and
    /// - a coordinate on the [`Strand::Negative`] moves negatively.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Number;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// // Interbase.
    ///
    /// let coordinate = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// let moved = coordinate
    ///     .into_move_forward(10)
    ///     .expect("coordinate to move");
    /// assert_eq!(moved.position().get(), 10);
    ///
    /// let coordinate = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// let moved = coordinate.into_move_forward(Number::MAX).unwrap();
    /// assert!(moved.into_move_forward(1).is_none());
    ///
    /// // Base.
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// let moved = coordinate
    ///     .into_move_forward(10)
    ///     .expect("coordinate to move");
    /// assert_eq!(moved.position().get(), 11);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new coordinate"]
    pub fn into_move_forward(self, magnitude: Number) -> Option<Coordinate<S>> {
        if magnitude == 0 {
            return Some(self);
        }

        match self.strand {
            Strand::Positive => self.position.checked_add(magnitude),
            Strand::Negative => self.position.checked_sub(magnitude),
        }
        .map(|position| Self::new(self.contig, self.strand, position))
    }

    /// Attempts to move the position backward by `magnitude` in place.
    ///
    /// This method is dependent on the strand of the coordinate:
    ///
    /// - a coordinate on the [`Strand::Positive`] moves negatively, and
    /// - a coordinate on the [`Strand::Negative`] moves positively.
    ///
    /// Returns `true` if the move succeeded, `false` if it would overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Number;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let value = format!("seq0:+:{}", Number::MAX);
    ///
    /// // Interbase.
    ///
    /// let mut coordinate = value.clone().parse::<Coordinate<Interbase>>()?;
    /// assert!(coordinate.move_backward(10));
    /// assert_eq!(coordinate.position().get(), Number::MAX - 10);
    ///
    /// let mut coordinate = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// assert!(!coordinate.move_backward(1));
    ///
    /// // Base.
    ///
    /// let mut coordinate = value.parse::<Coordinate<Base>>()?;
    /// assert!(coordinate.move_backward(10));
    /// assert_eq!(coordinate.position().get(), Number::MAX - 10);
    ///
    /// let mut coordinate = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// assert!(!coordinate.move_backward(1));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_backward(&mut self, magnitude: Number) -> bool {
        if magnitude == 0 {
            return true;
        }

        let result = match self.strand {
            Strand::Positive => self.position.checked_sub(magnitude),
            Strand::Negative => self.position.checked_add(magnitude),
        };

        match result {
            Some(position) => {
                self.position = position;
                true
            }
            None => false,
        }
    }

    /// Consumes `self` and attempts to move the position backwards by
    /// `magnitude`.
    ///
    /// This method is dependent on the strand of the coordinate:
    ///
    /// - a coordinate on the [`Strand::Positive`] moves negatively, and
    /// - a coordinate on the [`Strand::Negative`] moves positively.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Number;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let value = format!("seq0:+:{}", Number::MAX);
    ///
    /// // Interbase.
    ///
    /// let coordinate = value.clone().parse::<Coordinate<Interbase>>()?;
    /// let moved = coordinate
    ///     .into_move_backward(10)
    ///     .expect("coordinate to move");
    /// assert_eq!(moved.position().get(), Number::MAX - 10);
    ///
    /// let coordinate = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// assert!(coordinate.into_move_backward(1).is_none());
    ///
    /// // Base.
    ///
    /// let coordinate = value.parse::<Coordinate<Base>>()?;
    /// let moved = coordinate
    ///     .into_move_backward(10)
    ///     .expect("coordinate to move");
    /// assert_eq!(moved.position().get(), Number::MAX - 10);
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// assert!(coordinate.into_move_backward(1).is_none());
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new coordinate"]
    pub fn into_move_backward(self, magnitude: Number) -> Option<Coordinate<S>> {
        if magnitude == 0 {
            return Some(self);
        }

        match self.strand {
            Strand::Positive => self.position.checked_sub(magnitude),
            Strand::Negative => self.position.checked_add(magnitude),
        }
        .map(|position| Self::new(self.contig, self.strand, position))
    }

    /// Swaps the strand of the coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// //===========//
    /// // Interbase //
    /// //===========//
    ///
    /// let coordinate = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    /// let swapped = coordinate.swap_strand();
    /// assert_eq!(swapped.contig().as_str(), "seq0");
    /// assert_eq!(swapped.strand(), Strand::Negative);
    /// assert_eq!(swapped.position().get(), 10);
    ///
    /// //======//
    /// // Base //
    /// //======//
    ///
    /// let coordinate = Coordinate::<Base>::try_new("seq0", "-", 10)?;
    /// let swapped = coordinate.swap_strand();
    /// assert_eq!(swapped.contig().as_str(), "seq0");
    /// assert_eq!(swapped.strand(), Strand::Positive);
    /// assert_eq!(swapped.position().get(), 10);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new coordinate"]
    pub fn swap_strand(self) -> Coordinate<S> {
        let (contig, strand, position) = self.into_parts();
        Coordinate::new(contig, strand.complement(), position)
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl<S: System> std::fmt::Display for Coordinate<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if !f.alternate() {
            write!(f, "{}:{}:{}", self.contig, self.strand, self.position)
        } else {
            write!(
                f,
                "{}:{}:{} ({})",
                self.contig, self.strand, self.position, self.system
            )
        }
    }
}

impl<S: System> TryFrom<&str> for Coordinate<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    type Error = Error;

    /// Parses a coordinate from its `contig:strand:position` string form.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Base;
    ///
    /// let coordinate = Coordinate::<Base>::try_from("seq0:+:1")?;
    /// assert_eq!(coordinate.position().get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_from(s: &str) -> Result<Self> {
        s.parse()
    }
}

impl<S: System> std::str::FromStr for Coordinate<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let (prefix, position) =
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
        let position = position.parse::<Position<S>>().map_err(Error::Position)?;

        Ok(Self::new(contig, strand, position))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::Interbase;

    #[test]
    fn contig_with_colons_round_trips() -> Result<()> {
        let coordinate = Coordinate::<Interbase>::try_new("assembly:chr:1", "+", 10)?;
        let parsed = coordinate.to_string().parse::<Coordinate<Interbase>>()?;

        assert_eq!(parsed, coordinate);

        Ok(())
    }

    #[test]
    fn try_new_accepts_validated_parts() -> Result<()> {
        let contig = Contig::try_new("seq0")?;
        let coordinate = Coordinate::<Interbase>::try_new(contig.clone(), Strand::Positive, 10)?;

        assert_eq!(coordinate.contig(), &contig);
        assert_eq!(coordinate.strand(), Strand::Positive);
        assert_eq!(coordinate.position().get(), 10);

        Ok(())
    }
}
