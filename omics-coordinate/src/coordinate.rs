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

/// A [`Result`](std::result::Result) with an [`Error`].
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
        fn try_new(
            contig: impl TryInto<Contig, Error = contig::Error>,
            strand: impl TryInto<Strand, Error = strand::Error>,
            position: Number,
        ) -> Result<Self>;
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Coordinate
////////////////////////////////////////////////////////////////////////////////////////

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
    pub fn try_new(
        contig: impl TryInto<Contig, Error = contig::Error>,
        strand: impl TryInto<Strand, Error = strand::Error>,
        position: Number,
    ) -> Result<Self>
    where
        Self: r#trait::Coordinate<S>,
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
    /// assert_eq!(coordinate.into_contig().into_inner(), String::from("seq0"));
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
    /// assert_eq!(contig.into_inner(), String::from("seq0"));
    /// assert_eq!(strand, Strand::Positive);
    /// assert_eq!(position.get(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_parts(self) -> (Contig, Strand, Position<S>) {
        (self.contig, self.strand, self.position)
    }

    /// Consumes `self` and attempts to move the position forward by
    /// `magnitude`.
    ///
    /// This method is dependent on the strand of the coordinate:
    ///
    /// - a coordinate on the [`Strand::Positive`] moves positively, and
    /// - a coordinate on the [`Strand::Negative`] move negatively.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Position;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Number;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// // Interbase.
    ///
    /// let start = "seq0:+:0".parse::<Coordinate<Interbase>>()?;
    /// let coordinate = start.clone().move_forward(10).expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), 10);
    ///
    /// let coordinate = start.move_forward(Number::MAX).expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), Number::MAX);
    ///
    /// let coordinate = coordinate.move_forward(1);
    /// assert!(coordinate.is_none());
    ///
    /// // Base.
    ///
    /// let start = "seq0:+:1".parse::<Coordinate<Base>>()?;
    /// let coordinate = start.clone().move_forward(10).expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), 11);
    ///
    /// let coordinate = start
    ///     .move_forward(Number::MAX - 1)
    ///     .expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), Number::MAX);
    ///
    /// let coordinate = coordinate.move_forward(1);
    /// assert!(coordinate.is_none());
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new coordinate"]
    pub fn move_forward(self, magnitude: Number) -> Option<Coordinate<S>> {
        if magnitude == 0 {
            return Some(self);
        }

        match self.strand {
            Strand::Positive => self.position.checked_add(magnitude),
            Strand::Negative => self.position.checked_sub(magnitude),
        }
        .map(|position| Self::new(self.contig.clone(), self.strand, position))
    }

    /// Consumes `self` and attempts to move the position backwards by
    /// `magnitude`.
    ///
    /// This method is dependent on the strand of the coordinate:
    ///
    /// - a coordinate on the [`Strand::Positive`] moves negatively, and
    /// - a coordinate on the [`Strand::Negative`] move positively.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Contig;
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Position;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::position::Number;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let value = format!("seq0:+:{}", Number::MAX);
    ///
    /// // Interbase.
    ///
    /// let start = value.clone().parse::<Coordinate<Interbase>>()?;
    /// let coordinate = start.clone().move_backward(10).expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), Number::MAX - 10);
    ///
    /// let coordinate = start
    ///     .move_backward(Number::MAX)
    ///     .expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), 0);
    ///
    /// let coordinate = coordinate.move_backward(1);
    /// assert!(coordinate.is_none());
    ///
    /// // Base.
    ///
    /// let start = value.parse::<Coordinate<Base>>()?;
    /// let coordinate = start.clone().move_backward(10).expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), Number::MAX - 10);
    ///
    /// let coordinate = start
    ///     .move_backward(Number::MAX - 1)
    ///     .expect("coordinate to move");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.strand(), Strand::Positive);
    /// assert_eq!(coordinate.position().get(), 1);
    ///
    /// let coordinate = coordinate.move_backward(1);
    /// assert!(coordinate.is_none());
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use = "this method returns a new coordinate"]
    pub fn move_backward(self, magnitude: Number) -> Option<Coordinate<S>> {
        if magnitude == 0 {
            return Some(self);
        }

        match self.strand {
            Strand::Positive => self.position.checked_sub(magnitude),
            Strand::Negative => self.position.checked_add(magnitude),
        }
        .map(|position| Self::new(self.contig.clone(), self.strand, position))
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

impl<S: System> std::str::FromStr for Coordinate<S>
where
    Position<S>: position::r#trait::Position<S>,
{
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts = s.split(VARIANT_SEPARATOR).collect::<Vec<_>>();

        if parts.len() != 3 {
            return Err(Error::Parse(ParseError::Format {
                value: s.to_owned(),
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
        let position = parts
            .next()
            .unwrap()
            .parse::<Position<S>>()
            .map_err(Error::Position)?;

        Ok(Self::new(contig, strand, position))
    }
}
