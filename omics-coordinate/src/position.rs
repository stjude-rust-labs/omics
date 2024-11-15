//! Positions.

use std::fmt::Debug;

use crate::CheckedAdd;
use crate::CheckedSub;
use crate::position::value::Kind;
use crate::system::System;

pub mod one;
pub mod value;
pub mod zero;

use value::Number;
pub use value::Value;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// A error related to parsing a [`Position`].
#[derive(Debug)]
pub enum ParseError {
    /// Incompatible value for a given coordinate system.
    ///
    /// * The first argument is the name of the coordinate system.
    /// * The second argument is the value.
    IncompatibleValue(String, String),

    /// An invalid value was attempted to be parsed.
    Value(value::Error),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::IncompatibleValue(system, v) => {
                write!(f, "incompatible value for {system}: {v}")
            }
            ParseError::Value(v) => write!(f, "value error: {v}"),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a [`Position`].
#[derive(Debug)]
pub enum Error {
    /// Cannot construct a range for an interval.
    IntervalRange(String),

    /// A parse error.
    Parse(ParseError),

    /// An error with a [`Value`] was encountered.
    Value(value::Error),

    /// Could not create value from [`Number`].
    InvalidValue(Number),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::IntervalRange(interval) => {
                write!(f, "cannot construct range for interval: {interval}")
            }
            Error::Parse(err) => write!(f, "parse error: {err}"),
            Error::Value(err) => write!(f, "value error: {err}"),
            Error::InvalidValue(number) => write!(f, "invalid value: {number}"),
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// Positions
////////////////////////////////////////////////////////////////////////////////////////

/// A nucleotide offset from the start of a nucleic acid molecule.
///
/// For a more in-depth discussion on what positions are and the notations used
/// within this crate, please see [this section of the docs](crate#positions).
#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct Position<S: System> {
    /// The coordinate system.
    system: S,

    /// The value.
    pub(crate) inner: Value,
}

impl<S: System> Position<S> {
    /// Attempts to create a new [`Position`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    /// let value = Value::Usize(1);
    /// let position = Position::<One>::try_new(value)?;
    /// assert_eq!(position.inner().get(), Some(1));
    ///
    /// let value = Value::Usize(0);
    /// let err = Position::<One>::try_new(value).unwrap_err();
    /// assert_eq!(
    ///     err.to_string(),
    ///     "parse error: incompatible value for 1-based, fully-closed coordinate system: 0"
    /// );
    ///
    /// let value = Value::LowerBound;
    /// let err = Position::<One>::try_new(value).unwrap_err();
    /// assert_eq!(
    ///     err.to_string(),
    ///     String::from(
    ///         "parse error: incompatible value for 1-based, fully-closed coordinate system: ["
    ///     )
    /// );
    ///
    /// // 0-based position
    /// let value = Value::Usize(0);
    /// let position = Position::<Zero>::try_new(value)?;
    /// assert_eq!(position.inner().get(), Some(0));
    ///
    /// let value = Value::LowerBound;
    /// let position = Position::<Zero>::try_new(value)?;
    /// assert_eq!(position.inner().get(), None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(value: impl Into<Value>) -> Result<Self>
    where
        Self: r#trait::Position<S>,
    {
        <Self as r#trait::Position<S>>::try_new(value)
    }

    /// Gets the inner [`Value`] by reference.
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    ///
    /// let position = Position::<One>::try_from(1)?;
    /// assert_eq!(position.inner(), &Value::Usize(1));
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::try_from(0)?;
    /// assert_eq!(position.inner(), &Value::Usize(0));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn inner(&self) -> &Value {
        &self.inner
    }

    /// Consumes the [`Position`] and returns the inner [`Value`].
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    ///
    /// let position = Position::<One>::try_from(1)?;
    /// assert_eq!(position.into_inner(), Value::Usize(1));
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::try_from(0)?;
    /// assert_eq!(position.into_inner(), Value::Usize(0));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_inner(self) -> Value {
        self.inner
    }

    /// Gets the value of the [`Position`] as a [`Number`] wrapped in [`Some`]
    /// if the inner value is of type [`Value::Usize`]. Else, returns
    /// [`None`] (notably, for the lower bound).
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    ///
    /// let position = Position::<One>::try_from(1)?;
    /// assert_eq!(position.get(), Some(1));
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::try_from(0)?;
    /// assert_eq!(position.get(), Some(0));
    ///
    /// let position = Position::<Zero>::lower_bound();
    /// assert_eq!(position.get(), None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn get(&self) -> Option<Number> {
        self.inner.get()
    }

    /// Attempts to add the specified magnitude to inner value of the
    /// [`Position`].
    ///
    /// - If the addition occurs succesfully, then the new [`Position`] is
    ///   returned wrapped in [`Some`].
    /// - If the addition overflows, [`None`] is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    ///
    /// let position = Position::<One>::try_from(1)?;
    /// let result = position.checked_add(1);
    /// assert_eq!(result, Some(Position::<One>::try_from(2)?));
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::try_from(0)?;
    /// let result = position.checked_add(1);
    /// assert_eq!(result, Some(Position::<Zero>::try_from(1)?));
    ///
    /// let position = Position::<Zero>::lower_bound();
    /// let result = position.checked_add(1);
    /// assert_eq!(result, Some(Position::<Zero>::try_from(0)?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn checked_add<T>(&self, rhs: T) -> Option<Self>
    where
        Self: crate::CheckedAdd<T, Output = Self>,
    {
        <Self as CheckedAdd<T>>::checked_add(self, rhs)
    }

    /// Attempts to subtract the specified magnitude from inner value of the
    /// [`Position`].
    ///
    /// - If the subtraction occurs succesfully, then [`Some<Position>`] is
    ///   returned.
    /// - If the subtraction would overflow, [`None`] is returned.
    ///
    /// # Note
    ///
    /// This method does not take into account what strand the [`Position`]
    /// falls on. As such, it's probably not what you want to use unless you
    /// are working with very low-level code.
    ///
    /// If you're trying to move a [`Coordinate`](crate::Coordinate) backwards,
    /// use [`Coordinate::move_backward`](crate::Coordinate::move_backward)
    /// instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::One;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    ///
    /// let position = Position::<One>::try_from(2)?;
    /// let result = position.checked_sub(1);
    /// assert_eq!(result, Some(Position::<One>::try_from(1)?));
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::try_from(1)?;
    /// let result = position.checked_sub(1);
    /// assert_eq!(result, Some(Position::<Zero>::try_from(0)?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn checked_sub<T>(&self, rhs: T) -> Option<Self>
    where
        Self: crate::CheckedSub<T, Output = Self>,
    {
        <Self as CheckedSub<T>>::checked_sub(self, rhs)
    }

    /// Gets the magnitude of the distance between two [`Position<S>`]s.
    ///
    /// **Note**: this method calculates the magnitude of the distance between
    /// two [`Position<S>`]s that are assumed to be on the same strand and
    /// contig. Notably, **there is no check regarding strand or contig
    /// equivalence within this calculation**.
    pub fn distance_unchecked(&self, other: &Position<S>) -> Option<Number> {
        match (self.inner().kind(), other.inner().kind()) {
            (Kind::Numerical, Kind::Numerical) => {
                // SAFETY: because both kinds are numerical, they are required
                // to always unwrap here.
                let a = self.inner().get().unwrap();
                let b = self.inner().get().unwrap();

                if a > b {
                    a.checked_sub(b)
                } else {
                    b.checked_sub(a)
                }
            }
            (Kind::Numerical, Kind::LowerBound) => {
                // SAFETY: because the first member of the tuple is numerical,
                // it must always unwrap here.
                let a = self.inner().get().unwrap();
                a.checked_add(1)
            }
            (Kind::LowerBound, Kind::Numerical) => {
                // SAFETY: because the second member of the tuple is numerical,
                // it must always unwrap here.
                let b = other.inner().get().unwrap();
                b.checked_add(1)
            }
            (Kind::LowerBound, Kind::LowerBound) => Some(0),
        }
    }
}

impl<S: System> std::fmt::Display for Position<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if !f.alternate() {
            write!(f, "{}", self.inner)
        } else {
            write!(f, "{} ({})", self.inner, self.system)
        }
    }
}

impl<S: System> std::fmt::Debug for Position<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Position")
            .field("inner", &self.inner)
            .finish()
    }
}

impl<S: System> TryFrom<Position<S>> for Number {
    type Error = Error;

    fn try_from(value: Position<S>) -> std::result::Result<Self, Self::Error> {
        <Value as TryInto<Number>>::try_into(value.into_inner()).map_err(Error::Value)
    }
}

/// Traits related to a position.
pub mod r#trait {
    use super::*;

    /// Requirements to be a position.
    pub trait Position<S: System>:
        Eq
        + Ord
        + PartialEq
        + PartialOrd
        + std::str::FromStr<Err = Error>
        + CheckedAdd<Number, Output = Self>
        + CheckedSub<Number, Output = Self>
    {
        /// Attempts to create a new [`Position<S>`].
        fn try_new(value: impl Into<Value>) -> Result<Self>;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::One;
    use crate::system::Zero;

    #[test]
    fn it_converts_a_usize_valued_position_to_a_usize() {
        let position: Position<Zero> = 1.into();
        let value = <Position<Zero> as TryInto<Number>>::try_into(position)?;

        assert_eq!(value, 1);
    }

    #[test]
    fn it_fails_to_convert_a_lower_bound_valued_position_to_a_usize() {
        let position = Position::<Zero>::lower_bound();
        let err = <Position<Zero> as TryInto<Number>>::try_into(position).unwrap_err();

        assert!(matches!(err, Error::Value(value::Error::LowerBoundToUsize)));
    }

    #[test]
    fn it_calculates_the_distance_between_two_zero_based_positions() -> Result<()> {
        // Zero distance between two Number positions.
        let a = Position::<Zero>::try_new(10).unwrap();
        let b = Position::<Zero>::try_new(10).unwrap();

        assert_eq!(a.distance_unchecked(&b), Some(0));
        assert_eq!(b.distance_unchecked(&a), Some(0));

        // Non-zero distance between two Number positions.
        let a = Position::<Zero>::try_new(10).unwrap();
        let b = Position::<Zero>::try_new(5).unwrap();

        assert_eq!(a.distance_unchecked(&b), Some(5));
        assert_eq!(b.distance_unchecked(&a), Some(5));

        // distance between two Number positions.
        let a = Position::<Zero>::try_new(Number::MAX).unwrap();
        let b = Position::<Zero>::lower_bound();

        assert_eq!(a.distance_unchecked(&b), None);
        assert_eq!(b.distance_unchecked(&a), None);

        Ok(())
    }

    #[test]
    fn it_calculates_the_distance_between_two_one_based_positions() -> Result<()> {
        // Zero distance between two Number positions.
        let a = Position::<One>::try_from(10)?;
        let b = Position::<One>::try_from(10)?;

        assert_eq!(a.distance_unchecked(&b), Some(0));
        assert_eq!(b.distance_unchecked(&a), Some(0));

        // Non-zero distance between two Number positions.
        let a = Position::<One>::try_from(10)?;
        let b = Position::<One>::try_from(5)?;

        assert_eq!(a.distance_unchecked(&b), Some(5));
        assert_eq!(b.distance_unchecked(&a), Some(5));

        // Distance between two Number positions.
        let a = Position::<One>::try_from(Number::MAX)?;
        let b = Position::<One>::try_from(1)?;

        assert_eq!(a.distance_unchecked(&b), Some(Number::MAX - 1));
        assert_eq!(b.distance_unchecked(&a), Some(Number::MAX - 1));

        Ok(())
    }
}
