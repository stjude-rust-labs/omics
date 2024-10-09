//! 0-based, half-open [`Position`](crate::Position)s.

use crate::position;
use crate::position::Error;
use crate::position::ParseError;
use crate::position::Result;
use crate::position::Value;
use crate::system::Zero;

mod addition;
mod subtraction;

/// A 0-based, half-open [`Position`](crate::Position).
pub type Position = crate::Position<Zero>;

impl position::r#trait::Position<Zero> for Position {
    /// Creates a new [`Position`].
    ///
    /// Examples
    ///
    /// ```
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::position::zero::Position;
    ///
    /// // Usize
    ///
    /// let value = Value::Usize(0);
    /// let position = Position::try_new(value)?;
    ///
    /// // Lower bound
    ///
    /// let value = Value::LowerBound;
    /// let position = Position::try_new(value)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_new(value: impl Into<Value>) -> Result<Self> {
        Ok(Self {
            system: Zero,
            inner: value.into(),
        })
    }
}

impl Position {
    /// Creates a [`Position`] with a [`Value::LowerBound`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 1-based position
    ///
    /// let position = Position::<Zero>::lower_bound();
    /// assert_eq!(position.inner(), &Value::LowerBound);
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::lower_bound();
    /// assert_eq!(position.inner(), &Value::LowerBound);
    /// ```
    pub fn lower_bound() -> Self {
        Self {
            system: Zero,
            inner: Value::LowerBound,
        }
    }
}

impl std::str::FromStr for Position {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let value = s
            .parse::<Value>()
            .map_err(|err| Error::Parse(ParseError::ValueError(err)))?;

        Self::try_new(value)
    }
}

impl From<usize> for Position {
    fn from(value: usize) -> Self {
        // SAFETY: a zero-based position accepts any usize, so this will never
        // fail and can be safely unwrapped.
        Self::try_new(Value::Usize(value)).unwrap()
    }
}

impl From<Value> for Position {
    fn from(value: Value) -> Self {
        // SAFETY: a zero-based position accepts any [`Value`], so this will
        // never fail and can be safely unwrapped.
        Self::try_new(value).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_based_position_creation() {
        let position = Position::from(10);
        assert_eq!(position.inner(), &Value::Usize(10));
    }

    #[test]
    fn test_negative_position_creation() {
        let position = Position::lower_bound();
        assert_eq!(position.inner(), &Value::LowerBound);
    }

    #[test]
    fn test_ordering_of_positions() {
        // Ordering of two zero-based positions
        let a = Position::from(10);
        let b = Position::from(5);
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Greater);

        // Ordering of a zero-based position and a lower bound
        let a = Position::from(0);
        let b = Position::lower_bound();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Greater);

        // Ordering of a lower bound and a zero-based position
        let a = Position::lower_bound();
        let b = Position::from(0);
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Less);

        // Ordering of two lower bounds
        let a = Position::lower_bound();
        let b = Position::lower_bound();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Equal));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Equal);
    }

    #[test]
    fn test_usize_into_position() {
        let position: Position = 30usize.into();
        assert_eq!(position.inner(), &Value::Usize(30));
    }

    #[test]
    fn it_remains_the_size_of_the_inner_value() {
        assert_eq!(
            std::mem::size_of::<Position>(),
            std::mem::size_of::<Value>()
        )
    }
}
