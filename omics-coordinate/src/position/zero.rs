//! 0-based, interbase positions.

use crate::position;
use crate::position::value::Number;
use crate::position::Error;
use crate::position::ParseError;
use crate::position::Result;
use crate::position::Value;
use crate::system::Zero;

mod addition;
mod subtraction;

/// A 0-based, interbase position.
///
/// For a more in-depth discussion on what positions are and the notations used
/// within this crate, please see [this section of the docs](crate#positions).
pub type Position = crate::Position<Zero>;

impl position::r#trait::Position<Zero> for Position {
    /// Attempts to creates a new 0-based, interbase position from a [`Value`].
    fn try_new(value: impl Into<Value>) -> Result<Self> {
        Ok(Self {
            system: Zero,
            inner: value.into(),
        })
    }
}

impl Position {
    /// Creates a lower-bound [`Position`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Position;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::Zero;
    ///
    /// // 0-based position
    ///
    /// let position = Position::<Zero>::lower_bound();
    /// assert_eq!(position.inner(), &Value::LowerBound);
    /// ```
    pub fn lower_bound() -> Self {
        // SAFETY: a zero-based position accepts any [`Value`], so this will
        // never fail and can be safely unwrapped.
        Self::try_new(Value::lower_bound()).unwrap()
    }
}

impl TryFrom<Number> for Position {
    type Error = Error;

    fn try_from(value: Number) -> std::result::Result<Self, Self::Error> {
        let value = Value::try_new(value).ok_or(Error::InvalidNumbericalPosition(value))?;
        Self::try_new(value)
    }
}

impl std::str::FromStr for Position {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let value = s
            .parse::<Value>()
            .map_err(|err| Error::Parse(ParseError::Value(err)))?;

        Self::try_new(value)
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
        let position = Position::try_from(10).unwrap();
        assert_eq!(position.inner(), &Value::try_new(10).unwrap());
    }

    #[test]
    fn test_negative_position_creation() {
        let position = Position::lower_bound();
        assert_eq!(position.inner(), &Value::lower_bound());
    }

    #[test]
    fn test_ordering_of_positions() {
        // Ordering of two zero-based positions
        let a = Position::try_from(10).unwrap();
        let b = Position::try_from(5).unwrap();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Greater);

        // Ordering of a zero-based position and a lower bound
        let a = Position::try_from(0).unwrap();
        let b = Position::lower_bound();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Greater);

        // Ordering of a lower bound and a zero-based position
        let a = Position::lower_bound();
        let b = Position::try_from(0).unwrap();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Less);

        // Ordering of two lower bounds
        let a = Position::lower_bound();
        let b = Position::lower_bound();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Equal));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Equal);
    }

    #[test]
    fn test_number_into_position() {
        let position: Position = 30u32.try_into().unwrap();
        assert_eq!(position.inner(), &Value::try_new(30).unwrap());
    }

    #[test]
    fn it_remains_the_size_of_the_inner_value() {
        assert_eq!(
            std::mem::size_of::<Position>(),
            std::mem::size_of::<Value>()
        )
    }
}
