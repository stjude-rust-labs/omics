//! 1-based [`Position`](crate::Position)s.

use crate::position;
use crate::position::Error;
use crate::position::ParseError;
use crate::position::Result;
use crate::position::Value;
use crate::system::One;

mod addition;
mod subtraction;

/// A 1-based, base [`Position`](crate::Position).
pub type Position = crate::Position<One>;

impl position::r#trait::Position<One> for Position {
    /// Creates a new [`Position`].
    ///
    /// Examples
    ///
    /// ```
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::position::one::Position;
    ///
    /// // Non-zero usize
    ///
    /// let value = Value::Usize(1);
    /// let position = Position::try_new(value)?;
    ///
    /// // Zero usize
    ///
    /// let value = Value::Usize(0);
    /// let err = Position::try_new(value).unwrap_err();
    /// assert_eq!(
    ///     err.to_string(),
    ///     "parse error: incompatible value for 1-based, fully-closed coordinate system: 0"
    /// );
    ///
    /// // Lower bound
    ///
    /// let value = Value::LowerBound;
    /// let err = Position::try_new(value).unwrap_err();
    /// assert_eq!(
    ///     err.to_string(),
    ///     "parse error: incompatible value for 1-based, fully-closed coordinate system: ["
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_new(value: impl Into<Value>) -> Result<Self> {
        let value = match value.into() {
            Value::LowerBound => {
                return Err(Error::Parse(ParseError::IncompatibleValue(
                    One.to_string(),
                    Value::LowerBound.to_string(),
                )));
            }
            Value::Usize(0) => {
                return Err(Error::Parse(ParseError::IncompatibleValue(
                    One.to_string(),
                    Value::Usize(0).to_string(),
                )));
            }
            v => v,
        };

        Ok(Self {
            system: One,
            inner: value,
        })
    }
}

impl std::str::FromStr for Position {
    type Err = Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        let value = s
            .parse::<Value>()
            .map_err(|err| Error::Parse(ParseError::ValueError(err)))?;

        Self::try_new(value)
    }
}

impl TryFrom<usize> for Position {
    type Error = Error;

    fn try_from(value: usize) -> std::result::Result<Self, Self::Error> {
        Self::try_new(Value::Usize(value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_remains_the_size_of_the_inner_value() {
        assert_eq!(
            std::mem::size_of::<Position>(),
            std::mem::size_of::<Value>()
        )
    }
}
