//! 1-based, in-base [`Position`](crate::Position)s.

use crate::position;
use crate::position::value::Kind;
use crate::position::value::Number;
use crate::position::Error;
use crate::position::ParseError;
use crate::position::Result;
use crate::position::Value;
use crate::system::One;

mod addition;
mod subtraction;

/// A 1-based, in-base [`Position`](crate::Position).
///
/// For a more in-depth discussion on what positions are and the notations used
/// within this crate, please see [this section of the docs](crate#positions).
pub type Position = crate::Position<One>;

impl position::r#trait::Position<One> for Position {
    /// Attempts to creates a new 1-based, in-base position from a [`Value`].
    fn try_new(value: impl Into<Value>) -> Result<Self> {
        let value = value.into();

        match value.kind() {
            Kind::LowerBound => {
                return Err(Error::Parse(ParseError::IncompatibleValue(
                    One.to_string(),
                    value.to_string(),
                )));
            }
            Kind::Numerical => {
                if let Some(0) = value.get() {
                    return Err(Error::Parse(ParseError::IncompatibleValue(
                        One.to_string(),
                        value.to_string(),
                    )));
                }
            }
        };

        Ok(Self {
            system: One,
            inner: value,
        })
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

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        let value = s
            .parse::<Value>()
            .map_err(|err| Error::Parse(ParseError::Value(err)))?;

        Self::try_new(value)
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
