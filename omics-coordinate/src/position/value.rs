//! Values for a [`Position`](super::Position).

/// The serialized character that represents a lower bound.
const LOWER_BOUND_CHAR: &str = "[";

/// A error related to parsing a [`Value`].
#[derive(Debug)]
pub enum ParseError {
    /// An invalid value was attempted to be parsed.
    InvalidValue(String),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InvalidValue(v) => write!(f, "invalid value: {v}"),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a [`Value`].
#[derive(Debug)]
pub enum Error {
    /// A parse error.
    ParseError(ParseError),

    /// Attempted to convert a [`Value::LowerBound`] to a [`usize`].
    CannotConvertLowerBoundToUsize,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::ParseError(err) => write!(f, "parse error: {err}"),
            Error::CannotConvertLowerBoundToUsize => {
                write!(f, "lower bound cannot be converted to a usize")
            }
        }
    }
}

impl std::error::Error for Error {}

/// A value for a [`Position`](super::Position).
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Value {
    /// The lower bound.
    LowerBound,

    /// A [`usize`].
    Usize(usize),
}

impl Value {
    /// Attempts to get the value of the position as a [`usize`].
    ///
    /// - If the [`Position`](super::Position) has a [`Value::Usize`], a
    ///   [`Some(usize)`] is returned.
    /// - If the [`Position`](super::Position) has a [`Value::LowerBound`],
    ///   [`None`] is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Zero;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Zero>>()?;
    /// assert_eq!(coordinate.position().inner().get(), Some(1));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn get(&self) -> Option<usize> {
        match self {
            Value::Usize(v) => Some(*v),
            Value::LowerBound => None,
        }
    }
}

impl std::fmt::Display for Value {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Value::Usize(v) => write!(f, "{v}"),
            Value::LowerBound => write!(f, "{}", LOWER_BOUND_CHAR),
        }
    }
}

impl std::str::FromStr for Value {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == LOWER_BOUND_CHAR {
            return Ok(Value::LowerBound);
        }

        let value = s
            .parse::<usize>()
            .map_err(|_| Error::ParseError(ParseError::InvalidValue(s.to_owned())))?;

        Ok(Value::Usize(value))
    }
}

impl TryFrom<Value> for usize {
    type Error = Error;

    fn try_from(value: Value) -> Result<Self, Self::Error> {
        match value {
            Value::Usize(v) => Ok(v),
            Value::LowerBound => Err(Error::CannotConvertLowerBoundToUsize),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn it_deserializes_correctly() -> Result<(), Box<dyn std::error::Error>> {
        let value = "1".parse::<Value>()?;
        assert!(matches!(value, Value::Usize(1)));

        Ok(())
    }

    #[test]
    fn it_converts_a_usize_value_to_a_usize() -> Result<(), Box<dyn std::error::Error>> {
        let value = Value::Usize(1);
        assert_eq!(<Value as TryInto<usize>>::try_into(value)?, 1usize);

        Ok(())
    }

    #[test]
    fn it_errors_when_converting_a_lower_bound_to_a_usize() -> Result<(), Box<dyn std::error::Error>>
    {
        let value = Value::LowerBound;
        let err = <Value as TryInto<usize>>::try_into(value).unwrap_err();
        assert!(matches!(err, Error::CannotConvertLowerBoundToUsize));

        Ok(())
    }

    #[test]
    fn it_orders_values() -> Result<(), Box<dyn std::error::Error>> {
        assert!(Value::LowerBound < Value::Usize(0));
        assert!(Value::Usize(0) < Value::Usize(1));

        Ok(())
    }
}
