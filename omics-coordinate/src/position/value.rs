//! Values for a position.

use std::cmp::Ordering;

/// The serialized character that represents a lower-bound.
const LOWER_BOUND_CHAR: &str = "[";

/// The inner representation for a numerical position value.
#[cfg(not(feature = "position-u64"))]
pub type Number = u32;

/// The inner representation for a numerical position value.
#[cfg(feature = "position-u64")]
pub type Number = u64;

/// The special value that represents a lower-bound.
const LOWER_BOUND_VALUE: Number = Number::MAX;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

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

    /// Cannot convert a lower-bound to a [`usize`].
    LowerBoundToUsize,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::ParseError(err) => write!(f, "parse error: {err}"),
            Error::LowerBoundToUsize => {
                write!(f, "lower-bound cannot be converted to a usize")
            }
        }
    }
}

impl std::error::Error for Error {}

////////////////////////////////////////////////////////////////////////////////////////
// Kinds of values
////////////////////////////////////////////////////////////////////////////////////////

/// A kind of value.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Kind {
    /// A lower-bound position.
    ///
    /// The lower bound represents a non-resident value at -1. Obviously,
    /// positions cannot have a value of -1. That being said, intervals that are
    /// non-inclusive of the end position require such a value to ensure that,
    /// when there is an interval on the negative strand, zero can be included
    /// within the interval.
    LowerBound,

    /// A numerical position.
    Numerical,
}

impl std::fmt::Display for Kind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Kind::LowerBound => write!(f, "lower-bound"),
            Kind::Numerical => todo!(),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Values
////////////////////////////////////////////////////////////////////////////////////////

/// A value for a position.
///
/// # Notes
///
/// * The maximum value a numbered position can hold is the maximum of the inner
///   numerical type minus 1. This is because [`Number::MAX`] is used as a
///   special value to represent the lower-bound.
/// * As an end user, you shouldn't have to use this struct because you
///   shouldn't generally be constructing [`Value`]s directly—they are intended
///   to be an implementation detail of [`Position`]s. Instead, you should
///   probably just construct a position in the coordinate system of your
///   choosing directly.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Value(Number);

impl Ord for Value {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.0, other.0) {
            (LOWER_BOUND_VALUE, _) => Ordering::Less,
            (_, LOWER_BOUND_VALUE) => Ordering::Greater,
            // NOTE: the lower-bound to lower-bound check is handled here, so it
            // doesn't need it's own match clause.
            (a, b) => a.cmp(&b),
        }
    }
}

impl PartialOrd for Value {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Value {
    /// Attempts to create a value from a [`Number`].
    ///
    /// * If the [`Number`] can be turned into a [`Value`], a numerical position
    ///   wrapped in [`Some`] will be returned.
    /// * Else, [`None`] will be returned.
    ///
    /// # Notes
    ///
    /// See the note on [`Value`] regarding when this method should be used.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::position::value::Number;
    /// use omics_coordinate::system::Zero;
    ///
    /// let value = Value::try_new(0).unwrap();
    /// assert_eq!(value.get(), Some(0));
    ///
    /// let value = Value::try_new(Number::MAX);
    /// assert_eq!(value, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub(crate) fn try_new(value: Number) -> Option<Self> {
        match value {
            LOWER_BOUND_VALUE => None,
            v => Some(Self(v)),
        }
    }

    /// Creates a lower-bound value.
    ///
    /// See [`Kind::LowerBound`] for more information.
    ///
    /// # Notes
    ///
    /// See the note on [`Value`] regarding when this method should be used.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::position::Value;
    /// use omics_coordinate::system::Zero;
    ///
    /// let value = Value::lower_bound();
    /// assert_eq!(coordinate.position().inner().kind(), Some(1));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub(crate) fn lower_bound() -> Self {
        Self(LOWER_BOUND_VALUE)
    }

    /// Attempts to get the value of the position as a [`Number`].
    ///
    /// * If the value is a [`Kind::Numerical`], a [`Some`] is returned with the
    ///   numerical position enclosed.
    /// * If the value is a [`Kind::LowerBound`], [`None`] is returned.
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
    pub fn get(&self) -> Option<Number> {
        match self.0 {
            LOWER_BOUND_VALUE => None,
            v => Some(v),
        }
    }

    /// Gets the kind of the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::system::Zero;
    ///
    /// let coordinate = "seq0:+:1".parse::<Coordinate<Zero>>()?;
    /// assert_eq!(coordinate.position().inner().kind(), Some(1));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn kind(&self) -> Kind {
        match self.0 {
            LOWER_BOUND_VALUE => Kind::LowerBound,
            _ => Kind::Numerical,
        }
    }
}

impl std::str::FromStr for Value {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == LOWER_BOUND_CHAR {
            return Ok(Self::lower_bound());
        }

        let value = s
            .parse::<Number>()
            .map_err(|_| Error::ParseError(ParseError::InvalidValue(s.to_owned())))?;

        Ok(Value(value))
    }
}

impl std::fmt::Display for Value {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.kind() {
            Kind::LowerBound => write!(f, "{LOWER_BOUND_CHAR}"),
            Kind::Numerical => write!(f, "{}", self.0),
        }
    }
}

// NOTE: this is added for convenience—most values that are manually created
// will be lower than [`u16::MAX`] and, since that will always fit into the
// allowable range of a value, we can allow creation of a value from one without
// checking.
/// A macro to implement [`From`] for smaller integers.
macro_rules! from_small_int {
    ($from:ty) => {
        impl From<$from> for Value {
            fn from(value: $from) -> Self {
                Self(value as Number)
            }
        }
    };
}

#[cfg(feature = "position-u64")]
from_small_int!(u32);
from_small_int!(u16);
from_small_int!(u8);

impl TryFrom<Value> for Number {
    type Error = Error;

    fn try_from(value: Value) -> Result<Self, Self::Error> {
        match value.0 {
            LOWER_BOUND_VALUE => Err(Error::LowerBoundToUsize),
            v => Ok(v),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn memory_layout() {
        // A kind should fit in a single byte on all platforms.
        assert_eq!(std::mem::size_of::<Kind>(), 1);

        // A position should fit in a [`Number`] on all platforms.
        assert_eq!(std::mem::size_of::<Value>(), std::mem::size_of::<Number>());
    }

    #[test]
    fn deserialization() {
        let value = "1".parse::<Value>()?;
        assert!(matches!(value.get(), Some(1)));
    }

    #[test]
    fn from_and_try_into_number() {
        let value = Value::try_new(1).unwrap();
        assert_eq!(TryInto::<Number>::try_into(value)?, 1);
    }

    #[test]
    fn a_lower_bound_cannot_convert_to_number() {
        let value = Value::lower_bound();
        let err = TryInto::<Number>::try_into(value).unwrap_err();
        assert!(matches!(err, Error::LowerBoundToUsize));
    }

    #[test]
    fn values_are_ordered_correctly() {
        assert!(Value::lower_bound() == Value::lower_bound());
        assert!(Value::lower_bound() < Value::try_new(0).unwrap());
        assert!(Value::try_new(0).unwrap() > Value::lower_bound());
        assert!(Value::try_new(0).unwrap() < Value::try_new(1).unwrap());
        assert!(Value::try_new(1).unwrap() == Value::try_new(1).unwrap());
    }

    /// Since the lower-bound is a special value, we need to make sure a user
    /// can't accidentally create a lower-bound using [`Value::try_new`].
    #[test]
    fn lower_bound_cannot_be_created_using_try_new() {
        assert_eq!(
            Value::try_new(Number::MAX - 1).unwrap().get(),
            Some(Number::MAX - 1)
        );
        assert_eq!(Value::try_new(Number::MAX), None);
    }

    #[test]
    fn values_can_be_created_from_u32() {
        #[cfg(feature = "position-u64")]
        let v = Value::from(u32::MAX);
        #[cfg(feature = "position-u64")]
        assert_eq!(v.get(), Some(u32::MAX as Number));

        let v = Value::from(u16::MAX);
        assert_eq!(v.get(), Some(u16::MAX as Number));

        let v = Value::from(u8::MAX);
        assert_eq!(v.get(), Some(u8::MAX as Number));
    }

    #[test]
    fn value_kinds_are_ordered_correctly() {
        assert!(Kind::LowerBound < Kind::Numerical);
    }

    #[test]
    fn value_always_gets_when_the_kind_is_numerical() {
        for n in 0..Number::MAX {
            // SAFETY: since the range does not include [`Number::MAX`], this
            // should always unwrap.
            let value = Value::try_new(n).unwrap();

            match value.get() {
                Some(_) => assert!(
                    value.kind() == Kind::Numerical,
                    "any `Value` that `get()`s must be numerical"
                ),
                None => assert!(
                    value.kind() == Kind::LowerBound,
                    "the only `Value` that doesn't `get()` must be the lower bound"
                ),
            }
        }
    }
}
