//! Interbase positions.

use std::num::NonZero;

use crate::position::Error;
use crate::position::Number;
use crate::position::ParseError;
use crate::position::Result;
use crate::system::Interbase;

mod addition;
mod subtraction;

////////////////////////////////////////////////////////////////////////////////////////
// Assertions
////////////////////////////////////////////////////////////////////////////////////////

const _: () = {
    // Ensure that a value is only ever as big as the internal, numerical
    // representation.
    assert!(size_of::<Position>() == size_of::<Number>());

    /// A function to ensure that types are `Copy`.
    const fn is_copy<T: Copy>() {}
    is_copy::<Position>();
};

////////////////////////////////////////////////////////////////////////////////////////
// Position
////////////////////////////////////////////////////////////////////////////////////////

/// An interbase position.
///
/// Interbase positions start at zero (`0`).
pub type Position = crate::Position<Interbase>;

impl Position {
    /// Creates a new interbase position.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::position::interbase::Position;
    ///
    /// let position = Position::new(0);
    /// assert_eq!(position.get(), 0);
    /// ```
    pub const fn new(value: Number) -> Self {
        Self {
            system: Interbase,
            value,
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl super::r#trait::Position<Interbase> for Position {}

impl std::str::FromStr for Position {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Ok(Self::new(s.parse::<Number>().map_err(|error| {
            Error::Parse(ParseError::Int {
                inner: error,
                value: s.to_string(),
            })
        })?))
    }
}

impl TryFrom<Number> for Position {
    type Error = Error;

    fn try_from(value: Number) -> std::result::Result<Self, Self::Error> {
        Ok(Self::new(value))
    }
}

impl From<NonZero<Number>> for Position {
    fn from(value: NonZero<Number>) -> Self {
        Self::try_from(value.get()).unwrap()
    }
}

/// Creates implementations to convert from smaller numbers to a position.
macro_rules! position_from_smaller_number {
    ($from:ty) => {
        impl From<$from> for Position {
            fn from(value: $from) -> Self {
                Self::new(value as Number)
            }
        }

        impl From<NonZero<$from>> for Position {
            fn from(value: NonZero<$from>) -> Self {
                Self::new(value.get() as Number)
            }
        }
    };
}

#[cfg(feature = "position-u64")]
position_from_smaller_number!(u32);
position_from_smaller_number!(u16);
position_from_smaller_number!(u8);

#[cfg(test)]
mod tests {
    use std::num::NonZeroU8;
    use std::num::NonZeroU16;
    #[cfg(feature = "position-u64")]
    use std::num::NonZeroU32;

    use crate::Position;
    use crate::position::Number;
    use crate::position::Result;
    use crate::system::Interbase;

    #[test]
    fn try_from_number() {
        let position: Position<Interbase> = 1u32.try_into().unwrap();
        assert_eq!(position.get(), 1);
    }

    #[test]
    fn distance() {
        // Zero distance between two positions.
        let a = Position::<Interbase>::new(10);
        let b = Position::<Interbase>::new(10);
        assert_eq!(a.distance_unchecked(&b), 0);
        assert_eq!(b.distance_unchecked(&a), 0);

        // Non-zero distance between two Number positions.
        let a = Position::<Interbase>::new(Number::MAX);
        let b = Position::<Interbase>::new(0);
        assert_eq!(a.distance_unchecked(&b), Number::MAX);
        assert_eq!(b.distance_unchecked(&a), Number::MAX);
    }

    #[test]
    fn parse() {
        let zero = "0".parse::<Position<Interbase>>().unwrap();
        assert_eq!(zero.get(), 0);

        let position = "1".parse::<Position<Interbase>>().unwrap();
        assert_eq!(position.get(), 1);

        let err = "a".parse::<Position<Interbase>>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "parse error: invalid digit found in string: `a`"
        );
    }

    #[test]
    fn from_smaller_types() -> Result<()> {
        #[cfg(feature = "position-u64")]
        {
            // u32
            let position = Position::<Interbase>::from(0u32);
            assert_eq!(position.get(), 0);

            let position = Position::<Interbase>::from(NonZeroU32::new(1).unwrap());
            assert_eq!(position.get(), 1);
        }

        // u16
        let position = Position::<Interbase>::from(0u16);
        assert_eq!(position.get(), 0);

        let position = Position::<Interbase>::from(NonZeroU16::new(1).unwrap());
        assert_eq!(position.get(), 1);

        // u8
        let position = Position::<Interbase>::from(0u8);
        assert_eq!(position.get(), 0);

        let position = Position::<Interbase>::from(NonZeroU8::new(1).unwrap());
        assert_eq!(position.get(), 1);

        Ok(())
    }
}
