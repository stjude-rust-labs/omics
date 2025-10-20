//! Base positions.

mod addition;
mod subtraction;

use std::num::NonZero;

use crate::position::Error;
use crate::position::Number;
use crate::position::ParseError;
use crate::position::Result;
use crate::system::Base;

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

/// A base position.
///
/// Base positions start at one (`1`).
pub type Position = crate::Position<Base>;

impl Position {
    /// Create a new base position if the value is not zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::position::base::Position;
    ///
    /// let position = Position::try_new(1)?;
    /// assert_eq!(position.get(), 1);
    ///
    /// # Ok::<(), omics_coordinate::position::Error>(())
    /// ```
    pub const fn try_new(value: Number) -> Result<Self> {
        if value == 0 {
            return Err(Error::IncompatibleValue {
                system: Base::NAME,
                value,
            });
        }

        Ok(Self {
            system: Base,
            value,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl super::r#trait::Position<Base> for Position {}

impl std::str::FromStr for Position {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_new(s.parse::<Number>().map_err(|error| ParseError::Int {
            system: Base::NAME,
            inner: error,
            value: s.to_string(),
        })?)
    }
}

impl TryFrom<Number> for Position {
    type Error = Error;

    fn try_from(value: Number) -> Result<Self> {
        Self::try_new(value)
    }
}

impl From<NonZero<Number>> for Position {
    fn from(value: NonZero<Number>) -> Self {
        // SAFETY: because [`try_new()`] will only throw an error when zero
        // (`0`) is passed in and `value` here is a non-zero number, this will
        // always [`unwrap()`].
        Self::try_new(value.get()).unwrap()
    }
}

/// Creates implementations to convert from smaller numbers to a position.
macro_rules! position_from_smaller_number {
    ($from:ty) => {
        impl From<NonZero<$from>> for Position {
            fn from(value: NonZero<$from>) -> Self {
                // SAFETY: because [`try_from()`] will only throw an error when zero
                // (`0`) is passed in and `value` here is a non-zero number, this will
                // always [`unwrap()`].
                Self::try_new(value.get() as Number).unwrap()
            }
        }

        impl TryFrom<$from> for Position {
            type Error = Error;

            fn try_from(value: $from) -> Result<Self> {
                Self::try_new(value as Number)
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
    use crate::position::Error;
    use crate::position::Number;
    use crate::position::Result;
    use crate::system::Base;

    #[test]
    fn from_number() {
        let error: Result<Position<Base>> = 0u32.try_into();
        assert_eq!(
            error.unwrap_err(),
            Error::IncompatibleValue {
                system: Base::NAME,
                value: 0,
            }
        );

        let position: Position<Base> = 1u32.try_into().unwrap();
        assert_eq!(position.get(), 1);
    }

    #[test]
    fn distance() {
        // Zero distance between two positions.
        let a = Position::<Base>::try_from(10u8).unwrap();
        let b = Position::<Base>::try_from(10u8).unwrap();
        assert_eq!(a.distance_unchecked(&b), 0);
        assert_eq!(b.distance_unchecked(&a), 0);

        // Non-zero distance between two Number positions.
        let a = Position::<Base>::try_from(Number::MAX).unwrap();
        let b = Position::<Base>::try_from(1u8).unwrap();
        assert_eq!(a.distance_unchecked(&b), Number::MAX - 1);
        assert_eq!(b.distance_unchecked(&a), Number::MAX - 1);
    }

    #[test]
    fn parse() {
        let err = "0".parse::<Position<Base>>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "incompatible value for system \"base coordinate system\": `0`"
        );

        let position = "1".parse::<Position<Base>>().unwrap();
        assert_eq!(position.get(), 1);

        let err = "a".parse::<Position<Base>>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "parse error: failed to parse base coordinate system position from `a`: invalid digit \
             found in string"
        );
    }

    #[test]
    fn from_smaller_types() -> Result<()> {
        #[cfg(feature = "position-u64")]
        {
            // u32
            let position = Position::<Base>::try_from(1u32)?;
            assert_eq!(position.get(), 1);

            let position = Position::<Base>::from(NonZeroU32::new(1).unwrap());
            assert_eq!(position.get(), 1);
        }

        // u16
        let position = Position::<Base>::try_from(1u16)?;
        assert_eq!(position.get(), 1);

        let position = Position::<Base>::from(NonZeroU16::new(1).unwrap());
        assert_eq!(position.get(), 1);

        // u8
        let position = Position::<Base>::try_from(1u8)?;
        assert_eq!(position.get(), 1);

        let position = Position::<Base>::from(NonZeroU8::new(1).unwrap());
        assert_eq!(position.get(), 1);

        Ok(())
    }
}
