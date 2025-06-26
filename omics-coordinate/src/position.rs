//! Positions.

use std::num::ParseIntError;

use thiserror::Error;

use crate::System;
use crate::math::CheckedAdd;
use crate::math::CheckedSub;
use crate::system::Base;
use crate::system::Interbase;

pub mod base;
pub mod interbase;

////////////////////////////////////////////////////////////////////////////////////////
// Constants and Types
////////////////////////////////////////////////////////////////////////////////////////

/// The inner representation for a numerical position value.
///
/// Note that `u64` positions can be enabled by turning on the `position-u64`
/// feature for the crate.
#[cfg(not(feature = "position-u64"))]
pub type Number = u32;

/// The inner representation for a numerical position value.
///
/// Note that `u32` positions can be enabled by turning off the `position-u64`
/// feature for the crate.
#[cfg(feature = "position-u64")]
pub type Number = u64;

////////////////////////////////////////////////////////////////////////////////////////
// Assertions
////////////////////////////////////////////////////////////////////////////////////////

const _: () = {
    /// A function to ensure that types are `Copy`.
    const fn is_copy<T: Copy>() {}
    is_copy::<Number>();

    // Ensure that the types themselves are copy, as they should be able to be
    // passed around as such as well.
    is_copy::<Position<Interbase>>();
    is_copy::<Position<Base>>();
};

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// A error related to parsing a position.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// An integer parsing error.
    ///
    /// Occurs when an integer position value cannot be parsed.
    #[error("{inner}: `{value}`")]
    Int {
        /// The inner error.
        inner: ParseIntError,

        /// The value that was attempted to be parsed.
        value: String,
    },
}

/// A [`Result`](std::result::Result) with a [`ParseError`].
pub type ParseResult<T> = std::result::Result<T, ParseError>;

/// A position-related error.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// A parse error.
    #[error("parse error: {0}")]
    Parse(#[from] ParseError),

    /// Incompatible value.
    ///
    /// This error represents and incompatible value for a position that is
    /// placed within a coordinate system. For example, zero (`0`) is not a
    /// valid numerical position within a 1-based coordinate system.
    #[error("incompatible value for system \"{system}\": `{value}`")]
    IncompatibleValue {
        /// The system within when the value is incompatible.
        system: &'static str,

        /// The incompatible value.
        value: Number,
    },
}

/// A [`Result`](std::result::Result) with an [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

///////////////////////////////////////////////////////////////////////////////////////
// The `Position` trait
///////////////////////////////////////////////////////////////////////////////////////

/// Traits related to a position.
pub mod r#trait {
    use std::num::NonZero;

    use super::*;

    /// Requirements to be a position.
    pub trait Position<S: System>:
        std::fmt::Display
        + std::fmt::Debug
        + PartialEq
        + Eq
        + PartialOrd
        + Ord
        + std::str::FromStr<Err = Error>
        + CheckedAdd<Number, Output = Self>
        + CheckedSub<Number, Output = Self>
        + TryFrom<Number>
        + From<NonZero<Number>>
    where
        Self: Sized,
    {
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Positions
////////////////////////////////////////////////////////////////////////////////////////

/// An offset from the start of a molecule.
///
/// For a more in-depth discussion on what positions are and the notations used
/// within this crate, please see [this section of the docs](crate#positions).
#[derive(Copy, Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Position<S: System> {
    /// The coordinate system.
    system: S,

    /// The inner value.
    value: Number,
}

impl<S: System> Position<S> {
    /// Gets the numerical position.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let position = Position::<Interbase>::new(42);
    /// assert_eq!(position.get(), 42);
    /// ```
    pub fn get(&self) -> Number {
        self.value
    }

    /// Performs checked addition.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let position = Position::<Interbase>::new(42)
    ///     .checked_add(8)
    ///     .expect("addition to succeed");
    /// assert_eq!(position.get(), 50);
    /// ```
    pub fn checked_add(&self, rhs: Number) -> Option<Self>
    where
        Self: r#trait::Position<S>,
    {
        <Self as CheckedAdd<Number>>::checked_add(self, rhs)
    }

    /// Performs checked subtraction.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use omics_coordinate::Position;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let position = Position::<Interbase>::new(42)
    ///     .checked_sub(2)
    ///     .expect("subtraction to succeed");
    /// assert_eq!(position.get(), 40);
    /// ```
    pub fn checked_sub(&self, rhs: Number) -> Option<Self>
    where
        Self: r#trait::Position<S>,
    {
        <Self as CheckedSub<Number>>::checked_sub(self, rhs)
    }

    /// Gets the magnitude of the distance between two positions.
    ///
    /// # Note
    ///
    /// This method calculates the magnitude of distance between two positions
    /// that are assumed to be on the same number line (i.e., the same stand and
    /// contig). Notably, **there is not check** regarding strand or contig
    /// equivalent within this method.
    ///
    /// Because the use case of doing this is so niche, the method is currently
    /// only accessible within the crate.  If you're wanting to do this kind of
    /// thing, you're probably going to want to convert the position to
    /// coordinates and calculate the distance between the two coordinates. If
    /// you think you have a legitimate use case where this would be useful,
    /// please file an issue.
    pub(crate) fn distance_unchecked(&self, rhs: &Position<S>) -> Number {
        let a = self.get();
        let b = rhs.get();

        // SAFETY: because these are two unsigned numbers that are being
        // subtracted correctly (based on the `if` statement below, we
        // expect these to always unwrap).
        if a >= b {
            a.checked_sub(b).unwrap()
        } else {
            b.checked_sub(a).unwrap()
        }
    }
}

impl<S: System> std::fmt::Display for Position<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if !f.alternate() {
            write!(f, "{}", self.value)
        } else {
            write!(f, "{} ({})", self.value, self.system)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;

    use super::*;
    use crate::system::Interbase;

    #[test]
    fn serialize() {
        let position = Position::<Interbase>::from(0u8);

        let mut buffer = String::new();
        write!(&mut buffer, "{position}").unwrap();
        assert_eq!(buffer, "0");

        buffer.clear();
        write!(&mut buffer, "{position:#}").unwrap();
        assert_eq!(buffer, "0 (interbase coordinate system)");
    }
}
