//! The base coordinate system.

use crate::system::System;

////////////////////////////////////////////////////////////////////////////////////////
// Assertions
////////////////////////////////////////////////////////////////////////////////////////

const _: () = {
    use std::mem::size_of;

    // This should never take up any space.
    assert!(size_of::<Base>() == 0);
};

////////////////////////////////////////////////////////////////////////////////////////
// Base coordinate system
////////////////////////////////////////////////////////////////////////////////////////

/// The base coordinate system.
///
/// This coordinate system is also known as the "1-based, fully-closed"
/// coordinate system.
#[derive(Copy, Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct Base;

impl Base {
    /// A static string representing the name of the base coordinate system.
    pub const NAME: &str = "base coordinate system";
}

impl System for Base {}

impl std::fmt::Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Self::NAME)
    }
}
