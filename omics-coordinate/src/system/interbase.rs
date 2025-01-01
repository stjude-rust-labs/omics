//! The interbase coordinate system.

use crate::system::r#trait::System;

////////////////////////////////////////////////////////////////////////////////////////
// Assertions
////////////////////////////////////////////////////////////////////////////////////////

const _: () = {
    use std::mem::size_of;

    // This should never take up any space.
    assert!(size_of::<Interbase>() == 0);
};

////////////////////////////////////////////////////////////////////////////////////////
// Interbase coordinate system
////////////////////////////////////////////////////////////////////////////////////////

/// The interbase coordinate system.
///
/// This coordinate system is also known as the "0-based, half-open" coordinate
/// system.
#[derive(Copy, Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct Interbase;

impl Interbase {
    /// A static string representing the name of the interbase coordinate
    /// system.
    pub const NAME: &str = "interbase coordinate system";
}

impl System for Interbase {}

impl std::fmt::Display for Interbase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Self::NAME)
    }
}
