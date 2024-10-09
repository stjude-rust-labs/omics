//! A 0-based, half-open coordinate system.

use crate::system::System;

/// 0-based, half-open coordinate system.
#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct Zero;

impl System for Zero {}

impl std::fmt::Display for Zero {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "0-based, half-open coordinate system")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_remains_zero_sized() {
        assert_eq!(std::mem::size_of::<Zero>(), 0);
    }
}
