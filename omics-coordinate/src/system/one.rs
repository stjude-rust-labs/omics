//! A 1-based, fully-closed coordinate system.

use crate::system::System;

/// 1-based, fully-closed coordinate system.
#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct One;

impl System for One {}

impl std::fmt::Display for One {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "1-based, fully-closed coordinate system")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_remains_zero_sized() {
        assert_eq!(std::mem::size_of::<One>(), 0);
    }
}
