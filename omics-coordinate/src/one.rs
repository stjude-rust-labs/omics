//! Convenient access to a 1-based [`Coordinate`](crate::Coordinate).

use crate::system::One;

/// A 1-based [`Coordinate`](crate::Coordinate).
pub type Coordinate = crate::Coordinate<One>;

impl crate::r#trait::Coordinate<One> for Coordinate {}
