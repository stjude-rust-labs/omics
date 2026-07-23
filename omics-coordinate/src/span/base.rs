//! Base spans.

use crate::Position;
use crate::position::Number;
use crate::system::Base;

/// A base span.
pub type Span = crate::Span<Base>;

impl super::r#trait::Span<Base> for Span {
    fn contains_entity(&self, position: &Position<Base>) -> bool {
        self.contains_position(position)
    }

    fn count_entities(&self) -> Number {
        self.start().distance_unchecked(self.end()) + 1
    }

    fn is_empty(&self) -> bool {
        false
    }
}
