//! Interbase spans.

use crate::Position;
use crate::position::Number;
use crate::span::Direction;
use crate::system::Base;
use crate::system::Interbase;

/// An interbase span.
pub type Span = crate::Span<Interbase>;

impl super::r#trait::Span<Interbase> for Span {
    fn contains_entity(&self, position: &Position<Base>) -> bool {
        match self.direction() {
            Direction::Ascending => {
                self.start().get() < position.get() && position.get() <= self.end().get()
            }
            Direction::Stationary => false,
            Direction::Descending => {
                self.end().get() < position.get() && position.get() <= self.start().get()
            }
        }
    }

    fn count_entities(&self) -> Number {
        self.start().distance_unchecked(self.end())
    }

    fn is_empty(&self) -> bool {
        self.start() == self.end()
    }
}
