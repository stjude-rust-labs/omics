//! Addition for interbase positions.

use crate::math::CheckedAdd;
use crate::position::Number;
use crate::position::interbase::Position;

impl CheckedAdd<Number> for Position {
    type Output = Self;

    fn checked_add(&self, rhs: Number) -> Option<Self> {
        self.get().checked_add(rhs).map(Position::new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addition() {
        let zero = Position::new(0);

        let one = zero.checked_add(1).unwrap();
        assert_eq!(one.get(), 1);

        let two = one.checked_add(1).unwrap();
        assert_eq!(two.get(), 2);
    }

    #[test]
    fn overflow() {
        let max = Position::new(Number::MAX);

        let max = max.checked_add(0).unwrap();
        assert_eq!(max.get(), Number::MAX);

        assert!(max.checked_add(1).is_none());
    }
}
