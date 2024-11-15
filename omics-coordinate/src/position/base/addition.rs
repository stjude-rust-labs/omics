//! Addition for base positions.

use crate::math::CheckedAdd;
use crate::position::Number;
use crate::position::base::Position;

impl CheckedAdd<Number> for Position {
    type Output = Self;

    fn checked_add(&self, rhs: Number) -> Option<Self> {
        self.get().checked_add(rhs).map(Position::try_new)?.ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addition() {
        let one = Position::try_new(1).unwrap();

        let two = one.checked_add(1).unwrap();
        assert_eq!(two.get(), 2);

        let three = two.checked_add(1).unwrap();
        assert_eq!(three.get(), 3);
    }

    #[test]
    fn overflow() {
        let max = Position::try_from(Number::MAX).unwrap();

        let max = max.checked_add(0).unwrap();
        assert_eq!(max.get(), Number::MAX);

        assert!(max.checked_add(1).is_none());
    }
}
