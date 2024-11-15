//! Subtraction for base positions.

use crate::math::CheckedSub;
use crate::position::Number;
use crate::position::base::Position;

impl CheckedSub<Number> for Position {
    type Output = Self;

    fn checked_sub(&self, rhs: Number) -> Option<Self::Output> {
        self.get().checked_sub(rhs).map(Position::try_new)?.ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addition() {
        let ten = Position::try_new(10).unwrap();

        let nine = ten.checked_sub(1).unwrap();
        assert_eq!(nine.get(), 9);

        let one = nine.checked_sub(8).unwrap();
        assert_eq!(one.get(), 1);
    }

    #[test]
    fn overflow() {
        let max = Position::try_new(Number::MAX).unwrap();

        let one = max.checked_sub(Number::MAX - 1).unwrap();
        assert_eq!(one.get(), 1);

        assert!(one.checked_sub(1).is_none());
    }
}
