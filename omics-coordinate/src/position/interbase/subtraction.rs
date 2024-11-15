//! Subtraction for interbase positions.

use crate::math::CheckedSub;
use crate::position::Number;
use crate::position::interbase::Position;

impl CheckedSub<Number> for Position {
    type Output = Self;

    fn checked_sub(&self, rhs: Number) -> Option<Self::Output> {
        self.get().checked_sub(rhs).map(Position::new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addition() {
        let ten = Position::new(10);

        let nine = ten.checked_sub(1).unwrap();
        assert_eq!(nine.get(), 9);

        let zero = nine.checked_sub(9).unwrap();
        assert_eq!(zero.get(), 0);
    }

    #[test]
    fn overflow() {
        let max = Position::new(Number::MAX);

        let zero = max.checked_sub(Number::MAX).unwrap();
        assert_eq!(zero.get(), 0);

        assert!(zero.checked_sub(1).is_none());
    }
}
