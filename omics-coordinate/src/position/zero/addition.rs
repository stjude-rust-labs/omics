//! Addition for 0-based positions.

use crate::CheckedAdd;
use crate::position::Value;
use crate::position::zero::Position;

impl CheckedAdd<usize> for Position {
    type Output = Self;

    fn checked_add(&self, rhs: usize) -> Option<Self> {
        checked_add(self.inner(), rhs)
    }
}

/// Checked addition for 0-based positions.
fn checked_add(lhs: &Value, rhs: usize) -> Option<Position> {
    match (lhs, rhs) {
        // Adding two zero-based usize positions can simply be added using the
        // built-in `checked_add` method for [`usize`].
        (Value::Usize(lhs), rhs) => lhs.checked_add(rhs).map(Value::Usize),

        // Adding the a zero-based usize position to the lower bound will end up
        // as `-1 + position`, which is `position - 1`.  However, if the
        // zero-based usize position is zero, then the result would be the
        // identity of the `lhs` operand, which is the lower bound. Thus, that
        // case needs to be handled separately.
        (Value::LowerBound, rhs) => match rhs {
            0 => return Some(Position::lower_bound()),
            rhs => rhs.checked_sub(1).map(Value::Usize),
        },
    }
    .map(Position::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_adds_usize_positions_together_correctly()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let position = Position::from(1);
        let result = position.checked_add(2).unwrap();
        assert_eq!(result.inner(), &Value::Usize(3));
        assert_eq!(result.inner().get(), Some(3));

        let position = Position::from(usize::MAX);
        let result = position.checked_add(1);
        assert_eq!(result, None);

        Ok(())
    }

    #[test]
    fn it_adds_a_lower_bound_position_and_a_usize_position_together_correctly()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let position = Position::lower_bound();
        let result = position.checked_add(0).unwrap();
        assert_eq!(result.inner(), &Value::LowerBound);
        assert_eq!(result.inner().get(), None);

        let position = Position::lower_bound();
        let result = position.checked_add(1).unwrap();
        assert_eq!(result.inner(), &Value::Usize(0));
        assert_eq!(result.inner().get(), Some(0));

        let position = Position::lower_bound();
        let result = position.checked_add(usize::MAX).unwrap();
        assert_eq!(result.inner(), &Value::Usize(usize::MAX - 1));
        assert_eq!(result.inner().get(), Some(usize::MAX - 1));

        Ok(())
    }
}
