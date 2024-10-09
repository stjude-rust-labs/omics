//! Subtraction for 0-based positions.

use crate::CheckedSub;
use crate::position::Value;
use crate::position::zero::Position;

impl CheckedSub<usize> for Position {
    type Output = Self;

    fn checked_sub(&self, rhs: usize) -> Option<Self::Output> {
        checked_sub(self.inner(), rhs)
    }
}

/// Checked subtraction for 0-based positions.
fn checked_sub(lhs: &Value, rhs: usize) -> Option<Position> {
    match (lhs, rhs) {
        (Value::Usize(a), b) => match a {
            &usize::MAX => a.checked_sub(b).map(Value::Usize),
            a => {
                // SAFETY: we just checked to ensure that `lhs` was not the
                // maximum [`usize`], so this must succeed.
                let lhs_plus_one = a + 1;

                match lhs_plus_one.checked_sub(b) {
                    Some(0) => Some(Value::LowerBound),
                    lhs_plus_one => lhs_plus_one.map(|lhs_plus_one| Value::Usize(lhs_plus_one - 1)),
                }
            }
        },

        // Subtracting zero from the lower bound gives you the lower bound
        // again. Any other number subtracted from the lower bound would
        // produce an overflow (in the negative direction).
        (Value::LowerBound, rhs) => match rhs {
            0 => return Some(Position::lower_bound()),
            _ => None,
        },
    }
    .map(Position::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_subtracts_a_usize_position_from_a_usize_position_correctly()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // Standard subtraction.
        let position = Position::from(10);
        let result = position.checked_sub(5).unwrap();
        assert_eq!(result.inner(), &Value::Usize(5));
        assert_eq!(result.inner().get(), Some(5));

        // Lowest value possible.
        let position = Position::from(10);
        let result = position.checked_sub(10).unwrap();
        assert_eq!(result.inner(), &Value::Usize(0));
        assert_eq!(result.inner().get(), Some(0));

        // Lands on lower bound.
        let position = Position::from(10);
        let result = position.checked_sub(11).unwrap();
        assert_eq!(result.inner(), &Value::LowerBound);
        assert_eq!(result.inner().get(), None);

        // Overflow.
        let position = Position::from(10);
        let result = position.checked_sub(12);
        assert_eq!(result, None);

        Ok(())
    }

    #[test]
    fn it_subtracts_a_usize_position_from_a_lower_bound_position_correctly()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let position = Position::lower_bound();
        let result = position.checked_sub(0).unwrap();
        assert_eq!(result.inner(), &Value::LowerBound);
        assert_eq!(result.inner().get(), None);

        let a = Position::lower_bound();
        let result = a.checked_sub(1);
        assert_eq!(result, None);

        Ok(())
    }
}
