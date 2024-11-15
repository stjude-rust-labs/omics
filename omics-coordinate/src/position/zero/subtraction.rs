//! Subtraction for 0-based positions.

use crate::CheckedSub;
use crate::position::Value;
use crate::position::value::Kind;
use crate::position::value::Number;
use crate::position::zero::Position;

impl CheckedSub<Number> for Position {
    type Output = Self;

    fn checked_sub(&self, rhs: Number) -> Option<Self::Output> {
        checked_sub(self.inner(), rhs)
    }
}

/// Checked subtraction for 0-based positions.
fn checked_sub(lhs: &Value, rhs: Number) -> Option<Position> {
    match lhs.kind() {
        Kind::Numerical => {
            // SAFETY: if `kind()` is a [`Kind::Numerical`], this will always unwrap.
            let lhs = lhs.get().unwrap();

            match lhs {
                Number::MAX => lhs.checked_sub(rhs).and_then(Value::try_new),
                lhs => {
                    // SAFETY: we just checked to ensure that `lhs` was not the
                    // maximum [`Number`], so this must succeed.
                    let lhs_plus_one = lhs + 1;

                    match lhs_plus_one.checked_sub(rhs) {
                        Some(0) => Some(Value::lower_bound()),
                        lhs_plus_one_minus_rhs => lhs_plus_one_minus_rhs
                            .and_then(|lhs_plus_one| Value::try_new(lhs_plus_one - 1)),
                    }
                }
            }
        }
        // Subtracting zero from the lower bound gives you the lower bound
        // again. Any other number subtracted from the lower bound would
        // produce an overflow (in the negative direction).
        Kind::LowerBound => match rhs {
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
        assert_eq!(result.inner(), &Value::try_new(5).unwrap());
        assert_eq!(result.inner().get(), Some(5));

        // Lowest value possible.
        let position = Position::from(10);
        let result = position.checked_sub(10).unwrap();
        assert_eq!(result.inner(), &Value::try_new(0).unwrap());
        assert_eq!(result.inner().get(), Some(0));

        // Lands on lower bound.
        let position = Position::from(10);
        let result = position.checked_sub(11).unwrap();
        assert_eq!(result.inner(), &Value::lower_bound());
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
        assert_eq!(result.inner(), &Value::lower_bound());
        assert_eq!(result.inner().get(), None);

        let a = Position::lower_bound();
        let result = a.checked_sub(1);
        assert_eq!(result, None);

        Ok(())
    }
}
