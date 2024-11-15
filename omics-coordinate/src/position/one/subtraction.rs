//! Subtraction for 1-based positions.

use crate::CheckedSub;
use crate::position::Value;
use crate::position::one::Position;
use crate::position::value::Kind;
use crate::position::value::Number;

impl CheckedSub<Number> for Position {
    type Output = Self;

    fn checked_sub(&self, rhs: Number) -> Option<Self::Output> {
        checked_sub(self.inner(), rhs)
    }
}

/// Checked subtraction for 1-based positions.
fn checked_sub(lhs: &Value, rhs: Number) -> Option<Position> {
    match (lhs.kind(), rhs) {
        // Subtracting two one-based number positions can simply be computed using the
        // built-in `checked_sub` method for [`Number`].
        (Kind::Numerical, rhs) => {
            // SAFETY: if `kind()` is a [`Kind::Numerical`], this will always unwrap.
            let lhs = lhs.get().unwrap();
            lhs.checked_sub(rhs).and_then(Value::try_new)
        }

        // Subtracting the a one-based Number position from the lower bound position
        // will end up as `0 - position`. However, not that value cannot be
        // zero because that is not allowed in a one-based position. Thus, this
        // would always overflow in the negative direction.
        (Kind::LowerBound, _) => None,
    }
    .and_then(|v| Position::try_new(v).ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_subtracts_from_a_usize_position_correctly()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // Standard subtraction.
        let position = Position::try_from(10)?;
        let result = position.checked_sub(5).unwrap();
        assert_eq!(result.inner(), &Value::try_from(5).unwrap());
        assert_eq!(result.get(), Some(5));

        // Lowest value possible.
        let position = Position::try_from(10)?;
        let result = position.checked_sub(9).unwrap();
        assert_eq!(result.inner(), &Value::try_from(1).unwrap());
        assert_eq!(result.get(), Some(1));

        // Overflow.
        let position = Position::try_from(10)?;
        let result = position.checked_sub(10);
        assert_eq!(result, None);

        // Overflow.
        let position = Position::try_from(10)?;
        let result = position.checked_sub(11);
        assert_eq!(result, None);

        Ok(())
    }
}
