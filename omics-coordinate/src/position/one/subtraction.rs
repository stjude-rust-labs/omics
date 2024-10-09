//! Subtraction for 1-based positions.

use crate::CheckedSub;
use crate::position;
use crate::position::Value;
use crate::position::one::Position;

impl CheckedSub<usize> for Position {
    type Output = Self;

    fn checked_sub(&self, rhs: usize) -> Option<Self::Output> {
        checked_sub(self.inner(), rhs)
    }
}

/// Checked subtraction for 1-based positions.
fn checked_sub(lhs: &Value, rhs: usize) -> Option<Position> {
    let result = match (lhs, rhs) {
        // Subtracting two one-based usize positions can simply be computed using the
        // built-in `checked_sub` method for [`usize`].
        (Value::Usize(lhs), rhs) => lhs.checked_sub(rhs).map(Value::Usize),

        // Subtracting the a one-based usize position from the lower bound position
        // will end up as `0 - position`. However, not that value cannot be
        // zero because that is not allowed in a one-based position. Thus, this
        // would always overflow in the negative direction.
        (Value::LowerBound, _) => None,
    };

    // NOTE: at the time of writing, the only error that `Position::try_new()` is an
    // [`Error::IncompatibleValueForSystem`]. In this case, if we encounter that
    // error, we don't want to bubble up the error—instead, we just want to
    // treat the move as if it's out of bounds.
    //
    // Here, we treat this situation delicately to ensure any new errors introduced
    // into [`Position::try_new()`] would be handled correctly in the future. In
    // short, we explicitly return [`None`] when that specific error is
    // encountered. Otherwise, we unwrap.
    result
        .map(Position::try_new)
        .transpose()
        .or_else(|err| match err {
            position::Error::Parse(position::ParseError::IncompatibleValue(..)) => Ok(None),
            err => Err(err),
        })
        .unwrap()
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
        assert_eq!(result.inner(), &Value::Usize(5));
        assert_eq!(result.get(), Some(5));

        // Lowest value possible.
        let position = Position::try_from(10)?;
        let result = position.checked_sub(9).unwrap();
        assert_eq!(result.inner(), &Value::Usize(1));
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
