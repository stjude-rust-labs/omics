//! Addition for 1-based positions.

use crate::CheckedAdd;
use crate::position;
use crate::position::Value;
use crate::position::one::Position;

impl CheckedAdd<usize> for Position {
    type Output = Self;

    fn checked_add(&self, rhs: usize) -> Option<Self::Output> {
        checked_add(self.inner(), rhs)
    }
}

/// Checked addition for 1-based positions.
fn checked_add(lhs: &Value, rhs: usize) -> Option<Position> {
    let result = match (lhs, rhs) {
        // Adding two one-based usize positions can simply be computed using the
        // built-in `checked_add` method for [`usize`].
        (Value::Usize(lhs), rhs) => lhs.checked_add(rhs).map(Value::Usize),

        // Adding the a one-based usize position to the lower bound will end
        // up as `0 + position`, which is simply `position`.
        (Value::LowerBound, rhs) => Some(Value::Usize(rhs)),
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
    fn it_adds_to_a_usize_position_correctly() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        let position = Position::try_from(1)?;
        let result = position.checked_add(2).unwrap();
        assert_eq!(result.inner(), &Value::Usize(3));
        assert_eq!(result.get(), Some(3));

        let position = Position::try_from(usize::MAX)?;
        let result = position.checked_add(1);
        assert_eq!(result, None);

        Ok(())
    }
}
