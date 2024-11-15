//! Addition for 1-based positions.

use crate::CheckedAdd;
use crate::position::Value;
use crate::position::one::Position;
use crate::position::value::Kind;
use crate::position::value::Number;

impl CheckedAdd<Number> for Position {
    type Output = Self;

    fn checked_add(&self, rhs: Number) -> Option<Self::Output> {
        checked_add(self.inner(), rhs)
    }
}

/// Checked addition for 1-based positions.
fn checked_add(lhs: &Value, rhs: Number) -> Option<Position> {
    match (lhs.kind(), rhs) {
        // Adding two one-based number positions can simply be computed using the
        // built-in `checked_add` method for [`Number`].
        (Kind::Numerical, rhs) => {
            // SAFETY: if `kind()` is a [`Kind::Numerical`], this will always unwrap.
            let lhs = lhs.get().unwrap();
            lhs.checked_add(rhs).and_then(Value::try_new)
        }

        // Adding the a one-based Number position to the lower bound will end
        // up as `0 + position`, which is simply `position`.
        (Kind::LowerBound, rhs) => Value::try_new(rhs),
    }
    .and_then(|v| Position::try_new(v).ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_adds_to_a_usize_position_correctly() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        let position = Position::try_from(1)?;
        let result = position.checked_add(2).unwrap();
        assert_eq!(result.inner(), &Value::try_new(3).unwrap());
        assert_eq!(result.get(), Some(3));

        let position = Position::try_from(Number::MAX)?;
        let result = position.checked_add(1);
        assert_eq!(result, None);

        Ok(())
    }
}
