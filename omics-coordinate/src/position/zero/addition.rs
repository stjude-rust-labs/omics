//! Addition for 0-based positions.

use crate::CheckedAdd;
use crate::position::Value;
use crate::position::value::Kind;
use crate::position::value::Number;
use crate::position::zero::Position;

impl CheckedAdd<Number> for Position {
    type Output = Self;

    fn checked_add(&self, rhs: Number) -> Option<Self> {
        checked_add(self.inner(), rhs)
    }
}

/// Checked addition for 0-based positions.
fn checked_add(lhs: &Value, rhs: Number) -> Option<Position> {
    match (lhs.kind(), rhs) {
        // Adding two zero-based number positions can simply be added using the
        // built-in `checked_add` method for [`Number`].
        (Kind::Numerical, rhs) => {
            // SAFETY: if `kind()` is a [`Kind::Numerical`], this will always unwrap.
            let lhs = lhs.get().unwrap();
            lhs.checked_add(rhs).and_then(Value::try_new)
        }

        // Adding the a zero-based Number position to the lower bound will end up
        // as `-1 + position`, which is `position - 1`.  However, if the
        // zero-based number position is zero, then the result would be the
        // identity of the `lhs` operand, which is the lower bound. Thus, that
        // case needs to be handled separately.
        (Kind::LowerBound, rhs) => match rhs {
            0 => return Some(Position::lower_bound()),
            rhs => rhs.checked_sub(1).and_then(Value::try_new),
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
        assert_eq!(result.inner(), &Value::try_new(3).unwrap());
        assert_eq!(result.inner().get(), Some(3));

        let position = Position::from(Number::MAX);
        let result = position.checked_add(1);
        assert_eq!(result, None);

        Ok(())
    }

    #[test]
    fn it_adds_a_lower_bound_position_and_a_usize_position_together_correctly()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let position = Position::lower_bound();
        let result = position.checked_add(0).unwrap();
        assert_eq!(result.inner(), &Value::lower_bound());
        assert_eq!(result.inner().get(), None);

        let position = Position::lower_bound();
        let result = position.checked_add(1).unwrap();
        assert_eq!(result.inner(), &Value::try_new(0).unwrap());
        assert_eq!(result.inner().get(), Some(0));

        let position = Position::lower_bound();
        let result = position.checked_add(Number::MAX).unwrap();
        assert_eq!(result.inner(), &Value::try_new(Number::MAX - 1).unwrap());
        assert_eq!(result.inner().get(), Some(Number::MAX - 1));

        Ok(())
    }
}
