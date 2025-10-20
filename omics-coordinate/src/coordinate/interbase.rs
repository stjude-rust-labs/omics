//! Convenient access to an interbase [`Coordinate`](crate::Coordinate).

use crate::Strand;
use crate::contig;
use crate::coordinate::Error;
use crate::position;
use crate::position::interbase::Position;
use crate::strand;
use crate::system::Base;
use crate::system::Interbase;

/// An interbase coordinate.
pub type Coordinate = crate::Coordinate<Interbase>;

impl Coordinate {
    /// Consumes `self` and attempts converts the coordinate to the next in-base
    /// coordinate.
    ///
    /// This is only useful when converting a coordinate to the other coordinate
    /// system, which is not a common operation. If, instead, you wish to move
    /// the coordinate forward within the _same_ coordinate system, you're
    /// almost certainly looking for
    /// [`move_forward()`](crate::Coordinate::move_forward).
    pub fn nudge_forward(self) -> Option<crate::Coordinate<Base>> {
        let (contig, strand, position) = self.into_parts();

        let position = match strand {
            Strand::Positive => position
                .get()
                .checked_add(1)
                .and_then(|value| crate::Position::<Base>::try_new(value).ok()),
            Strand::Negative => crate::Position::<Base>::try_new(position.get()).ok(),
        }?;

        Some(crate::Coordinate::new(contig, strand, position))
    }

    /// Consumes `self` and converts the coordinate to the previous in-base
    /// coordinate.
    ///
    /// This is only useful when converting a coordinate to the other coordinate
    /// system, which is not a common operation. If, instead, you wish to move
    /// the coordinate backward within the _same_ coordinate system, you're
    /// almost certainly looking for
    /// [`move_backward()`](crate::Coordinate::move_backward).
    pub fn nudge_backward(self) -> Option<crate::Coordinate<Base>> {
        let (contig, strand, position) = self.into_parts();

        let position = match strand {
            Strand::Positive => crate::Position::<Base>::try_new(position.get()).ok(),
            Strand::Negative => position
                .get()
                .checked_add(1)
                .and_then(|value| crate::Position::<Base>::try_new(value).ok()),
        }?;

        Some(crate::Coordinate::new(contig, strand, position))
    }
}

impl crate::coordinate::r#trait::Coordinate<Interbase> for Coordinate {
    fn try_new(
        contig: impl TryInto<crate::Contig, Error = contig::Error>,
        strand: impl TryInto<crate::Strand, Error = strand::Error>,
        position: position::Number,
    ) -> super::Result<Self> {
        let contig = contig.try_into().map_err(Error::Contig)?;
        let strand = strand.try_into().map_err(Error::Strand)?;
        let position = Position::new(position);

        Ok(Self {
            system: Interbase,
            contig,
            strand,
            position,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base;
    use crate::position::Number;

    fn create_coordinate(contig: &str, strand: &str, position: Number) -> Coordinate {
        Coordinate::try_new(contig, strand, position).unwrap()
    }

    fn create_base_coordinate(contig: &str, strand: &str, position: Number) -> base::Coordinate {
        base::Coordinate::try_new(contig, strand, position).unwrap()
    }

    #[test]
    fn nudge_forward() {
        let coordinate = create_coordinate("seq0", "+", 0);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_base_coordinate("seq0", "+", 1)
        );

        let coordinate = create_coordinate("seq0", "+", 1);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_base_coordinate("seq0", "+", 2)
        );

        let coordinate = create_coordinate("seq0", "+", 10);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_base_coordinate("seq0", "+", 11)
        );

        let coordinate = create_coordinate("seq0", "+", Number::MAX);
        assert!(coordinate.nudge_forward().is_none());

        let coordinate = create_coordinate("seq0", "-", 0);
        assert!(coordinate.nudge_forward().is_none());

        let coordinate = create_coordinate("seq0", "-", 1);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_base_coordinate("seq0", "-", 1)
        );

        let coordinate = create_coordinate("seq0", "-", 10);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_base_coordinate("seq0", "-", 10)
        );

        let coordinate = create_coordinate("seq0", "-", Number::MAX);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_base_coordinate("seq0", "-", Number::MAX)
        );
    }

    #[test]
    fn nudge_backward() {
        let coordinate = create_coordinate("seq0", "+", 0);
        assert!(coordinate.nudge_backward().is_none());

        let coordinate = create_coordinate("seq0", "+", 1);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_base_coordinate("seq0", "+", 1)
        );

        let coordinate = create_coordinate("seq0", "+", 10);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_base_coordinate("seq0", "+", 10)
        );

        let coordinate = create_coordinate("seq0", "+", Number::MAX);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_base_coordinate("seq0", "+", Number::MAX)
        );

        let coordinate = create_coordinate("seq0", "-", 0);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_base_coordinate("seq0", "-", 1)
        );

        let coordinate = create_coordinate("seq0", "-", 1);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_base_coordinate("seq0", "-", 2)
        );

        let coordinate = create_coordinate("seq0", "-", 10);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_base_coordinate("seq0", "-", 11)
        );

        let coordinate = create_coordinate("seq0", "-", Number::MAX);
        assert!(coordinate.nudge_backward().is_none());
    }
}
