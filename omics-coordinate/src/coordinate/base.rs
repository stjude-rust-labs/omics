//! Convenient access to a base [`Coordinate`](crate::Coordinate).

use crate::Strand;
use crate::contig;
use crate::coordinate::Error;
use crate::position;
use crate::position::base::Position;
use crate::strand;
use crate::system::Base;
use crate::system::Interbase;

/// A base coordinate.
pub type Coordinate = crate::Coordinate<Base>;

impl Coordinate {
    /// Consumes `self` and attempts to convert the coordinate to the next
    /// interbase coordinate.
    ///
    /// This is only useful when converting a coordinate to the other coordinate
    /// system, which is not a common operation. If, instead, you wish to move
    /// the coordinate forward within the _same_ coordinate system, you're
    /// almost certainly looking for
    /// [`move_forward()`](crate::Coordinate::move_forward).
    pub fn nudge_forward(self) -> Option<crate::Coordinate<Interbase>> {
        let (contig, strand, position) = self.into_parts();

        let position = match strand {
            Strand::Positive => crate::Position::<Interbase>::new(position.get()),
            Strand::Negative => position
                .get()
                .checked_sub(1)
                .map(crate::Position::<Interbase>::new)?,
        };

        Some(crate::Coordinate::<Interbase>::new(
            contig, strand, position,
        ))
    }

    /// Consumes `self` and attempts to convert the coordinate to the previous
    /// interbase coordinate.
    ///
    /// This is only useful when converting a coordinate to the other coordinate
    /// system, which is not a common operation. If, instead, you wish to move
    /// the coordinate backward within the _same_ coordinate system, you're
    /// almost certainly looking for
    /// [`move_backward()`](crate::Coordinate::move_backward).
    pub fn nudge_backward(self) -> Option<crate::Coordinate<Interbase>> {
        let (contig, strand, position) = self.into_parts();

        let position = match strand {
            Strand::Positive => position
                .get()
                .checked_sub(1)
                .map(crate::Position::<Interbase>::new)?,
            Strand::Negative => crate::Position::<Interbase>::new(position.get()),
        };

        Some(crate::Coordinate::<Interbase>::new(
            contig, strand, position,
        ))
    }
}

impl crate::coordinate::r#trait::Coordinate<Base> for Coordinate {
    fn try_new(
        contig: impl TryInto<crate::Contig, Error = contig::Error>,
        strand: impl TryInto<crate::Strand, Error = strand::Error>,
        position: position::Number,
    ) -> super::Result<Self> {
        let contig = contig.try_into().map_err(Error::Contig)?;
        let strand = strand.try_into().map_err(Error::Strand)?;
        let position = Position::try_new(position).map_err(Error::Position)?;

        Ok(Self {
            system: Base,
            contig,
            strand,
            position,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interbase;
    use crate::position::Number;

    fn create_coordinate(contig: &str, strand: &str, position: Number) -> Coordinate {
        Coordinate::try_new(contig, strand, position).unwrap()
    }

    fn create_interbase_coordinate(
        contig: &str,
        strand: &str,
        position: Number,
    ) -> interbase::Coordinate {
        interbase::Coordinate::try_new(contig, strand, position).unwrap()
    }

    #[test]
    fn nudge_forward() {
        let coordinate = create_coordinate("seq0", "+", 1);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_interbase_coordinate("seq0", "+", 1)
        );

        let coordinate = create_coordinate("seq0", "+", 10);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_interbase_coordinate("seq0", "+", 10)
        );

        let coordinate = create_coordinate("seq0", "+", Number::MAX);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_interbase_coordinate("seq0", "+", Number::MAX)
        );

        let coordinate = create_coordinate("seq0", "-", 1);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_interbase_coordinate("seq0", "-", 0)
        );

        let coordinate = create_coordinate("seq0", "-", 10);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_interbase_coordinate("seq0", "-", 9)
        );

        let coordinate = create_coordinate("seq0", "-", Number::MAX);
        assert_eq!(
            coordinate.nudge_forward().unwrap(),
            create_interbase_coordinate("seq0", "-", Number::MAX - 1)
        );
    }

    #[test]
    fn nudge_backward() {
        let coordinate = create_coordinate("seq0", "+", 1);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_interbase_coordinate("seq0", "+", 0)
        );

        let coordinate = create_coordinate("seq0", "+", 10);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_interbase_coordinate("seq0", "+", 9)
        );

        let coordinate = create_coordinate("seq0", "+", Number::MAX);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_interbase_coordinate("seq0", "+", Number::MAX - 1)
        );

        let coordinate = create_coordinate("seq0", "-", 1);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_interbase_coordinate("seq0", "-", 1)
        );

        let coordinate = create_coordinate("seq0", "-", 10);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_interbase_coordinate("seq0", "-", 10)
        );

        let coordinate = create_coordinate("seq0", "-", Number::MAX);
        assert_eq!(
            coordinate.nudge_backward().unwrap(),
            create_interbase_coordinate("seq0", "-", Number::MAX)
        );
    }
}
