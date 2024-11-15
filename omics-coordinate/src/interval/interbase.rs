//! Interbase intervals.

use crate::Strand;
use crate::base;
use crate::interbase::Coordinate;
use crate::interval::Number;
use crate::interval::r#trait;
use crate::system::Interbase;

////////////////////////////////////////////////////////////////////////////////////////
// Intervals
////////////////////////////////////////////////////////////////////////////////////////

/// An interbase interval.
///
/// Interbase intervals consist of two interbase positions. The range is
/// represented by the interval `[start, end]`.
///
/// It is worth noting that, historically, the translation of interbase
/// positions to the nucleotide encoded in the range have included a heuristic
/// whereby the following nucleotide was associated with the interbase position.
/// This is why interbase intervals are often described as _exclusively_ bounded
/// at the end, as the nucleotide _following_ the interbase position is not
/// included in the range (and, generally, what one is most interested in is the
/// nucleotides or other entities contained therein).
///
/// Because this crate does not co-mingle the idea of interbase and in-base
/// representations, no such gymnastics are required.
pub type Interval = crate::Interval<Interbase>;

impl Interval {
    /// Checks whether the interbase interval contains the entity _after_ the
    /// specified interbase coordinate.
    ///
    /// This method returns an [`Option`] because the next coordinate may or may
    /// not be a valid position.
    pub fn contains_next_entity(&self, coordinate: Coordinate) -> Option<bool> {
        let coordinate = coordinate.nudge_forward()?;
        Some(self.contains_entity(&coordinate))
    }

    /// Checks whether the interbase interval contains the entity _before_ the
    /// specified interbase coordinate.
    ///
    /// This method returns an [`Option`] because the next coordinate may or may
    /// not be a valid position.
    pub fn contains_prev_entity(&self, coordinate: Coordinate) -> Option<bool> {
        let coordinate = coordinate.nudge_backward()?;
        Some(self.contains_entity(&coordinate))
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl r#trait::Interval<Interbase> for Interval {
    fn contains_entity(&self, coordinate: &base::Coordinate) -> bool {
        if self.contig() != coordinate.contig() {
            return false;
        }

        if self.strand() != coordinate.strand() {
            return false;
        }

        match self.strand() {
            Strand::Positive => {
                self.start().position().get() < coordinate.position().get()
                    && self.end().position().get() >= coordinate.position().get()
            }
            Strand::Negative => {
                self.start().position().get() >= coordinate.position().get()
                    && self.end().position().get() < coordinate.position().get()
            }
        }
    }

    /// Gets the number of entities within the interval.
    fn count_entities(&self) -> Number {
        self.start()
            .position()
            .distance_unchecked(self.end().position())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Coordinate;
    use crate::system::Base;

    fn create_coordinate(
        contig: &str,
        strand: &str,
        position: Number,
    ) -> crate::Coordinate<Interbase> {
        Coordinate::try_new(contig, strand, position).unwrap()
    }

    fn create_base_coordinate(
        contig: &str,
        strand: &str,
        position: Number,
    ) -> crate::Coordinate<Base> {
        Coordinate::try_new(contig, strand, position).unwrap()
    }

    fn create_interval(contig: &str, strand: &str, start: Number, end: Number) -> Interval {
        Interval::try_new(
            create_coordinate(contig, strand, start),
            create_coordinate(contig, strand, end),
        )
        .unwrap()
    }

    #[test]
    fn contains() {
        let interval = create_interval("seq0", "+", 10, 20);

        // An interval contains the coordinate representing its start position.
        assert!(interval.contains_coordinate(interval.start()));

        // An interval contains the coordinate representing its end position.
        assert!(interval.contains_coordinate(interval.end()));

        // An interval contains a coordinate in the middle of its range.
        assert!(interval.contains_coordinate(&create_coordinate("seq0", "+", 15)));

        // An interval does not contain the position _before_ its start
        // position.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "+", 9)));

        // An interval does not contain the position _after_ its end position.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "+", 21)));

        // An interval does not contain a random other position.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "+", 1000)));

        // An interval does not contain a coordinate on another contig.
        assert!(!interval.contains_coordinate(&create_coordinate("seq1", "+", 15)));

        // An interval does not contain a coordinate on another strand.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "-", 15)));

        let interval = create_interval("seq0", "-", 20, 10);

        // An interval contains the coordinate representing its start position.
        assert!(interval.contains_coordinate(interval.start()));

        // An interval contains the coordinate representing its end position.
        assert!(interval.contains_coordinate(interval.end()));

        // An interval contains a coordinate in the middle of its range.
        assert!(interval.contains_coordinate(&create_coordinate("seq0", "-", 15)));

        // An interval does not contain the position _before_ its start
        // position.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "-", 21)));

        // An interval does not contain the position _after_ its end position.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "-", 9)));

        // An interval does not contain a random other position.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "-", 0)));

        // An interval does not contain a coordinate on another contig.
        assert!(!interval.contains_coordinate(&create_coordinate("seq1", "-", 15)));

        // An interval does not contain a coordinate on another strand.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "+", 15)));
    }

    #[test]
    fn contains_entity() {
        let interval = create_interval("seq0", "+", 10, 20);

        // An interbase interval does not contain the entity at the half-step
        // _before_ its start position when on the positive strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq0", "+", 10)));

        // An interbase interval does contain the entity at the half-step
        // _after_ its start position when on the positive strand.
        assert!(interval.contains_entity(&create_base_coordinate("seq0", "+", 11)));

        // An interbase interval does contain the entity at the half-step
        // _before_ its end position when on the positive strand.
        assert!(interval.contains_entity(&create_base_coordinate("seq0", "+", 20)));

        // An interbase interval does not contain the entity at the half-step
        // _after_ its end position when on the positive strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq0", "+", 21)));

        // An interbase interval does contain an entity midway through its range
        // on the positive strand.
        assert!(interval.contains_entity(&create_base_coordinate("seq0", "+", 15)));

        // An interbase interval does not contain an entity on a different
        // contig on the positive strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq1", "+", 15)));

        // An interbase interval does not contain an entity on a different
        // strand on the positive strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq0", "-", 15)));

        let interval = create_interval("seq0", "-", 20, 10);

        // An interbase interval does not contain the entity at the half-step
        // _before_ its start position when on the negative strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq0", "-", 21)));

        // An interbase interval does contain the entity at the half-step
        // _after_ its start position when on the negative strand.
        assert!(interval.contains_entity(&create_base_coordinate("seq0", "-", 20)));

        // An interbase interval does contain the entity at the half-step
        // _before_ its end position when on the negative strand.
        assert!(interval.contains_entity(&create_base_coordinate("seq0", "-", 11)));

        // An interbase interval does not contain the entity at the half-step
        // _after_ its end position when on the negative strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq0", "-", 10)));

        // An interbase interval does contain an entity midway through its range
        // on the negative strand.
        assert!(interval.contains_entity(&create_base_coordinate("seq0", "-", 15)));

        // An interbase interval does not contain an entity on a different
        // contig on the negative strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq1", "-", 15)));

        // An interbase interval does not contain an entity on a different
        // strand on the negative strand.
        assert!(!interval.contains_entity(&create_base_coordinate("seq0", "+", 15)));
    }

    #[test]
    fn contains_next_entity() {
        let interval = create_interval("seq0", "+", 10, 20);

        // An interbase interval does not contain the next entity after its start
        // position is moved backwards one on the positive strand.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq0", "+", 9))
                .unwrap()
        );

        // An interbase interval does contain the next entity after its start
        // position on the positive strand.
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "+", 10))
                .unwrap()
        );

        // An interbase interval does contain the next entity after its end
        // position is moved backwards one on the positive strand.
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "+", 19))
                .unwrap()
        );

        // An interbase interval does not contain the next entity after its end
        // position on the positive strand.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq0", "+", 20))
                .unwrap()
        );

        // An interbase interval does contain the next entity after a position
        // in the middle of its range.
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "+", 15))
                .unwrap()
        );

        // An interbase interval does not contain the next entity for a
        // coordinate on a different contig.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq1", "+", 15))
                .unwrap()
        );

        // An interbase interval does not contain the next entity for a
        // coordinate on a different strand.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq0", "-", 15))
                .unwrap()
        );

        let interval = create_interval("seq0", "-", 20, 10);

        // An interbase interval does not contain the next entity after its start
        // position is moved backwards one on the negative strand.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq0", "-", 21))
                .unwrap()
        );

        // An interbase interval does contain the next entity after its start
        // position on the negative strand.
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "-", 20))
                .unwrap()
        );

        // An interbase interval does contain the next entity after its end
        // position is moved backwards on the negative strand.
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "-", 11))
                .unwrap()
        );

        // An interbase interval does not contain the next entity after its end
        // position on the negative strand.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq0", "-", 10))
                .unwrap()
        );

        // An interbase interval does contain the next entity after a position
        // in the middle of its range.
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "-", 15))
                .unwrap()
        );

        // An interbase interval does not contain the next entity for a
        // coordinate on a different contig.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq1", "-", 15))
                .unwrap()
        );

        // An interbase interval does not contain the next entity for a
        // coordinate on a different strand.
        assert!(
            !interval
                .contains_next_entity(create_coordinate("seq0", "+", 15))
                .unwrap()
        );

        let interval = create_interval("seq0", "+", Number::MAX - 10, Number::MAX);

        // This position should fail to bump forward a half-step, so the entire
        // operation should return a [`None`].
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "+", Number::MAX))
                .is_none()
        );

        let interval = create_interval("seq0", "-", 10, 0);

        // This position should fail to bump forward a half-step, so the entire
        // operation should return a [`None`].
        assert!(
            interval
                .contains_next_entity(create_coordinate("seq0", "-", 0))
                .is_none()
        );
    }

    #[test]
    fn contains_prev_entity() {
        let interval = create_interval("seq0", "+", 10, 20);

        // An interbase interval does not contain the previous entity before its start
        // position on the positive strand.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq0", "+", 10))
                .unwrap()
        );

        // An interbase interval does contain the previous entity before its start
        // position is moved forward by one on the positive strand.
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "+", 11))
                .unwrap()
        );

        // An interbase interval does contain the previous entity before its end
        // position on the positive strand.
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "+", 20))
                .unwrap()
        );

        // An interbase interval does not contain the previous entity before its end
        // position is moved forward by one on the positive strand.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq0", "+", 21))
                .unwrap()
        );

        // An interbase interval does contain the previous entity before a position
        // in the middle of its range.
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "+", 15))
                .unwrap()
        );

        // An interbase interval does not contain the previous entity for a
        // coordinate on a different contig.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq1", "+", 15))
                .unwrap()
        );

        // An interbase interval does not contain the previous entity for a
        // coordinate on a different strand.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq0", "-", 15))
                .unwrap()
        );

        let interval = create_interval("seq0", "-", 20, 10);

        // An interbase interval does not contain the previous entity before its start
        // position on the negative strand.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq0", "-", 20))
                .unwrap()
        );

        // An interbase interval does contain the previous entity before its start
        // position is moved forward by one on the negative strand.
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "-", 19))
                .unwrap()
        );

        // An interbase interval does contain the previous entity before its end
        // position on the negative strand.
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "-", 10))
                .unwrap()
        );

        // An interbase interval does not contain the previous entity before its end
        // position is moved forward by one on the negative strand.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq0", "-", 9))
                .unwrap()
        );

        // An interbase interval does contain the previous entity before a position
        // in the middle of its range.
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "-", 15))
                .unwrap()
        );

        // An interbase interval does not contain the previous entity for a
        // coordinate on a different contig.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq1", "-", 15))
                .unwrap()
        );

        // An interbase interval does not contain the previous entity for a
        // coordinate on a different strand.
        assert!(
            !interval
                .contains_prev_entity(create_coordinate("seq0", "+", 15))
                .unwrap()
        );

        let interval = create_interval("seq0", "+", 0, 10);

        // This position should fail to bump backwards a half-step, so the entire
        // operation should return a [`None`].
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "+", 0))
                .is_none()
        );

        let interval = create_interval("seq0", "-", Number::MAX, Number::MAX - 10);

        // This position should fail to bump forward a half-step, so the entire
        // operation should return a [`None`].
        assert!(
            interval
                .contains_prev_entity(create_coordinate("seq0", "-", Number::MAX))
                .is_none()
        );
    }
}
