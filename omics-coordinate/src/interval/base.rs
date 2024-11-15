//! Base intervals.

use crate::base::Coordinate;
use crate::interval::r#trait;
use crate::position::Number;
use crate::system::Base;

////////////////////////////////////////////////////////////////////////////////////////
// Intervals
////////////////////////////////////////////////////////////////////////////////////////

/// A base interval.
///
/// Base intervals consist of two in-base positions. The range is represented by
/// the interval `[start, end]`.
pub type Interval = crate::Interval<Base>;

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl r#trait::Interval<Base> for Interval {
    fn contains_entity(&self, coordinate: &Coordinate) -> bool {
        // NOTE: for in-base positions whether or not the entity is contained in
        // the interval matches the implementation of whether or not the
        // coordinate is contained within the interval, as entities and
        // coordinates are essentially the same thing in this system.
        self.contains_coordinate(coordinate)
    }

    /// Gets the number of entities within the interval.
    fn count_entities(&self) -> Number {
        self.start()
            .position()
            .distance_unchecked(self.end().position())
            + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Coordinate;
    use crate::system::Base;

    fn create_coordinate(contig: &str, strand: &str, position: Number) -> crate::Coordinate<Base> {
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
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "-", 1)));

        // An interval does not contain a coordinate on another contig.
        assert!(!interval.contains_coordinate(&create_coordinate("seq1", "-", 15)));

        // An interval does not contain a coordinate on another strand.
        assert!(!interval.contains_coordinate(&create_coordinate("seq0", "+", 15)));
    }

    #[test]
    fn contains_entity() {
        let interval = create_interval("seq0", "+", 10, 20);

        // A base interval does not contain the entity at the full-step _before_
        // its start position when on the positive strand.
        assert!(!interval.contains_entity(&create_coordinate("seq0", "+", 9)));

        // A base interval does contain the entity at it start position when on
        // the positive strand.
        assert!(interval.contains_entity(&create_coordinate("seq0", "+", 10)));

        // A base interval does contain the entity at the its end position when
        // on the positive strand.
        assert!(interval.contains_entity(&create_coordinate("seq0", "+", 20)));

        // A base interval does not contain the entity at the full-step _after_
        // its end position when on the positive strand.
        assert!(!interval.contains_entity(&create_coordinate("seq0", "+", 21)));

        // A base interval does contain an entity midway through its range on
        // the positive strand.
        assert!(interval.contains_entity(&create_coordinate("seq0", "+", 15)));

        // A base interval does not contain an entity on a different contig on
        // the positive strand.
        assert!(!interval.contains_entity(&create_coordinate("seq1", "+", 15)));

        // A base interval does not contain an entity on a different strand on
        // the positive strand.
        assert!(!interval.contains_entity(&create_coordinate("seq0", "-", 15)));

        let interval = create_interval("seq0", "-", 20, 10);

        // A base interval does not contain the entity at the full-step _before_
        // its start position when on the negative strand.
        assert!(!interval.contains_entity(&create_coordinate("seq0", "-", 21)));

        // A base interval does contain the entity at its start position when on
        // the negative strand.
        assert!(interval.contains_entity(&create_coordinate("seq0", "-", 20)));

        // A base interval does contain the entity at its end position when on
        // the negative strand.
        assert!(interval.contains_entity(&create_coordinate("seq0", "-", 10)));

        // A base interval does not contain the entity at the full-step _after_
        // its end position when on the negative strand.
        assert!(!interval.contains_entity(&create_coordinate("seq0", "-", 9)));

        // A base interval does contain an entity midway through its range on
        // the negative strand.
        assert!(interval.contains_entity(&create_coordinate("seq0", "-", 15)));

        // A base interval does not contain an entity on a different contig on
        // the negative strand.
        assert!(!interval.contains_entity(&create_coordinate("seq1", "-", 15)));

        // A base interval does not contain an entity on a different strand on
        // the negative strand.
        assert!(!interval.contains_entity(&create_coordinate("seq0", "+", 15)));
    }
}
