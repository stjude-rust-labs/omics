//! Base intervals.

use crate::base::Coordinate;
use crate::interval::r#trait;
use crate::position::Number;
use crate::span::Span;
use crate::system::Base;
use crate::system::Interbase;

////////////////////////////////////////////////////////////////////////////////////////
// Intervals
////////////////////////////////////////////////////////////////////////////////////////

/// A base interval.
///
/// Base intervals consist of two in-base positions. The range is represented by
/// the interval `[start, end]`.
pub type Interval = crate::Interval<Base>;

impl Interval {
    /// Consumes `self` and returns the equivalent interbase interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::system::Base;
    /// use omics_coordinate::system::Interbase;
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<Base>>()?;
    /// let equivalent = interval.into_equivalent_interbase();
    ///
    /// assert_eq!("seq0:+:0-1000".parse::<Interval<Interbase>>()?, equivalent);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_equivalent_interbase(self) -> crate::interval::Interval<Interbase> {
        let (start, end) = self.into_coordinates();

        // SAFETY: every valid base interval endpoint has the required interbase
        // boundary.
        let start = start.nudge_backward().unwrap();
        // SAFETY: every valid base interval endpoint has the required interbase
        // boundary.
        let end = end.nudge_forward().unwrap();
        let (contig, strand, start) = start.into_parts();
        let (_, _, end) = end.into_parts();
        let span = Span::from((start, end));
        let strand = match strand {
            crate::Strand::Positive => "+",
            crate::Strand::Negative => "-",
        };

        // SAFETY: the converted span retains the direction required by the strand.
        crate::interval::Interval::<Interbase>::try_new(contig.as_str(), strand, span).unwrap()
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Trait implementations
////////////////////////////////////////////////////////////////////////////////////////

impl r#trait::Interval<Base> for Interval {
    fn contains_entity(&self, coordinate: &Coordinate) -> bool {
        self.contig() == coordinate.contig()
            && self.strand() == coordinate.strand()
            && self.span().contains_entity(coordinate.position())
    }

    /// Gets the number of entities within the interval.
    fn count_entities(&self) -> Number {
        self.span().count_entities()
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
        let interval = Interval::try_from((
            create_coordinate(contig, strand, start),
            create_coordinate(contig, strand, end),
        ));
        // SAFETY: the test helper receives endpoints compatible with the requested
        // strand.
        interval.unwrap()
    }

    #[test]
    fn contains() {
        let interval = create_interval("seq0", "+", 10, 20);

        // An interval contains the coordinate representing its start position.
        assert!(interval.contains_coordinate(&interval.start().into_owned()));

        // An interval contains the coordinate representing its end position.
        assert!(interval.contains_coordinate(&interval.end().into_owned()));

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
        assert!(interval.contains_coordinate(&interval.start().into_owned()));

        // An interval contains the coordinate representing its end position.
        assert!(interval.contains_coordinate(&interval.end().into_owned()));

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
