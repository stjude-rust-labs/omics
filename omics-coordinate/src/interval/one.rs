//! A 1-based, fully-closed interval consisting of a start and stop
//! [`Coordinate`].
//!
//! Intervals can be either on the positive strand or the negative strand of a
//! contig. Below is an example of both a positive- and negative-stranded,
//! [`One`]-based interval on `seq0`.
//!
//! ```text
//! ================ seq0 ===============
//!
//! | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
//! -------------------------------------
//! |   |   | X | X | X | X | X |   |   |  <= seq0:+:3-7
//! |   |   |   |   | X | X | X | X | X |  <= seq0:-:9-5
//! ```
//! In this figure, `X` represents a position that is included in the range
//! while `O` represents the non-inclusive bound:
//!
//! - The first interval above (`seq0:3-7`) is a forward-stranded interval from
//!   3 up to and including 7.
//! - The second interval above (`seq0:9-5`) is a reverse-stranded interval from
//!   9 down to and including 5.

use crate::Strand;
use crate::interval;
use crate::interval::Error;
use crate::one::Coordinate;
use crate::system::One;

/// A 1-based, fully-closed [`Interval`](crate::Interval).
pub type Interval = crate::Interval<One>;

impl interval::r#trait::Interval<One> for Interval {
    /// Attempts to create a new 1-based, fully-closed [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::One;
    ///
    /// // Positive-stranded interval
    ///
    /// let start = Coordinate::<One>::try_new("seq0", Strand::Positive, 1)?;
    /// let end = Coordinate::<One>::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// // Negative-stranded interval
    ///
    /// let start = Coordinate::<One>::try_new("seq0", Strand::Negative, 1000)?;
    /// let end = Coordinate::<One>::try_new("seq0", Strand::Negative, 1)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_new(start: Coordinate, end: Coordinate) -> Result<Interval, Error> {
        interval::try_new(start, end)
    }

    /// Indicates whether a [`Coordinate`] falls within an [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::One;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:+:1-1000".parse::<Interval<One>>()?;
    ///
    /// assert!(interval.contains(&"seq0:+:1".parse::<Coordinate<One>>()?));
    /// assert!(interval.contains(&"seq0:+:1000".parse::<Coordinate<One>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:+:1".parse::<Coordinate<One>>()?));
    /// assert!(!interval.contains(&"seq0:-:50".parse::<Coordinate<One>>()?));
    /// assert!(!interval.contains(&"seq0:+:1001".parse::<Coordinate<One>>()?));
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:-:1000-1".parse::<Interval<One>>()?;
    ///
    /// assert!(interval.contains(&"seq0:-:1000".parse::<Coordinate<One>>()?));
    /// assert!(interval.contains(&"seq0:-:1".parse::<Coordinate<One>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:-:1000".parse::<Coordinate<One>>()?));
    /// assert!(!interval.contains(&"seq0:+:1000".parse::<Coordinate<One>>()?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn contains(&self, coordinate: &Coordinate) -> bool {
        if self.contig() != coordinate.contig() {
            return false;
        }

        if self.strand() != coordinate.strand() {
            return false;
        }

        match self.strand() {
            Strand::Positive => {
                if self.start().position() > coordinate.position() {
                    return false;
                }

                if self.end().position() < coordinate.position() {
                    return false;
                }
            }
            Strand::Negative => {
                if self.start().position() < coordinate.position() {
                    return false;
                }

                if self.end().position() > coordinate.position() {
                    return false;
                }
            }
        }

        true
    }

    /// Gets the length of the [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::One;
    ///
    /// // Positive-stranded interval
    ///
    /// let start = "seq0:+:1".parse::<Coordinate<One>>()?;
    /// let end = "seq0:+:1000".parse::<Coordinate<One>>()?;
    ///
    /// let interval = Interval::<One>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// // Negative-stranded interval
    ///
    /// let start = "seq0:-:1000".parse::<Coordinate<One>>()?;
    /// let end = "seq0:-:1".parse::<Coordinate<One>>()?;
    ///
    /// let interval = Interval::<One>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn len(&self) -> usize {
        // SAFETY: For the first unwrap, because of the bounds checks we do when
        // creating an [`Interval`], the resulting [`Position`] should always
        // unwrap.
        self.start()
            .position()
            // NOTE: for 1-based, fully closed intervals, you must add one after
            // subtracting the start from the end.
            .distance_unchecked(self.end().position())
            .and_then(|distance| distance.checked_add(1))
            .unwrap()
    }

    fn complement(self) -> Result<Option<interval::Interval<One>>, Error> {
        // Because 1-based, fully-closed intervals are inclusive, the inversion
        // is quite simple: one can simply switch the start and end coordinates.
        let (start, end) = self.into_coordinates();

        Interval::try_new(
            end.swap_strand().map_err(Error::Coordinate)?,
            start.swap_strand().map_err(Error::Coordinate)?,
        )
        .map(Some)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interval::r#trait::Interval as _;
    use crate::position::Value;

    #[test]
    fn a_zero_sized_interval_evaluates_correctly() -> Result<(), Box<dyn std::error::Error>> {
        let position = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let interval = Interval::try_new(position.clone(), position.clone())?;

        assert_eq!(interval.start(), &position);
        assert_eq!(interval.end(), &position);
        assert_eq!(interval.len(), 1);

        Ok(())
    }

    #[test]
    fn complement_a_one_based_interval_works_correctly() -> Result<(), Box<dyn std::error::Error>> {
        // A positive-stranded, zero-sized interval.
        let position = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let interval = Interval::try_new(position.clone(), position.clone())?;
        let complemented = interval.complement()?.unwrap();

        assert_eq!(complemented.contig().inner(), "seq0");
        assert_eq!(complemented.strand(), &Strand::Negative);
        assert_eq!(complemented.start().position().inner(), &Value::Usize(1));
        assert_eq!(complemented.end().position().inner(), &Value::Usize(1));

        // A positive-stranded interval with some magnitude.
        let start = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let end = Coordinate::try_new("seq0", Strand::Positive, 10)?;
        let complemented = Interval::try_new(start, end)?.complement()?.unwrap();

        assert_eq!(complemented.contig().inner(), "seq0");
        assert_eq!(complemented.strand(), &Strand::Negative);
        assert_eq!(complemented.start().position().inner(), &Value::Usize(10));
        assert_eq!(complemented.end().position().inner(), &Value::Usize(1));

        // A negative-stranded, zero-sized interval.
        let position = Coordinate::try_new("seq0", Strand::Negative, 1)?;
        let interval = Interval::try_new(position.clone(), position.clone())?;
        let complemented = interval.complement()?.unwrap();

        assert_eq!(complemented.contig().inner(), "seq0");
        assert_eq!(complemented.strand(), &Strand::Positive);
        assert_eq!(complemented.start().position().inner(), &Value::Usize(1));
        assert_eq!(complemented.end().position().inner(), &Value::Usize(1));

        // A positive-stranded interval with some magnitude.
        let start = Coordinate::try_new("seq0", Strand::Negative, 10)?;
        let end = Coordinate::try_new("seq0", Strand::Negative, 1)?;
        let complemented = Interval::try_new(start, end)?.complement()?.unwrap();

        assert_eq!(complemented.contig().inner(), "seq0");
        assert_eq!(complemented.strand(), &Strand::Positive);
        assert_eq!(complemented.start().position().inner(), &Value::Usize(1));
        assert_eq!(complemented.end().position().inner(), &Value::Usize(10));

        // The maximum positive-stranded interval.
        let start = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let end = Coordinate::try_new("seq0", Strand::Positive, usize::MAX)?;
        let complemented = Interval::try_new(start, end)?.complement()?.unwrap();

        assert_eq!(complemented.contig().inner(), "seq0");
        assert_eq!(complemented.strand(), &Strand::Negative);
        assert_eq!(
            complemented.start().position().inner(),
            &Value::Usize(usize::MAX)
        );
        assert_eq!(complemented.end().position().inner(), &Value::Usize(1));

        // The maximum negative-stranded interval.
        let start = Coordinate::try_new("seq0", Strand::Negative, usize::MAX)?;
        let end = Coordinate::try_new("seq0", Strand::Negative, 1)?;
        let complemented = Interval::try_new(start, end)?.complement()?.unwrap();

        assert_eq!(complemented.contig().inner(), "seq0");
        assert_eq!(complemented.strand(), &Strand::Positive);
        assert_eq!(complemented.start().position().inner(), &Value::Usize(1));
        assert_eq!(
            complemented.end().position().inner(),
            &Value::Usize(usize::MAX)
        );

        Ok(())
    }
}
