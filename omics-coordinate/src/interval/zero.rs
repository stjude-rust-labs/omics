//! A 0-based, half-open interval consisting of a start and stop [`Coordinate`].
//!
//! Intervals can be either on the positive strand or the negative strand of a
//! contig. Below is an example of both a positive- and negative-stranded,
//! [`Zero`]-based interval on `seq0`.
//!
//! ```text
//! ================== seq0 =================
//!
//! | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
//! -----------------------------------------
//! |   |   |   | X | X | X | X | O |   |   |  <= seq0:+:3-7
//! |   |   |   |   |   | O | X | X | X | X |  <= seq0:-:9-5
//! ```
//! In this figure, `X` represents a position that is included in the range
//! while `O` represents the non-inclusive bound:
//!
//! - The first interval above (`seq0:3-7`) is a forward-stranded interval from
//!   3 up until (but not including) 7.
//! - The second interval above (`seq0:9-5`) is a reverse-stranded interval from
//!   9 down until (but not including) 5.

use crate::Strand;
use crate::interval;
use crate::interval::Error;
use crate::position::Value;
use crate::system::Zero;
use crate::zero::Coordinate;

/// A 0-based, half-open [`Interval`](crate::Interval).
pub type Interval = crate::Interval<Zero>;

impl crate::interval::r#trait::Interval<Zero> for Interval {
    /// Attempts to create a new 0-based, half-open [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::interval::zero::Interval;
    /// use omics_coordinate::zero::Coordinate;
    ///
    /// // Positive-stranded interval
    ///
    /// let start = Coordinate::try_new("seq0", Strand::Positive, 0)?;
    /// let end = Coordinate::try_new("seq0", Strand::Positive, 1000)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// // Negative-stranded interval
    ///
    /// let start = Coordinate::try_new("seq0", Strand::Negative, 1000)?;
    /// let end = Coordinate::try_new("seq0", Strand::Negative, 0)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    fn try_new(start: Coordinate, end: Coordinate) -> Result<Interval, Error> {
        let result = interval::try_new(start, end)?;

        // For 0-based intervals only, the start and end coordinate cannot be
        // the same. This is because the end coordinate is non-inclusive, which
        // conflicts with the start coordinate being inclusive.
        if result.start() == result.end() {
            return Err(Error::ZeroSizedInterval);
        }

        Ok(result)
    }

    /// Indicates whether a [`Coordinate`] falls within an [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Interval;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:+:0-1000".parse::<Interval<Zero>>()?;
    ///
    /// assert!(interval.contains(&"seq0:+:0".parse::<Coordinate<Zero>>()?));
    /// assert!(interval.contains(&"seq0:+:999".parse::<Coordinate<Zero>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:+:0".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:-:50".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:+:1000".parse::<Coordinate<Zero>>()?));
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:-:1000-0".parse::<Interval<Zero>>()?;
    ///
    /// assert!(interval.contains(&"seq0:-:1000".parse::<Coordinate<Zero>>()?));
    /// assert!(interval.contains(&"seq0:-:1".parse::<Coordinate<Zero>>()?));
    ///
    /// assert!(!interval.contains(&"seq1:-:1000".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:+:1000".parse::<Coordinate<Zero>>()?));
    /// assert!(!interval.contains(&"seq0:-:0".parse::<Coordinate<Zero>>()?));
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

                if self.end().position() <= coordinate.position() {
                    return false;
                }
            }
            Strand::Negative => {
                if self.start().position() < coordinate.position() {
                    return false;
                }

                if self.end().position() >= coordinate.position() {
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
    /// use omics_coordinate::system::Zero;
    ///
    /// // Positive-stranded interval
    ///
    /// let start = "seq0:+:0".parse::<Coordinate<Zero>>()?;
    /// let end = "seq0:+:1000".parse::<Coordinate<Zero>>()?;
    ///
    /// let interval = Interval::<Zero>::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// // Negative-stranded interval
    ///
    /// let start = "seq0:-:1000".parse::<Coordinate<Zero>>()?;
    /// let end = "seq0:-:0".parse::<Coordinate<Zero>>()?;
    ///
    /// let interval = Interval::<Zero>::try_new(start, end)?;
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
            // NOTE: for 0-based, half-open intervals, the distance is simply
            // the end minus the start—no need to add one after that
            // calculation.
            .distance_unchecked(self.end().position())
            .unwrap()
    }

    fn complement(self) -> Result<Option<interval::Interval<Zero>>, Error> {
        // Because 0-based, half-open intervals are not inclusive of the last
        // element, both the start and end positions needs to (a) swap to the
        // opposite strand and (b) move forward by a magnitude of one to shift
        // the bounding coordinates to the correct location.
        let (current_start, current_end) = self.into_coordinates();

        let new_start = match current_end.position().inner() {
            // If the current end's position is the lower bound, then we know
            // the complemented position will be zero on the positive strand.
            //
            // This case has to be handled specially because moving to the
            // opposite strand and moving forward by a magnitude of one are not
            // an atomic operation. As such, the code will error if it tries to
            // complement the lower bound before it can move forward.
            //
            // Notably, if you change the order of operations so that the
            // position moves before complementing the strand, other cases
            // break.

            // SAFETY: this coordinate is hand crafted to always unwrap.
            Value::LowerBound => {
                Coordinate::try_new(current_end.contig().clone(), Strand::Positive, 0).unwrap()
            }
            // SAFETY: for the first unwrap, the only case that can fail
            // swapping strands is handled separately above.
            Value::Usize(_) => {
                let coordinate = current_end
                    .swap_strand()
                    .map_err(Error::Coordinate)?
                    .move_forward(1)
                    .map_err(Error::Coordinate)?;

                match coordinate {
                    Some(coordinate) => coordinate,
                    None => return Ok(None),
                }
            }
        };

        let new_end = match current_start
            .swap_strand()
            .map_err(Error::Coordinate)?
            .move_forward(1)
            .map_err(Error::Coordinate)?
        {
            Some(end) => end,
            None => return Ok(None),
        };

        Interval::try_new(new_start, new_end).map(Some)
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::interval::r#trait::Interval as _;

    #[test]
    fn complementing_a_positive_interval_with_two_usize_positions_works()
    -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let end = Coordinate::try_new("seq0", Strand::Positive, 6)?;

        // This represents the interval [1, 6), which should be complemented to
        // (0, 5] when it's all said and done.
        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.start().contig().inner(), "seq0");
        assert_eq!(interval.start().strand(), &Strand::Positive);
        assert_eq!(interval.start().position().inner(), &Value::Usize(1));
        assert_eq!(interval.end().contig().inner(), "seq0");
        assert_eq!(interval.end().strand(), &Strand::Positive);
        assert_eq!(interval.end().position().inner(), &Value::Usize(6));

        let complement = interval.complement()?.unwrap();
        assert_eq!(complement.start().contig().inner(), "seq0");
        assert_eq!(complement.start().strand(), &Strand::Negative);
        assert_eq!(complement.start().position().inner(), &Value::Usize(5));
        assert_eq!(complement.end().contig().inner(), "seq0");
        assert_eq!(complement.end().strand(), &Strand::Negative);
        assert_eq!(complement.end().position().inner(), &Value::Usize(0));

        // Testing converting to lower bound.
        let start = Coordinate::try_new("seq0", Strand::Positive, 0)?;
        let end = Coordinate::try_new("seq0", Strand::Positive, 6)?;

        // This represents the interval [0, 6), which should be complemented to
        // (LB, 5] when it's all said and done.
        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.start().contig().inner(), "seq0");
        assert_eq!(interval.start().strand(), &Strand::Positive);
        assert_eq!(interval.start().position().inner(), &Value::Usize(0));
        assert_eq!(interval.end().contig().inner(), "seq0");
        assert_eq!(interval.end().strand(), &Strand::Positive);
        assert_eq!(interval.end().position().inner(), &Value::Usize(6));

        let complement = interval.complement()?.unwrap();
        assert_eq!(complement.start().contig().inner(), "seq0");
        assert_eq!(complement.start().strand(), &Strand::Negative);
        assert_eq!(complement.start().position().inner(), &Value::Usize(5));
        assert_eq!(complement.end().contig().inner(), "seq0");
        assert_eq!(complement.end().strand(), &Strand::Negative);
        assert_eq!(complement.end().position().inner(), &Value::LowerBound);

        // Testing the case where it's the maximum positive interval.
        let start = Coordinate::try_new("seq0", Strand::Positive, 0)?;
        let end = Coordinate::try_new("seq0", Strand::Positive, usize::MAX)?;

        // This represents the interval [0, usize::MAX), which should be complemented to
        // (LB, usize::MAX - 1] when it's all said and done.
        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.start().contig().inner(), "seq0");
        assert_eq!(interval.start().strand(), &Strand::Positive);
        assert_eq!(interval.start().position().inner(), &Value::Usize(0));
        assert_eq!(interval.end().contig().inner(), "seq0");
        assert_eq!(interval.end().strand(), &Strand::Positive);
        assert_eq!(interval.end().position().inner(), &Value::Usize(usize::MAX));

        let complement = interval.complement()?.unwrap();
        assert_eq!(complement.start().contig().inner(), "seq0");
        assert_eq!(complement.start().strand(), &Strand::Negative);
        assert_eq!(
            complement.start().position().inner(),
            &Value::Usize(usize::MAX - 1)
        );
        assert_eq!(complement.end().contig().inner(), "seq0");
        assert_eq!(complement.end().strand(), &Strand::Negative);
        assert_eq!(complement.end().position().inner(), &Value::LowerBound);

        Ok(())
    }

    #[test]
    fn complementing_a_negative_interval_with_two_usize_positions_works()
    -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", Strand::Negative, 5)?;
        let end = Coordinate::try_new("seq0", Strand::Negative, 0)?;

        // This represents the interval (0, 5], which should be complemented to
        // [1, 6) when it's all said and done.
        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.start().contig().inner(), "seq0");
        assert_eq!(interval.start().strand(), &Strand::Negative);
        assert_eq!(interval.start().position().inner(), &Value::Usize(5));
        assert_eq!(interval.end().contig().inner(), "seq0");
        assert_eq!(interval.end().strand(), &Strand::Negative);
        assert_eq!(interval.end().position().inner(), &Value::Usize(0));

        let complement = interval.complement()?.unwrap();
        assert_eq!(complement.start().contig().inner(), "seq0");
        assert_eq!(complement.start().strand(), &Strand::Positive);
        assert_eq!(complement.start().position().inner(), &Value::Usize(1));
        assert_eq!(complement.end().contig().inner(), "seq0");
        assert_eq!(complement.end().strand(), &Strand::Positive);
        assert_eq!(complement.end().position().inner(), &Value::Usize(6));

        // Testing converting from lower bound.
        let start = Coordinate::try_new("seq0", Strand::Negative, 5)?;
        let end = Coordinate::lower_bound("seq0");

        // This represents the interval (LB, 5], which should be complemented to
        // [0, 6) when it's all said and done.
        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.start().contig().inner(), "seq0");
        assert_eq!(interval.start().strand(), &Strand::Negative);
        assert_eq!(interval.start().position().inner(), &Value::Usize(5));
        assert_eq!(interval.end().contig().inner(), "seq0");
        assert_eq!(interval.end().strand(), &Strand::Negative);
        assert_eq!(interval.end().position().inner(), &Value::LowerBound);

        let complement = interval.complement()?.unwrap();
        assert_eq!(complement.start().contig().inner(), "seq0");
        assert_eq!(complement.start().strand(), &Strand::Positive);
        assert_eq!(complement.start().position().inner(), &Value::Usize(0));
        assert_eq!(complement.end().contig().inner(), "seq0");
        assert_eq!(complement.end().strand(), &Strand::Positive);
        assert_eq!(complement.end().position().inner(), &Value::Usize(6));

        // Testing the case where it's the maximum negative interval.
        let start = Coordinate::try_new("seq0", Strand::Negative, usize::MAX)?;
        let end = Coordinate::lower_bound("seq0");

        // This represents the interval (LB, usize::MAX], which cannot be
        // converted and should return `None`.
        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.start().contig().inner(), "seq0");
        assert_eq!(interval.start().strand(), &Strand::Negative);
        assert_eq!(
            interval.start().position().inner(),
            &Value::Usize(usize::MAX)
        );
        assert_eq!(interval.end().contig().inner(), "seq0");
        assert_eq!(interval.end().strand(), &Strand::Negative);
        assert_eq!(interval.end().position().inner(), &Value::LowerBound);

        let err = interval.complement()?;
        assert!(err.is_none());

        Ok(())
    }
}
