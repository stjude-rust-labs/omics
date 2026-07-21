//! Validated pairwise alignments with lossless per-operation traversal.
//!
//! # Quick start
//!
//! ```
//! use omics_alignment::Alignment;
//! use omics_alignment::cigar::Cigar;
//! use omics_coordinate::interbase::Coordinate;
//!
//! let reference_start = "ref:+:0".parse::<Coordinate>()?;
//! let query_start = "query:+:0".parse::<Coordinate>()?;
//! let cigar = "3M1I1M".parse::<Cigar>()?;
//!
//! let alignment = Alignment::try_new(reference_start, query_start, cigar)?;
//! assert_eq!(alignment.reference_end().to_string(), "ref:+:4");
//! assert_eq!(alignment.query_end().to_string(), "query:+:5");
//! assert_eq!(alignment.steps().count(), 3);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use omics_coordinate::interbase::Coordinate;
use omics_coordinate::interval::interbase::Interval;
use omics_coordinate::position::Number;
use thiserror::Error;

use crate::cigar::Axis;
use crate::cigar::Cigar;
use crate::cigar::Operation;
use crate::step::Step;

/// An error constructing an [`Alignment`].
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum Error {
    /// An operation moved a coordinate out of representable bounds.
    #[error(
        "operation {operation_index} ({operation}) moved the {axis} coordinate out of bounds at \
         {coordinate}"
    )]
    OutOfBounds {
        /// Zero-based index of the offending operation within the CIGAR.
        operation_index: usize,
        /// The offending operation.
        operation: Operation,
        /// The axis on which the movement failed.
        axis: Axis,
        /// The coordinate at which the movement failed. This is the
        /// coordinate as it stood before the failed move, since a failed
        /// move never mutates its input.
        coordinate: Coordinate,
    },
}

/// A validated pairwise alignment between a reference sequence and a query
/// sequence.
///
/// An [`Alignment`] pairs a starting reference coordinate, a starting query
/// coordinate, and a [`Cigar`] describing the operations that transform the
/// query into alignment with the reference. Construction eagerly validates
/// that every operation's movement stays within representable coordinate
/// bounds, so [`Alignment::steps`] is infallible.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alignment {
    /// The starting reference coordinate.
    reference_start: Coordinate,
    /// The starting query coordinate.
    query_start: Coordinate,
    /// The validated CIGAR describing the alignment operations.
    cigar: Cigar,
    /// The ending reference coordinate, reached after consuming every
    /// reference-consuming operation in `cigar`.
    reference_end: Coordinate,
    /// The ending query coordinate, reached after consuming every
    /// query-consuming operation in `cigar`.
    query_end: Coordinate,
}

impl Alignment {
    /// Constructs a new [`Alignment`], eagerly validating that every
    /// operation in `cigar` can be applied to `reference_start` and
    /// `query_start` without moving a coordinate out of bounds.
    ///
    /// For each operation, the reference axis is checked (and, if the
    /// operation consumes it, moved) before the query axis is touched. If
    /// the reference move fails, the query pointer for that operation is
    /// never mutated, so the [`Error::OutOfBounds`] reported always
    /// reflects the first axis that failed.
    ///
    /// # Errors
    ///
    /// Returns [`Error::OutOfBounds`] if any operation would move the
    /// reference or query coordinate out of representable bounds.
    pub fn try_new(
        reference_start: Coordinate,
        query_start: Coordinate,
        cigar: Cigar,
    ) -> Result<Self, Error> {
        let mut reference = reference_start.clone();
        let mut query = query_start.clone();

        for (operation_index, operation) in cigar.iter().enumerate() {
            let kind = operation.kind();
            let length = operation.length();

            if kind.consumes_reference() && !reference.move_forward(length) {
                return Err(Error::OutOfBounds {
                    operation_index,
                    operation: *operation,
                    axis: Axis::Reference,
                    coordinate: reference,
                });
            }

            if kind.consumes_query() && !query.move_forward(length) {
                return Err(Error::OutOfBounds {
                    operation_index,
                    operation: *operation,
                    axis: Axis::Query,
                    coordinate: query,
                });
            }
        }

        Ok(Self {
            reference_start,
            query_start,
            cigar,
            reference_end: reference,
            query_end: query,
        })
    }

    /// Returns the [`Cigar`] describing this alignment's operations.
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns the starting reference coordinate.
    pub fn reference_start(&self) -> &Coordinate {
        &self.reference_start
    }

    /// Returns the starting query coordinate.
    pub fn query_start(&self) -> &Coordinate {
        &self.query_start
    }

    /// Returns the ending reference coordinate, reached after consuming
    /// every reference-consuming operation.
    pub fn reference_end(&self) -> &Coordinate {
        &self.reference_end
    }

    /// Returns the ending query coordinate, reached after consuming every
    /// query-consuming operation.
    pub fn query_end(&self) -> &Coordinate {
        &self.query_end
    }

    /// Returns the total number of reference bases consumed by this
    /// alignment's [`Cigar`].
    pub fn reference_length(&self) -> Number {
        self.cigar.reference_length()
    }

    /// Returns the total number of query bases consumed by this alignment's
    /// [`Cigar`].
    pub fn query_length(&self) -> Number {
        self.cigar.query_length()
    }

    /// Returns an iterator yielding one [`Step`] per operation in this
    /// alignment's [`Cigar`], in order.
    ///
    /// Every operation is preserved, including hard clips and padding, which
    /// yield a [`Step`] with no reference or query interval. No item in this
    /// iterator is fallible; [`Alignment::try_new`] has already proven that
    /// every consumed move stays within bounds.
    pub fn steps(&self) -> impl Iterator<Item = Step> + '_ {
        let mut reference = self.reference_start.clone();
        let mut query = self.query_start.clone();

        self.cigar.iter().copied().map(move |operation| {
            let kind = operation.kind();
            let length = operation.length();

            let reference_interval = kind
                .consumes_reference()
                .then(|| Self::advance(&mut reference, length));
            let query_interval = kind
                .consumes_query()
                .then(|| Self::advance(&mut query, length));

            Step::new(operation, reference_interval, query_interval)
        })
    }

    /// Moves `pointer` forward by `length` and returns the interbase
    /// interval spanning its position before and after the move.
    ///
    /// This is only ever called for an axis that the current operation
    /// consumes, so `length` is positive by construction of [`Operation`].
    ///
    /// [`Coordinate::move_forward`](omics_coordinate::Coordinate::move_forward)
    /// never changes a coordinate's contig or strand; it advances the
    /// position on the positive strand and retreats it on the negative
    /// strand. Either way, the coordinate observed before the move and the
    /// coordinate observed after the move are ordered exactly as
    /// [`Interval::try_new`] requires of a start and an end on that strand
    /// (start at or before end on the positive strand, end at or before
    /// start on the negative strand).
    ///
    /// [`Alignment::try_new`] has already replayed this identical sequence of
    /// operations from these identical starting coordinates and proven that
    /// every consumed move stays within representable bounds, so the move
    /// performed here cannot overflow or underflow and the interval
    /// constructed from it cannot be rejected.
    fn advance(pointer: &mut Coordinate, length: Number) -> Interval {
        let start = pointer.clone();
        // The move below cannot fail; see the proof on this function.
        let moved = pointer.move_forward(length);
        debug_assert!(moved, "movement already validated by Alignment::try_new");

        // SAFETY: the proof on this function guarantees `start` and
        // `pointer` are ordered as `Interval::try_new` requires for this
        // strand.
        Interval::try_new(start, pointer.clone()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::interbase::Coordinate;
    use omics_coordinate::position::Number;

    use super::*;
    use crate::cigar::Axis;
    use crate::cigar::Cigar;
    use crate::cigar::OperationKind;

    #[test]
    fn placement_and_steps() -> Result<(), Box<dyn std::error::Error>> {
        let reference_start = "ref:+:0".parse::<Coordinate>()?;
        let query_start = "query:-:5".parse::<Coordinate>()?;
        let cigar = "3M1I1M".parse::<Cigar>()?;

        let alignment = Alignment::try_new(reference_start, query_start, cigar)?;

        let steps = alignment.steps().collect::<Vec<_>>();
        assert_eq!(steps.len(), 3);
        assert_eq!(steps[0].reference().unwrap().to_string(), "ref:+:0-3");
        assert_eq!(steps[0].query().unwrap().to_string(), "query:-:5-2");
        assert!(steps[1].reference().is_none());
        assert_eq!(steps[1].query().unwrap().to_string(), "query:-:2-1");
        assert_eq!(steps[2].reference().unwrap().to_string(), "ref:+:3-4");
        assert_eq!(steps[2].query().unwrap().to_string(), "query:-:1-0");
        assert_eq!(alignment.reference_end().to_string(), "ref:+:4");
        assert_eq!(alignment.query_end().to_string(), "query:-:0");

        Ok(())
    }

    #[test]
    fn every_operation_kind_yields_one_step() -> Result<(), Box<dyn std::error::Error>> {
        // Clipping-valid nine-kind CIGAR (H and S at the terminal positions
        // allowed by the SAM specification).
        let reference_start = "ref:+:100".parse::<Coordinate>()?;
        let query_start = "query:+:200".parse::<Coordinate>()?;
        let cigar = "1H1S1M1I1=1X1D1N1P1S1H".parse::<Cigar>()?;

        let alignment = Alignment::try_new(reference_start, query_start, cigar)?;
        let steps = alignment.steps().collect::<Vec<_>>();

        let expected: [(OperationKind, Option<&str>, Option<&str>); 11] = [
            (OperationKind::HardClip, None, None),
            (OperationKind::SoftClip, None, Some("query:+:200-201")),
            (
                OperationKind::Match,
                Some("ref:+:100-101"),
                Some("query:+:201-202"),
            ),
            (OperationKind::Insertion, None, Some("query:+:202-203")),
            (
                OperationKind::SequenceMatch,
                Some("ref:+:101-102"),
                Some("query:+:203-204"),
            ),
            (
                OperationKind::SequenceMismatch,
                Some("ref:+:102-103"),
                Some("query:+:204-205"),
            ),
            (OperationKind::Deletion, Some("ref:+:103-104"), None),
            (OperationKind::ReferenceSkip, Some("ref:+:104-105"), None),
            (OperationKind::Padding, None, None),
            (OperationKind::SoftClip, None, Some("query:+:205-206")),
            (OperationKind::HardClip, None, None),
        ];

        assert_eq!(steps.len(), expected.len());
        for (step, (kind, reference, query)) in steps.iter().zip(expected.iter()) {
            assert_eq!(step.operation().kind(), *kind);
            assert_eq!(
                step.reference().map(|i| i.to_string()),
                reference.map(String::from)
            );
            assert_eq!(step.query().map(|i| i.to_string()), query.map(String::from));
        }

        assert_eq!(alignment.reference_end().to_string(), "ref:+:105");
        assert_eq!(alignment.query_end().to_string(), "query:+:206");

        Ok(())
    }

    #[test]
    fn positive_strand_reference_overflow_is_rejected() -> Result<(), Box<dyn std::error::Error>> {
        let reference_start = Coordinate::try_new("ref", "+", Number::MAX)?;
        let query_start = "query:+:0".parse::<Coordinate>()?;
        let cigar = "1M".parse::<Cigar>()?;

        let error = Alignment::try_new(reference_start, query_start, cigar).unwrap_err();
        assert!(matches!(
            error,
            Error::OutOfBounds {
                operation_index: 0,
                axis: Axis::Reference,
                ..
            }
        ));

        Ok(())
    }

    #[test]
    fn negative_strand_query_underflow_is_rejected() -> Result<(), Box<dyn std::error::Error>> {
        let reference_start = "ref:+:0".parse::<Coordinate>()?;
        let query_start = Coordinate::try_new("query", "-", 0)?;
        let cigar = "1I".parse::<Cigar>()?;

        let error = Alignment::try_new(reference_start, query_start, cigar).unwrap_err();
        assert!(matches!(
            error,
            Error::OutOfBounds {
                operation_index: 0,
                axis: Axis::Query,
                ..
            }
        ));

        Ok(())
    }

    #[test]
    fn aligned_steps_preserve_operation_boundaries() -> Result<(), Box<dyn std::error::Error>> {
        let reference_start = "ref:+:0".parse::<Coordinate>()?;
        let query_start = "query:+:0".parse::<Coordinate>()?;
        let cigar = "2M1=3X1I4M1D5M1P6M".parse::<Cigar>()?;

        let alignment = Alignment::try_new(reference_start, query_start, cigar)?;
        assert_eq!(
            alignment
                .steps()
                .filter(Step::is_aligned)
                .map(|step| step.operation().length())
                .collect::<Vec<_>>(),
            vec![2, 1, 3, 4, 5, 6]
        );

        Ok(())
    }

    #[test]
    fn aligned_steps_match_chainfile_geometry() -> Result<(), Box<dyn std::error::Error>> {
        let reference_start = "ref:+:0".parse::<Coordinate>()?;
        let query_start = "query:-:5".parse::<Coordinate>()?;
        let cigar = "3M1I1M".parse::<Cigar>()?;

        let alignment = Alignment::try_new(reference_start, query_start, cigar)?;
        let steps = alignment
            .steps()
            .filter(Step::is_aligned)
            .map(|step| {
                (
                    step.reference().map(ToString::to_string),
                    step.query().map(ToString::to_string),
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(
            steps,
            vec![
                (
                    Some("ref:+:0-3".to_string()),
                    Some("query:-:5-2".to_string())
                ),
                (
                    Some("ref:+:3-4".to_string()),
                    Some("query:-:1-0".to_string())
                ),
            ]
        );

        Ok(())
    }
}
