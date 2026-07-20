//! A single lossless step of an [`Alignment`](crate::Alignment) traversal.

use omics_coordinate::interval::interbase::Interval;

use crate::cigar::Operation;

/// A single step of an alignment, corresponding to exactly one CIGAR
/// operation.
///
/// Every operation is preserved, including hard clips (`H`) and padding (`P`)
/// which consume neither the reference nor the query axis. When an axis is
/// not consumed by the operation, the corresponding interval is `None`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Step {
    /// The CIGAR operation this step corresponds to.
    operation: Operation,
    /// The reference interbase interval consumed by this step, if any.
    reference: Option<Interval>,
    /// The query interbase interval consumed by this step, if any.
    query: Option<Interval>,
}

impl Step {
    /// Constructs a new [`Step`] from an operation and its optional reference
    /// and query intervals.
    pub(crate) fn new(
        operation: Operation,
        reference: Option<Interval>,
        query: Option<Interval>,
    ) -> Self {
        Self {
            operation,
            reference,
            query,
        }
    }

    /// Returns the CIGAR operation this step corresponds to.
    pub const fn operation(&self) -> Operation {
        self.operation
    }

    /// Returns the reference interval consumed by this step.
    ///
    /// This is `None` when the operation does not consume the reference
    /// axis, such as an insertion, a hard clip, or padding.
    pub fn reference(&self) -> Option<&Interval> {
        self.reference.as_ref()
    }

    /// Returns the query interval consumed by this step.
    ///
    /// This is `None` when the operation does not consume the query axis,
    /// such as a deletion, a reference skip, or padding.
    pub fn query(&self) -> Option<&Interval> {
        self.query.as_ref()
    }

    /// Returns `true` if this step contributes to the pairwise base
    /// alignment (its operation kind is `M`, `=`, or `X`).
    pub const fn is_aligned(&self) -> bool {
        self.operation.kind().is_aligned()
    }
}
