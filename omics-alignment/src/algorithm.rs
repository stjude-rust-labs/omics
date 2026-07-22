//! Pairwise global and local sequence-alignment algorithms.

use std::collections::TryReserveError;
use std::ops::Range;

use thiserror::Error;

use crate::cigar::Axis;
use crate::cigar::Cigar;
use crate::cigar::OperationError;

/// A score used by pairwise alignment algorithms.
pub type Score = i64;

/// Validated scores for pairwise alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Scoring {
    /// Score for equal symbols.
    match_score: Score,
    /// Score for unequal symbols.
    mismatch_score: Score,
    /// Score for the first symbol in a gap.
    gap_open_score: Score,
    /// Score for each additional symbol in a gap.
    gap_extend_score: Score,
}

impl Scoring {
    /// Constructs scoring values accepted by global and local alignment.
    pub fn try_new(
        match_score: Score,
        mismatch_score: Score,
        gap_open_score: Score,
        gap_extend_score: Score,
    ) -> Result<Self, Error> {
        if match_score <= 0 {
            return Err(Error::InvalidMatchScore { score: match_score });
        }
        if mismatch_score > 0 {
            return Err(Error::InvalidMismatchScore {
                score: mismatch_score,
            });
        }
        if gap_open_score > 0 {
            return Err(Error::InvalidGapOpenScore {
                score: gap_open_score,
            });
        }
        if gap_extend_score > 0 {
            return Err(Error::InvalidGapExtendScore {
                score: gap_extend_score,
            });
        }

        Ok(Self {
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        })
    }

    /// Returns the score for equal symbols.
    pub const fn match_score(self) -> Score {
        self.match_score
    }

    /// Returns the score for unequal symbols.
    pub const fn mismatch_score(self) -> Score {
        self.mismatch_score
    }

    /// Returns the score for the first symbol in a gap.
    pub const fn gap_open_score(self) -> Score {
        self.gap_open_score
    }

    /// Returns the score for each additional symbol in a gap.
    pub const fn gap_extend_score(self) -> Score {
        self.gap_extend_score
    }
}

/// The result of one non-empty pairwise alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Outcome {
    /// Alignment score.
    score: Score,
    /// Canonical alignment operations.
    cigar: Cigar,
    /// Covered half-open reference input range.
    reference_range: Range<usize>,
    /// Covered half-open query input range.
    query_range: Range<usize>,
}

impl Outcome {
    /// Returns the alignment score.
    pub const fn score(&self) -> Score {
        self.score
    }

    /// Returns the canonical CIGAR.
    pub const fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns the half-open reference input range.
    pub const fn reference_range(&self) -> &Range<usize> {
        &self.reference_range
    }

    /// Returns the half-open query input range.
    pub const fn query_range(&self) -> &Range<usize> {
        &self.query_range
    }
}

/// An error produced while configuring or computing an alignment.
#[derive(Debug, Error)]
pub enum Error {
    /// The match score was not positive.
    #[error("match score must be positive, got {score}")]
    InvalidMatchScore {
        /// Invalid score supplied by the caller.
        score: Score,
    },
    /// The mismatch score was positive.
    #[error("mismatch score must be non-positive, got {score}")]
    InvalidMismatchScore {
        /// Invalid score supplied by the caller.
        score: Score,
    },
    /// The gap-open score was positive.
    #[error("gap-open score must be non-positive, got {score}")]
    InvalidGapOpenScore {
        /// Invalid score supplied by the caller.
        score: Score,
    },
    /// The gap-extension score was positive.
    #[error("gap-extension score must be non-positive, got {score}")]
    InvalidGapExtendScore {
        /// Invalid score supplied by the caller.
        score: Score,
    },
    /// Both inputs to global alignment were empty.
    #[error("a global alignment of two empty inputs has no representable CIGAR")]
    EmptyGlobal,
    /// An input length did not fit the active coordinate number type.
    #[error("{axis} input length {length} does not fit the configured CIGAR length type")]
    LengthOutOfRange {
        /// Input axis whose length is not representable.
        axis: Axis,
        /// Unrepresentable slice length.
        length: usize,
    },
    /// Matrix row and column multiplication overflowed.
    #[error("alignment matrix dimensions overflow usize")]
    MatrixSizeOverflow,
    /// Matrix memory could not be reserved.
    #[error("failed to allocate alignment matrix")]
    MatrixAllocation {
        /// Allocation failure reported by `Vec`.
        #[source]
        source: TryReserveError,
    },
    /// Checked score arithmetic overflowed.
    #[error("alignment score arithmetic overflowed i64")]
    ScoreOverflow,
    /// A reachable endpoint produced an empty or incomplete traceback.
    #[error("alignment traceback invariant failed")]
    TracebackInvariant,
    /// A traceback operation could not be constructed.
    #[error("failed to construct a traceback operation")]
    Operation {
        /// Checked operation-construction failure.
        #[from]
        source: OperationError,
    },
    /// The checked operation list could not form a CIGAR.
    #[error("failed to construct the alignment CIGAR")]
    Cigar {
        /// Checked CIGAR-construction failure.
        #[from]
        source: crate::cigar::Error,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scoring_accepts_valid_boundaries() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, 0, 0, -1)?;
        assert_eq!(scoring.match_score(), 2);
        assert_eq!(scoring.mismatch_score(), 0);
        assert_eq!(scoring.gap_open_score(), 0);
        assert_eq!(scoring.gap_extend_score(), -1);
        Ok(())
    }

    #[test]
    fn scoring_rejects_invalid_signs() {
        assert!(matches!(
            Scoring::try_new(0, -1, -2, -1),
            Err(Error::InvalidMatchScore { score: 0 })
        ));
        assert!(matches!(
            Scoring::try_new(1, 1, -2, -1),
            Err(Error::InvalidMismatchScore { score: 1 })
        ));
        assert!(matches!(
            Scoring::try_new(1, -1, 1, -1),
            Err(Error::InvalidGapOpenScore { score: 1 })
        ));
        assert!(matches!(
            Scoring::try_new(1, -1, -2, 1),
            Err(Error::InvalidGapExtendScore { score: 1 })
        ));
    }

    #[test]
    fn outcome_exposes_immutable_result_parts() -> Result<(), Box<dyn std::error::Error>> {
        let outcome = Outcome {
            score: 7,
            cigar: "2=1X".parse()?,
            reference_range: 3..6,
            query_range: 5..8,
        };
        assert_eq!(outcome.score(), 7);
        assert_eq!(outcome.cigar().to_string(), "2=1X");
        assert_eq!(outcome.reference_range(), &(3..6));
        assert_eq!(outcome.query_range(), &(5..8));
        Ok(())
    }
}
