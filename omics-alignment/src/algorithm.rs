//! Deterministic global and local affine-gap sequence alignment.
//!
//! This module computes pairwise alignments over generic symbol slices
//! `&[T]` where `T: Eq`. Call [`crate::algorithm::global`] to align both
//! complete inputs. Call [`crate::algorithm::local`] to align the
//! highest-scoring positive subsequences and receive `Ok(None)` when either
//! input is empty or no positive local alignment exists.
//!
//! # Scoring
//!
//! [`crate::algorithm::Scoring`] stores four signed scores. A gap of length
//! `L` contributes `gap_open_score + (L - 1) * gap_extend_score`.
//!
//! # Output
//!
//! [`crate::algorithm::Outcome`] stores the alignment score, a canonical
//! [`crate::cigar::Cigar`], and half-open input ranges. These algorithms
//! emit only `=`, `X`, `I`, and `D`. `=` records equal symbols. `X` records
//! unequal symbols. `I` records a query-only gap. `D` records a
//! reference-only gap. They never emit `M`. Adjacent identical operations
//! coalesce into one run.
//!
//! The returned ranges are zero-based indexes into the original inputs.
//! [`crate::algorithm::global`] always returns
//! `0..reference.len()` and `0..query.len()`. [`crate::algorithm::local`]
//! returns the best-scoring subsequence ranges, and the returned CIGAR starts
//! at `reference_range().start` and `query_range().start` within the original
//! inputs. Callers that need genomic placement can advance the reference and
//! query coordinates for input index `0` by those starts and then pass the
//! returned CIGAR to [`crate::Alignment::try_new`].
//!
//! # Determinism
//!
//! Equal scores resolve through a fixed policy. Aligned predecessors prefer
//! aligned, then deletion, then insertion. Insertion predecessors prefer
//! insertion, then aligned, then deletion. Deletion predecessors prefer
//! deletion, then aligned, then insertion. Global endpoint selection prefers
//! aligned, then deletion, then insertion. Local endpoint selection keeps the
//! first positive maximum in row-major order over exclusive
//! `(reference_end, query_end)` pairs.
//!
//! # Example
//!
//! ```
//! use omics_alignment::algorithm::Scoring;
//! use omics_alignment::algorithm::local;
//!
//! let scoring = Scoring::try_new(2, -3, -2, -1)?;
//! let outcome = local(b"GGACGT", b"TTACGA", scoring)?.unwrap();
//!
//! assert_eq!(outcome.score(), 6);
//! assert_eq!(outcome.cigar().to_string(), "3=");
//! assert_eq!(outcome.reference_range(), &(2..5));
//! assert_eq!(outcome.query_range(), &(2..5));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

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
    /// Constructs scoring values accepted by [`global`] and [`local`].
    ///
    /// `match_score` must be positive. `mismatch_score`, `gap_open_score`,
    /// and `gap_extend_score` must be non-positive. A gap of length `L`
    /// contributes `gap_open_score + (L - 1) * gap_extend_score`.
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

    /// Returns the match or mismatch score for two symbols.
    pub(crate) fn substitution<T: Eq>(self, reference: &T, query: &T) -> Score {
        if reference == query {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The result of one deterministic non-empty pairwise alignment.
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
    ///
    /// This start index places the first CIGAR operation on the original
    /// reference input.
    pub const fn reference_range(&self) -> &Range<usize> {
        &self.reference_range
    }

    /// Returns the half-open query input range.
    ///
    /// This start index places the first CIGAR operation on the original
    /// query input.
    pub const fn query_range(&self) -> &Range<usize> {
        &self.query_range
    }

    /// Constructs an outcome from a checked traceback.
    pub(crate) fn new(
        score: Score,
        cigar: Cigar,
        reference_range: Range<usize>,
        query_range: Range<usize>,
    ) -> Self {
        Self {
            score,
            cigar,
            reference_range,
            query_range,
        }
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

/// Dynamic-programming and traceback implementation.
mod engine;

/// Computes a deterministic global affine-gap alignment over two complete
/// inputs.
///
/// This function accepts generic symbol slices `&[T]` where `T: Eq`. The
/// returned ranges always span both complete inputs. If both inputs are
/// empty, the function returns [`Error::EmptyGlobal`] because the CIGAR model
/// cannot represent a zero-operation traversal.
pub fn global<T: Eq>(reference: &[T], query: &[T], scoring: Scoring) -> Result<Outcome, Error> {
    engine::global(reference, query, scoring)
}

/// Computes the highest-scoring positive local affine-gap alignment.
///
/// This function accepts generic symbol slices `&[T]` where `T: Eq`.
/// `Ok(None)` means that either input is empty or no positive-scoring local
/// alignment exists. Successful results report half-open subsequence ranges in
/// the original inputs, and their CIGAR begins and ends with `=` or `X`.
pub fn local<T: Eq>(
    reference: &[T],
    query: &[T],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    engine::local(reference, query, scoring)
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

    #[test]
    fn global_aligns_matches_and_one_deletion() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, -3, -2, -1)?;
        let outcome = global(b"ACGT", b"AGT", scoring)?;
        assert_eq!(outcome.score(), 4);
        assert_eq!(outcome.cigar().to_string(), "1=1D2=");
        assert_eq!(outcome.reference_range(), &(0..4));
        assert_eq!(outcome.query_range(), &(0..3));
        Ok(())
    }

    #[test]
    fn global_prefers_two_gaps_to_an_expensive_mismatch() -> Result<(), Error> {
        let scoring = Scoring::try_new(1, -5, -1, -1)?;
        let outcome = global(b"A", b"B", scoring)?;
        assert_eq!(outcome.score(), -2);
        assert_eq!(outcome.cigar().to_string(), "1I1D");
        Ok(())
    }

    #[test]
    fn global_handles_one_empty_input() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, -3, -2, -1)?;
        let outcome = global(b"", b"AAA", scoring)?;
        assert_eq!(outcome.score(), -4);
        assert_eq!(outcome.cigar().to_string(), "3I");
        Ok(())
    }

    #[test]
    fn global_rejects_two_empty_inputs() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, -3, -2, -1)?;
        assert!(matches!(
            global::<u8>(b"", b"", scoring),
            Err(Error::EmptyGlobal)
        ));
        Ok(())
    }

    #[test]
    fn global_reports_intermediate_score_overflow() -> Result<(), Error> {
        let scoring = Scoring::try_new(1, -1, Score::MIN, -1)?;
        assert!(matches!(
            global(b"", b"AA", scoring),
            Err(Error::ScoreOverflow)
        ));
        Ok(())
    }

    #[test]
    fn local_returns_best_positive_subsequences() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, -3, -2, -1)?;
        let outcome = local(b"GGACGTCC", b"TTACGAAA", scoring)?.unwrap();
        assert_eq!(outcome.score(), 6);
        assert_eq!(outcome.cigar().to_string(), "3=");
        assert_eq!(outcome.reference_range(), &(2..5));
        assert_eq!(outcome.query_range(), &(2..5));
        Ok(())
    }

    #[test]
    fn local_returns_none_without_a_positive_alignment() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, -3, -2, -1)?;
        assert!(local(b"AAAA", b"CCCC", scoring)?.is_none());
        assert!(local::<u8>(b"", b"", scoring)?.is_none());
        assert!(local(b"", b"AAAA", scoring)?.is_none());
        assert!(local(b"AAAA", b"", scoring)?.is_none());
        Ok(())
    }

    #[test]
    fn local_trims_zero_cost_leading_gaps() -> Result<(), Error> {
        let scoring = Scoring::try_new(1, 0, 0, 0)?;
        let outcome = local(b"A", b"BA", scoring)?.unwrap();
        assert_eq!(outcome.score(), 1);
        assert_eq!(outcome.cigar().to_string(), "1=");
        assert_eq!(outcome.reference_range(), &(0..1));
        assert_eq!(outcome.query_range(), &(1..2));
        Ok(())
    }

    #[test]
    fn local_chooses_the_earliest_equal_endpoint() -> Result<(), Error> {
        let scoring = Scoring::try_new(1, 0, 0, 0)?;
        let outcome = local(b"AA", b"A", scoring)?.unwrap();
        assert_eq!(outcome.cigar().to_string(), "1=");
        assert_eq!(outcome.reference_range(), &(0..1));
        assert_eq!(outcome.query_range(), &(0..1));
        Ok(())
    }
}
