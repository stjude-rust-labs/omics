//! Narrow score lanes and the per-state containers built from them.
//!
//! [`LaneScore`] abstracts over the `i16` and `i32` lane types and reserves
//! each type's minimum as the unreachable sentinel. `ByState`,
//! [`ByNeighbour`], and [`ByDelta`] are small fixed-field containers that
//! let the engine index by recurrence vocabulary, for example
//! `neighbours.get(state.neighbour())`, instead of repeating three-way
//! matches.

use crate::algorithm::Error;
use crate::algorithm::Score;
use crate::algorithm::Scoring;
use crate::algorithm::recurrence::Delta;
use crate::algorithm::recurrence::Neighbour;
use crate::algorithm::recurrence::State;

/// A score lane that reserves its minimum value for unreachable states.
pub(in crate::algorithm::simd) trait LaneScore: Copy + Eq + Ord {
    /// The sentinel used for an unreachable score.
    const UNREACHABLE: Self;
    /// The greatest score that can remain reachable in this lane.
    const MAX: Self;
    /// The score used for a reachable origin or local reset.
    const ZERO: Self;

    /// Converts a scalar score when it fits without using the sentinel.
    fn try_from_score(score: Score) -> Option<Self>;

    /// Converts this lane value to the scalar score representation.
    fn into_score(self) -> Score;

    /// Converts a packed predecessor byte into a lane value.
    fn from_code(code: u8) -> Self;

    /// Returns the low byte of a lane that holds a packed predecessor byte.
    fn into_code(self) -> u8;
}

impl LaneScore for i16 {
    const MAX: Self = i16::MAX;
    const UNREACHABLE: Self = i16::MIN;
    const ZERO: Self = 0;

    #[inline]
    fn try_from_score(score: Score) -> Option<Self> {
        let score = Self::try_from(score).ok()?;
        (score != i16::MIN).then_some(score)
    }

    #[inline]
    fn into_score(self) -> Score {
        Score::from(self)
    }

    #[inline]
    fn from_code(code: u8) -> Self {
        Self::from(code)
    }

    #[inline]
    fn into_code(self) -> u8 {
        self as u8
    }
}

impl LaneScore for i32 {
    const MAX: Self = i32::MAX;
    const UNREACHABLE: Self = i32::MIN;
    const ZERO: Self = 0;

    #[inline]
    fn try_from_score(score: Score) -> Option<Self> {
        let score = Self::try_from(score).ok()?;
        (score != i32::MIN).then_some(score)
    }

    #[inline]
    fn into_score(self) -> Score {
        Score::from(self)
    }

    #[inline]
    fn from_code(code: u8) -> Self {
        Self::from(code)
    }

    #[inline]
    fn into_code(self) -> u8 {
        self as u8
    }
}

/// One value for each recurrence state.
#[derive(Clone, Copy, Debug)]
pub(super) struct ByState<T> {
    /// The aligned state's value.
    pub(super) aligned: T,
    /// The insertion state's value.
    pub(super) insertion: T,
    /// The deletion state's value.
    pub(super) deletion: T,
}

impl<T> ByState<T> {
    /// Returns the value for `state`.
    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    pub(super) fn get(&self, state: State) -> &T {
        match state {
            State::Aligned => &self.aligned,
            State::Insertion => &self.insertion,
            State::Deletion => &self.deletion,
        }
    }
}

/// One value for each neighbouring cell.
pub(super) struct ByNeighbour<T> {
    /// The diagonal neighbour's value.
    pub(super) diagonal: T,
    /// The left neighbour's value.
    pub(super) left: T,
    /// The up neighbour's value.
    pub(super) up: T,
}

impl<T> ByNeighbour<T> {
    /// Returns the value for `neighbour`.
    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    pub(super) fn get(&self, neighbour: Neighbour) -> &T {
        match neighbour {
            Neighbour::Diagonal => &self.diagonal,
            Neighbour::Left => &self.left,
            Neighbour::Up => &self.up,
        }
    }
}

/// One value for each kind of recurrence delta.
pub(super) struct ByDelta<T> {
    /// The substitution score.
    pub(super) substitution: T,
    /// The gap-open score.
    pub(super) gap_open: T,
    /// The gap-extension score.
    pub(super) gap_extend: T,
}

impl<T> ByDelta<T> {
    /// Returns the value for `delta`.
    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    pub(super) fn get(&self, delta: Delta) -> &T {
        match delta {
            Delta::Substitution => &self.substitution,
            Delta::GapOpen => &self.gap_open,
            Delta::GapExtend => &self.gap_extend,
        }
    }
}

/// All three state scores for one retained anti-diagonal.
pub(super) type DiagonalScores<S> = ByState<Vec<S>>;

impl<S: LaneScore> DiagonalScores<S> {
    /// Allocates unreachable scores for one maximum-width anti-diagonal plus
    /// the slack read and written by a final partial vector.
    #[inline]
    pub(super) fn try_new(length: usize) -> Result<Self, Error> {
        Ok(Self {
            aligned: allocated_scores(length)?,
            insertion: allocated_scores(length)?,
            deletion: allocated_scores(length)?,
        })
    }
}

/// Allocates one unreachable score array.
#[inline]
pub(super) fn allocated_scores<S: LaneScore>(length: usize) -> Result<Vec<S>, Error> {
    let mut scores = Vec::new();
    scores
        .try_reserve_exact(length)
        .map_err(|source| Error::MatrixAllocation { source })?;
    scores.resize(length, S::UNREACHABLE);
    Ok(scores)
}

/// The scoring parameters converted to one lane's score type.
pub(super) struct LaneScoring<S> {
    /// Score for an equal byte pair.
    pub(super) match_score: S,
    /// Score for an unequal byte pair.
    pub(super) mismatch_score: S,
    /// Score for opening a gap.
    pub(super) gap_open_score: S,
    /// Score for extending a gap.
    pub(super) gap_extend_score: S,
}

impl<S: LaneScore> LaneScoring<S> {
    /// Converts all scalar scoring values into narrow lanes.
    #[inline]
    pub(super) fn try_new(scoring: Scoring) -> Result<Self, Error> {
        Ok(Self {
            match_score: narrow(scoring.match_score())?,
            mismatch_score: narrow(scoring.mismatch_score())?,
            gap_open_score: narrow(scoring.gap_open_score())?,
            gap_extend_score: narrow(scoring.gap_extend_score())?,
        })
    }
}

/// Converts one scalar score into a usable narrow lane.
#[inline]
pub(super) fn narrow<S: LaneScore>(score: Score) -> Result<S, Error> {
    S::try_from_score(score).ok_or(Error::ScoreOverflow)
}
