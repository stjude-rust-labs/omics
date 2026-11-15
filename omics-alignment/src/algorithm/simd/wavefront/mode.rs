//! The global and local alignment modes plugged into the shared sweep.
//!
//! [`AlignmentMode`] isolates everything that differs between the modes: the
//! row-zero and column-zero boundary cells, the local clamp of aligned scores
//! at zero, and endpoint selection. [`Global`] ends at the bottom-right cell;
//! [`Local`] keeps the first positive maximum in row-major order. The scalar
//! helpers here compute boundary cells with the same tie-break order as the
//! vector path.

use std::cmp::Reverse;

use super::kernel::Kernel;
use super::lane::DiagonalScores;
use super::lane::LaneScore;
use super::lane::LaneScoring;
use super::lane::narrow;
use super::layout::Diagonal;
use super::layout::Layout;
use super::traceback::pack_predecessors;
use crate::algorithm::Error;
use crate::algorithm::Score;
use crate::algorithm::recurrence::Delta;
use crate::algorithm::recurrence::State;

/// Records the terminal state and coordinate selected by an alignment mode.
#[derive(Clone, Copy)]
pub(super) struct Endpoint {
    /// Terminal alignment score.
    pub(super) score: Score,
    /// Terminal recurrence state.
    pub(super) state: State,
    /// Exclusive reference coordinate.
    pub(super) row: usize,
    /// Exclusive query coordinate.
    pub(super) column: usize,
}

/// The global or local behaviour plugged into the shared wavefront sweep.
pub(super) trait AlignmentMode {
    /// Writes the boundary cells (row zero and column zero) of one diagonal.
    fn write_boundaries<S: LaneScore>(
        diagonal_index: usize,
        layout: &Layout,
        current: &mut DiagonalScores<S>,
        previous: &DiagonalScores<S>,
        lane_scoring: &LaneScoring<S>,
        traceback: &mut [u8],
    ) -> Result<(), Error>;

    /// Applies the mode's aligned-state rule to one vector.
    fn clamp_aligned<K: Kernel>(
        kernel: K,
        zero: K::Vector,
        aligned: K::Vector,
        predecessor: K::Vector,
    ) -> (K::Vector, K::Vector);

    /// Updates the endpoint after one diagonal's interior is computed.
    fn track_endpoint<S: LaneScore>(
        endpoint: &mut Option<Endpoint>,
        diagonal_index: usize,
        diagonal: Diagonal,
        current: &DiagonalScores<S>,
        interior: (usize, usize),
    );

    /// Returns the endpoint after the final diagonal is computed.
    fn finish<S: LaneScore>(
        endpoint: Option<Endpoint>,
        layout: &Layout,
        last: &DiagonalScores<S>,
    ) -> Option<Endpoint>;
}

/// Aligns both complete input slices.
pub(super) struct Global;

impl AlignmentMode for Global {
    #[inline]
    fn write_boundaries<S: LaneScore>(
        diagonal_index: usize,
        layout: &Layout,
        current: &mut DiagonalScores<S>,
        previous: &DiagonalScores<S>,
        lane_scoring: &LaneScoring<S>,
        traceback: &mut [u8],
    ) -> Result<(), Error> {
        let diagonal = layout.diagonals[diagonal_index];

        if diagonal_index == 0 {
            set_cell(
                current,
                0,
                (S::ZERO, None),
                (S::UNREACHABLE, None),
                (S::UNREACHABLE, None),
                traceback,
                diagonal.offset,
            );
            return Ok(());
        }

        let previous_diagonal = layout.diagonals[diagonal_index - 1];

        if diagonal_index < layout.columns {
            let (insertion, predecessor) =
                choose_boundary(State::Insertion, previous, 0, lane_scoring)?;
            set_cell(
                current,
                0,
                (S::UNREACHABLE, None),
                (insertion, predecessor),
                (S::UNREACHABLE, None),
                traceback,
                diagonal.offset,
            );
        }

        if diagonal_index < layout.rows {
            let score_index = diagonal_index - diagonal.first_row;
            let previous_index = diagonal_index - 1 - previous_diagonal.first_row;
            let (deletion, predecessor) =
                choose_boundary(State::Deletion, previous, previous_index, lane_scoring)?;
            set_cell(
                current,
                score_index,
                (S::UNREACHABLE, None),
                (S::UNREACHABLE, None),
                (deletion, predecessor),
                traceback,
                diagonal.offset + score_index,
            );
        }

        Ok(())
    }

    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    fn clamp_aligned<K: Kernel>(
        _kernel: K,
        _zero: K::Vector,
        aligned: K::Vector,
        predecessor: K::Vector,
    ) -> (K::Vector, K::Vector) {
        (aligned, predecessor)
    }

    #[inline]
    fn track_endpoint<S: LaneScore>(
        _endpoint: &mut Option<Endpoint>,
        _diagonal_index: usize,
        _diagonal: Diagonal,
        _current: &DiagonalScores<S>,
        _interior: (usize, usize),
    ) {
    }

    #[inline]
    fn finish<S: LaneScore>(
        _endpoint: Option<Endpoint>,
        layout: &Layout,
        last: &DiagonalScores<S>,
    ) -> Option<Endpoint> {
        let row = layout.rows - 1;
        let column = layout.columns - 1;
        let score_index = row - layout.diagonals[row + column].first_row;
        let (score, state) = select_score(
            State::GLOBAL_ENDPOINT_ORDER.map(|state| (last.get(state)[score_index], state)),
        );

        state.map(|state| Endpoint {
            score: score.into_score(),
            state,
            row,
            column,
        })
    }
}

/// Aligns the highest-scoring positive subsequences.
pub(super) struct Local;

impl AlignmentMode for Local {
    #[inline]
    fn write_boundaries<S: LaneScore>(
        diagonal_index: usize,
        layout: &Layout,
        current: &mut DiagonalScores<S>,
        _previous: &DiagonalScores<S>,
        _values: &LaneScoring<S>,
        traceback: &mut [u8],
    ) -> Result<(), Error> {
        let diagonal = layout.diagonals[diagonal_index];

        if diagonal_index < layout.columns {
            set_cell(
                current,
                0,
                (S::ZERO, None),
                (S::UNREACHABLE, None),
                (S::UNREACHABLE, None),
                traceback,
                diagonal.offset,
            );
        }

        if diagonal_index != 0 && diagonal_index < layout.rows {
            let score_index = diagonal_index - diagonal.first_row;
            set_cell(
                current,
                score_index,
                (S::ZERO, None),
                (S::UNREACHABLE, None),
                (S::UNREACHABLE, None),
                traceback,
                diagonal.offset + score_index,
            );
        }

        Ok(())
    }

    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    fn clamp_aligned<K: Kernel>(
        kernel: K,
        zero: K::Vector,
        aligned: K::Vector,
        predecessor: K::Vector,
    ) -> (K::Vector, K::Vector) {
        let positive = kernel.greater_than(aligned, zero);

        (
            kernel.select(positive, aligned, zero),
            kernel.select(positive, predecessor, zero),
        )
    }

    #[inline]
    fn track_endpoint<S: LaneScore>(
        endpoint: &mut Option<Endpoint>,
        diagonal_index: usize,
        diagonal: Diagonal,
        current: &DiagonalScores<S>,
        (interior_start, interior_end): (usize, usize),
    ) {
        for score_index in interior_start..interior_end {
            let row = diagonal.first_row + score_index;
            let column = diagonal_index - row;
            update_local_endpoint(
                endpoint,
                current.aligned[score_index].into_score(),
                row,
                column,
            );
        }
    }

    #[inline]
    fn finish<S: LaneScore>(
        endpoint: Option<Endpoint>,
        _layout: &Layout,
        _last: &DiagonalScores<S>,
    ) -> Option<Endpoint> {
        endpoint
    }
}

/// Updates the local endpoint using score, row, and column tie order.
#[inline]
pub(super) fn update_local_endpoint(
    endpoint: &mut Option<Endpoint>,
    score: Score,
    row: usize,
    column: usize,
) {
    if score <= 0 {
        return;
    }

    let candidate = (score, Reverse(row), Reverse(column));
    let replace = endpoint
        .as_ref()
        .map(|current| candidate > (current.score, Reverse(current.row), Reverse(current.column)))
        .unwrap_or(true);

    if replace {
        *endpoint = Some(Endpoint {
            score,
            state: State::Aligned,
            row,
            column,
        });
    }
}

/// Writes scores and predecessor fields for one cell in the current diagonal.
#[inline]
pub(super) fn set_cell<S: LaneScore>(
    scores: &mut DiagonalScores<S>,
    score_index: usize,
    aligned: (S, Option<State>),
    insertion: (S, Option<State>),
    deletion: (S, Option<State>),
    traceback: &mut [u8],
    traceback_index: usize,
) {
    scores.aligned[score_index] = aligned.0;
    scores.insertion[score_index] = insertion.0;
    scores.deletion[score_index] = deletion.0;
    traceback[traceback_index] = pack_predecessors(aligned.1, insertion.1, deletion.1);
}

/// Chooses the best candidate for a gap `state` in one boundary cell whose
/// neighbour is at `index` in `source`.
#[inline]
pub(super) fn choose_boundary<S: LaneScore>(
    state: State,
    source: &DiagonalScores<S>,
    index: usize,
    lane_scoring: &LaneScoring<S>,
) -> Result<(S, Option<State>), Error> {
    choose_scalar(state.candidates().map(|(predecessor, delta)| {
        let delta = match delta {
            Delta::GapOpen => lane_scoring.gap_open_score,
            Delta::GapExtend => lane_scoring.gap_extend_score,
            Delta::Substitution => unreachable!("boundary cells hold only gap states"),
        };
        (source.get(predecessor)[index], delta, predecessor)
    }))
}

/// Evaluates three scalar candidates in the required deterministic order.
#[inline]
pub(super) fn choose_scalar<S: LaneScore>(
    candidates: [(S, S, State); 3],
) -> Result<(S, Option<State>), Error> {
    let candidates = candidates.map(|(predecessor, delta, state)| {
        scalar_candidate(predecessor, delta).map(|score| (score, state))
    });
    let [first, second, third] = candidates;

    Ok(select_score([first?, second?, third?]))
}

/// Adds a delta to a reachable scalar lane.
#[inline]
pub(super) fn scalar_candidate<S: LaneScore>(predecessor: S, delta: S) -> Result<S, Error> {
    if predecessor == S::UNREACHABLE {
        return Ok(S::UNREACHABLE);
    }

    let score = predecessor
        .into_score()
        .checked_add(delta.into_score())
        .ok_or(Error::ScoreOverflow)?;
    narrow(score)
}

/// Selects the greatest candidate and preserves candidate order on ties.
#[inline]
pub(super) fn select_score<S: LaneScore>(candidates: [(S, State); 3]) -> (S, Option<State>) {
    let mut best_score = S::UNREACHABLE;
    let mut best_predecessor = None;

    for (candidate, predecessor) in candidates {
        if candidate != S::UNREACHABLE && (best_score == S::UNREACHABLE || candidate > best_score) {
            best_score = candidate;
            best_predecessor = Some(predecessor);
        }
    }

    (best_score, best_predecessor)
}
