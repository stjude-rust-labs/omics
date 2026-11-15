//! The lane-wise kernel contract and the vectorized recurrence step.
//!
//! A [`Kernel`] supplies eight lane-wise operations (load, store, splat, add,
//! compare, select, and substitution) for one vector width. The engine builds
//! the whole recurrence from them: [`choose_vector`] walks
//! [`State::candidates`] in tie-break order and keeps, per lane, the best
//! score together with its predecessor code already shifted into the
//! destination state's traceback field.

use super::lane::ByDelta;
use super::lane::ByNeighbour;
use super::lane::ByState;
use super::lane::DiagonalScores;
use super::lane::LaneScore;
use super::lane::LaneScoring;
use super::traceback::field_code;
use crate::algorithm::recurrence::State;

/// A lane-wise backend for the anti-diagonal affine-gap recurrence.
///
/// # Safety
///
/// Implementations use `LANES` greater than zero and preserve the lane-wise
/// semantics of every method. The load and store methods access exactly
/// `LANES` scores. The substitution method reads exactly `LANES` ascending
/// bytes from each of its two pointers. Addition wraps on overflow, because
/// lanes past the end of a diagonal hold arbitrary values (see
/// [`Wavefront::fill_vector`](super::Wavefront::fill_vector)).
pub(in crate::algorithm::simd) unsafe trait Kernel: Copy {
    /// The narrow score representation used by each lane.
    type Score: LaneScore;
    /// The vector or scalar-array representation containing every lane.
    type Vector: Copy;

    /// The number of independently processed lanes.
    const LANES: usize;

    /// Loads `LANES` contiguous score values.
    ///
    /// # Safety
    ///
    /// `pointer` points to at least `Self::LANES` readable scores.
    unsafe fn load(self, pointer: *const Self::Score) -> Self::Vector;

    /// Stores `LANES` contiguous score values.
    ///
    /// # Safety
    ///
    /// `pointer` points to at least `Self::LANES` writable scores.
    unsafe fn store(self, pointer: *mut Self::Score, value: Self::Vector);

    /// Broadcasts one score into every lane.
    fn splat(self, value: Self::Score) -> Self::Vector;

    /// Adds corresponding lanes with wrapping arithmetic.
    fn add(self, left: Self::Vector, right: Self::Vector) -> Self::Vector;

    /// Returns a lane mask for strict greater-than comparisons.
    fn greater_than(self, left: Self::Vector, right: Self::Vector) -> Self::Vector;

    /// Returns a lane mask for equality comparisons.
    fn equal(self, left: Self::Vector, right: Self::Vector) -> Self::Vector;

    /// Selects one lane from `if_true` or `if_false` for each mask lane.
    fn select(
        self,
        mask: Self::Vector,
        if_true: Self::Vector,
        if_false: Self::Vector,
    ) -> Self::Vector;

    /// Returns `match_score` in each lane whose reference and reversed-query
    /// bytes are equal and `mismatch_score` in every other lane.
    ///
    /// # Safety
    ///
    /// `reference` and `query_reversed` each point to `Self::LANES` readable
    /// bytes.
    unsafe fn substitution(
        self,
        reference: *const u8,
        query_reversed: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector;
}

/// Holds recurrence constants broadcast across every kernel lane.
pub(super) struct BroadcastScoring<K: Kernel> {
    /// Unreachable lane score.
    pub(super) unreachable: K::Vector,
    /// Reachable zero lane score.
    pub(super) zero: K::Vector,
    /// Gap-open score.
    pub(super) gap_open_score: K::Vector,
    /// Gap-extension score.
    pub(super) gap_extend_score: K::Vector,
    /// Predecessor codes indexed by destination state and then predecessor
    /// state, each already shifted into the destination's traceback field.
    ///
    /// Because the three fields do not overlap, adding the chosen code of
    /// each destination state yields the packed traceback byte.
    pub(super) predecessor_codes: ByState<ByState<K::Vector>>,
}

impl<K: Kernel> BroadcastScoring<K> {
    /// Broadcasts narrow recurrence constants before the diagonal loop.
    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    pub(super) fn new(kernel: K, lane_scoring: &LaneScoring<K::Score>) -> Self {
        let codes = |destination| {
            let splat = |predecessor| {
                kernel.splat(K::Score::from_code(field_code(predecessor, destination)))
            };
            ByState {
                aligned: splat(State::Aligned),
                insertion: splat(State::Insertion),
                deletion: splat(State::Deletion),
            }
        };

        Self {
            unreachable: kernel.splat(K::Score::UNREACHABLE),
            zero: kernel.splat(K::Score::ZERO),
            gap_open_score: kernel.splat(lane_scoring.gap_open_score),
            gap_extend_score: kernel.splat(lane_scoring.gap_extend_score),
            predecessor_codes: ByState {
                aligned: codes(State::Aligned),
                insertion: codes(State::Insertion),
                deletion: codes(State::Deletion),
            },
        }
    }
}

/// Forms one candidate score while retaining unreachable lanes.
#[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
#[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
pub(super) fn vector_candidate<K: Kernel>(
    kernel: K,
    predecessor: K::Vector,
    delta: K::Vector,
    unreachable: K::Vector,
) -> K::Vector {
    let candidate = kernel.add(predecessor, delta);

    kernel.select(
        kernel.equal(predecessor, unreachable),
        unreachable,
        candidate,
    )
}

/// Chooses the best candidate for `state` in one vector, following the
/// tie-break order of [`State::candidates`], and returns the score together
/// with the chosen predecessor's code in `state`'s traceback field.
#[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
#[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
pub(super) fn choose_vector<K: Kernel>(
    kernel: K,
    state: State,
    neighbours: &ByNeighbour<ByState<K::Vector>>,
    deltas: &ByDelta<K::Vector>,
    broadcast: &BroadcastScoring<K>,
) -> (K::Vector, K::Vector) {
    let sources = neighbours.get(state.neighbour());
    let codes = broadcast.predecessor_codes.get(state);
    let mut best_score = broadcast.unreachable;
    let mut best_predecessor = broadcast.zero;

    for (predecessor, delta) in state.candidates() {
        let candidate = vector_candidate(
            kernel,
            *sources.get(predecessor),
            *deltas.get(delta),
            broadcast.unreachable,
        );
        let replace = kernel.greater_than(candidate, best_score);
        best_score = kernel.select(replace, candidate, best_score);
        best_predecessor = kernel.select(replace, *codes.get(predecessor), best_predecessor);
    }

    (best_score, best_predecessor)
}

/// Loads one vector of all three state scores starting at `index`.
///
/// # Safety
///
/// Each of `scores`' three arrays holds at least `index + K::LANES` scores.
#[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
#[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
pub(super) unsafe fn load_states<K: Kernel>(
    kernel: K,
    scores: &DiagonalScores<K::Score>,
    index: usize,
) -> ByState<K::Vector> {
    // SAFETY: The caller guarantees `K::LANES` readable scores at `index`.
    unsafe {
        ByState {
            aligned: kernel.load(scores.aligned.as_ptr().add(index)),
            insertion: kernel.load(scores.insertion.as_ptr().add(index)),
            deletion: kernel.load(scores.deletion.as_ptr().add(index)),
        }
    }
}
