//! The anti-diagonal wavefront engine shared by every SIMD backend.
//!
//! # The recurrence
//!
//! Affine-gap alignment fills a matrix with one row per reference byte and
//! one column per query byte, plus row zero and column zero for the empty
//! prefixes. Cell `(row, column)` holds three scores for aligning
//! `reference[..row]` with `query[..column]`, one per [`State`] the
//! alignment ends in:
//!
//! - aligned: `reference[row - 1]` against `query[column - 1]`, extending the
//!   diagonal neighbour `(row - 1, column - 1)`;
//! - insertion: `query[column - 1]` against a gap, extending the left neighbour
//!   `(row, column - 1)`;
//! - deletion: `reference[row - 1]` against a gap, extending the up neighbour
//!   `(row - 1, column)`.
//!
//! [`State::candidates`] lists each state's predecessors in tie-break order.
//!
//! # Anti-diagonals
//!
//! Anti-diagonal `d` holds every cell with `row + column == d`. A cell's left
//! and up neighbours lie on diagonal `d - 1` and its diagonal neighbour on
//! `d - 2`, never on `d` itself, so all cells of one anti-diagonal are
//! independent and can be computed together. The engine sweeps the diagonals
//! in order and computes `K::LANES` cells of a diagonal per vector. Below,
//! each cell shows its diagonal number; `X` is on diagonal 5, `L` and `U` are
//! on diagonal 4, and `D` is on diagonal 3.
//!
//! ```text
//!              column (query) ->
//!               0   1   2   3   4
//!             +---+---+---+---+---+
//!   row     0 | 0 | 1 | 2 | 3 | 4 |
//!   (ref)     +---+---+---+---+---+
//!    |      1 | 1 | 2 | D | U | 5 |
//!    v        +---+---+---+---+---+
//!           2 | 2 | 3 | L | X | 6 |
//!             +---+---+---+---+---+
//! ```
//!
//! Only three diagonals of scores are kept: `current` (being computed),
//! `previous` (the left and up neighbours), and `two_back` (the diagonal
//! neighbours). [`Wavefront::rotate`] shifts them after every diagonal and
//! reuses the oldest buffer.
//!
//! # Diagonal-major layout
//!
//! A diagonal stores its cells in increasing row order, starting at its first
//! valid row. Moving one cell down a diagonal adds one to the row and
//! subtracts one from the column, which moves every neighbour one index
//! forward too. A vector's left, up, and diagonal neighbours are therefore
//! contiguous on their own diagonals, and [`Cursor`] computes their indexes
//! once per diagonal. The traceback matrix concatenates all diagonals in the
//! same order ([`Layout`]), so a vector's traceback bytes are contiguous
//! too.
//!
//! # Reversed query
//!
//! Lane `k` of a vector starting at `(row, column)` needs
//! `reference[row - 1 + k]` but `query[column - 1 - k]`, which descends.
//! The engine stores the query reversed, so the second byte is
//! `query_reversed[query.len() - column + k]`, and both sequences are read
//! with ascending contiguous loads.
//!
//! # Narrow lanes and the unreachable sentinel
//!
//! Every path through the matrix takes at most
//! `reference.len() + query.len()` steps and adds one scoring parameter per
//! step. No reachable score's magnitude can exceed that step count times the
//! largest parameter magnitude. [`lane_width`] picks the narrowest lane
//! (`i16`, then `i32`) whose maximum covers this bound and leaves every
//! other alignment to the scalar engine, so reachable scores never overflow.
//!
//! The lane's minimum value is reserved as [`LaneScore::UNREACHABLE`] for
//! states that no path reaches, such as an insertion in column zero. A
//! candidate built from an unreachable predecessor stays unreachable instead
//! of adding a delta. Every reachable score is at least `-MAX`, which is
//! greater than the sentinel, so the sentinel loses every comparison with a
//! reachable score.
//!
//! # Padding instead of a scalar tail
//!
//! A diagonal's interior rarely fills a whole number of vectors. Instead of
//! finishing with scalar code, the final vector computes all `K::LANES`
//! lanes. Every diagonal buffer has `K::LANES - 1` slack scores and both
//! sequence copies have `K::LANES - 1` padding bytes, so the overhanging
//! lanes load and store in bounds. Their scores are arbitrary (the kernel's
//! addition wraps), but valid lanes never read them and only the valid lanes'
//! traceback bytes are written. The overhang can overwrite the column-zero
//! cell at the end of a diagonal, so the boundary cells are written after
//! the interior.
//!
//! # Modes
//!
//! [`compute`] is the sweep shared by both modes. [`Global`] and [`Local`]
//! implement [`AlignmentMode`]: they write the row-zero and column-zero
//! boundary cells, clamp aligned scores at zero (local only), and track the
//! endpoint. [`global`] and [`local`] run the sweep and then follow the
//! traceback from the endpoint.
//!
//! # Inlining
//!
//! The hot generic helpers live in submodules, which rustc may compile in
//! other codegen units than [`compute`]. They are marked `#[inline]` so that
//! the sweep still compiles to one loop. On Linux `x86_64` they are
//! `#[inline(always)]` instead, because they must be inlined into the
//! `avx2`-enabled entry point to use AVX2 instructions.

#![cfg_attr(
    all(
        not(test),
        not(any(
            all(target_arch = "aarch64", target_os = "macos"),
            all(target_arch = "x86_64", target_os = "linux")
        ))
    ),
    expect(
        dead_code,
        reason = "targets without a SIMD backend use only the scalar engine outside tests"
    )
)]

/// The lane-wise kernel contract and the vectorized recurrence step.
mod kernel;

/// Narrow score lanes and the per-state containers built from them.
mod lane;

/// The diagonal-major matrix layout and the per-vector index cursor.
mod layout;

/// The global and local alignment modes plugged into the shared sweep.
mod mode;

/// Differential and fixture tests run against portable array kernels
/// and, where available, the native kernels.
#[cfg(test)]
mod tests;

/// The packed traceback byte and the diagonal-major traceback matrix.
mod traceback;

use kernel::BroadcastScoring;
pub(super) use kernel::Kernel;
use kernel::choose_vector;
use kernel::load_states;
use lane::ByDelta;
use lane::ByNeighbour;
use lane::DiagonalScores;
pub(super) use lane::LaneScore;
use lane::LaneScoring;
use lane::allocated_scores;
use layout::Cursor;
use layout::Layout;
use layout::interior_range;
use mode::AlignmentMode;
use mode::Endpoint;
use mode::Global;
use mode::Local;
use traceback::TracebackMatrix;

use crate::algorithm::Error;
use crate::algorithm::Outcome;
use crate::algorithm::Scoring;
use crate::algorithm::engine;
use crate::algorithm::recurrence::State;
use crate::algorithm::recurrence::traceback;
use crate::cigar::OperationKind;

/// The narrowest score lane proven safe for one alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum LaneWidth {
    /// A sixteen-bit signed lane is sufficient.
    I16,
    /// A thirty-two-bit signed lane is sufficient.
    I32,
    /// No SIMD lane proof is available.
    Scalar,
}

impl LaneWidth {
    /// Returns the lane width's name for structured log fields.
    pub(super) const fn name(self) -> &'static str {
        match self {
            Self::I16 => "i16",
            Self::I32 => "i32",
            Self::Scalar => "scalar",
        }
    }
}

/// Chooses the lane width for one alignment and traces the decision.
///
/// `mode` is `"global"` or `"local"`. This is [`lane_width`] plus logging for
/// the backend dispatchers.
pub(super) fn select_lane_width(
    mode: &'static str,
    reference_len: usize,
    query_len: usize,
    scoring: Scoring,
) -> LaneWidth {
    let width = lane_width(reference_len, query_len, scoring);
    tracing::trace!(
        mode,
        reference_len,
        query_len,
        width = width.name(),
        "selected SIMD lane width"
    );
    if width == LaneWidth::Scalar {
        tracing::trace!(
            mode,
            reason = "score bound exceeds every SIMD lane",
            "falling back to the scalar engine"
        );
    }
    width
}

/// Chooses the narrowest lane width whose full-path bound fits.
pub(super) fn lane_width(reference_len: usize, query_len: usize, scoring: Scoring) -> LaneWidth {
    let Some(bound) = score_bound(reference_len, query_len, scoring) else {
        return LaneWidth::Scalar;
    };

    if bound <= i16::MAX as u64 {
        LaneWidth::I16
    } else if bound <= i32::MAX as u64 {
        LaneWidth::I32
    } else {
        LaneWidth::Scalar
    }
}

/// Returns the checked full-path magnitude bound for one alignment.
fn score_bound(reference_len: usize, query_len: usize, scoring: Scoring) -> Option<u64> {
    let steps = reference_len.checked_add(query_len)?;
    let magnitude = [
        scoring.match_score().unsigned_abs(),
        scoring.mismatch_score().unsigned_abs(),
        scoring.gap_open_score().unsigned_abs(),
        scoring.gap_extend_score().unsigned_abs(),
    ]
    .into_iter()
    .fold(0u64, u64::max);
    let Ok(steps) = u64::try_from(steps) else {
        return None;
    };
    steps.checked_mul(magnitude)
}

/// Returns whether the full-path proof permits this narrow score type.
fn kernel_is_eligible<S: LaneScore>(
    reference_len: usize,
    query_len: usize,
    scoring: Scoring,
) -> bool {
    score_bound(reference_len, query_len, scoring)
        .map(|bound| bound <= S::MAX.into_score().unsigned_abs())
        .unwrap_or(false)
}

/// Computes a checked global alignment through a generic wavefront kernel.
#[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
pub(super) fn global<K: Kernel>(
    kernel: K,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Outcome, Error> {
    if reference.is_empty() && query.is_empty() {
        return Err(Error::EmptyGlobal);
    }
    if !kernel_is_eligible::<K::Score>(reference.len(), query.len(), scoring) {
        tracing::trace!(
            mode = "global",
            reason = "score bound exceeds the kernel lane",
            "falling back to the scalar engine"
        );
        return engine::global(reference, query, scoring);
    }

    let (matrix, endpoint) = compute::<K, Global>(kernel, reference, query, scoring)?;
    let Some(endpoint) = endpoint else {
        unreachable!("a non-empty global matrix always has a reachable endpoint")
    };
    let (path, reference_start, query_start) = traceback(
        &matrix,
        reference,
        query,
        endpoint.state,
        endpoint.row,
        endpoint.column,
    )?;
    debug_assert_eq!(reference_start, 0);
    debug_assert_eq!(query_start, 0);
    let cigar = engine::path_to_cigar(&path)?;

    Ok(Outcome::new(
        endpoint.score,
        cigar,
        0..reference.len(),
        0..query.len(),
    ))
}

/// Computes a checked local alignment through a generic wavefront kernel.
#[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
pub(super) fn local<K: Kernel>(
    kernel: K,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    if reference.is_empty() || query.is_empty() {
        return Ok(None);
    }
    if !kernel_is_eligible::<K::Score>(reference.len(), query.len(), scoring) {
        tracing::trace!(
            mode = "local",
            reason = "score bound exceeds the kernel lane",
            "falling back to the scalar engine"
        );
        return engine::local(reference, query, scoring);
    }

    let (matrix, endpoint) = compute::<K, Local>(kernel, reference, query, scoring)?;
    let Some(endpoint) = endpoint else {
        return Ok(None);
    };
    let (path, reference_start, query_start) = traceback(
        &matrix,
        reference,
        query,
        endpoint.state,
        endpoint.row,
        endpoint.column,
    )?;
    debug_assert!(matches!(
        path.first(),
        Some(OperationKind::SequenceMatch | OperationKind::SequenceMismatch)
    ));
    debug_assert!(matches!(
        path.last(),
        Some(OperationKind::SequenceMatch | OperationKind::SequenceMismatch)
    ));
    let cigar = engine::path_to_cigar(&path)?;

    Ok(Some(Outcome::new(
        endpoint.score,
        cigar,
        reference_start..endpoint.row,
        query_start..endpoint.column,
    )))
}

/// Fills all anti-diagonals and returns their packed traceback matrix.
#[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
#[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
fn compute<K: Kernel, M: AlignmentMode>(
    kernel: K,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<(TracebackMatrix, Option<Endpoint>), Error> {
    engine::validate_lengths(reference.len(), query.len())?;
    let (rows, columns) = engine::matrix_dimensions(reference.len(), query.len())?;
    let lane_scoring = LaneScoring::<K::Score>::try_new(scoring)?;
    let mut traceback = TracebackMatrix::try_new(Layout::try_new(rows, columns)?)?;
    let mut wavefront =
        Wavefront::try_new(kernel, reference, query, lane_scoring, rows.min(columns))?;
    let mut endpoint = None;

    for diagonal_index in 0..traceback.layout.diagonals.len() {
        let diagonal = traceback.layout.diagonals[diagonal_index];
        let interior = interior_range(diagonal_index, diagonal);

        if let Some((interior_start, interior_end)) = interior {
            let mut cursor = Cursor::first(&traceback.layout, diagonal_index, interior_start);

            while cursor.current < interior_end {
                let valid = K::LANES.min(interior_end - cursor.current);
                wavefront.fill_vector::<M>(cursor, valid, &mut traceback.entries);
                cursor = cursor.advance(K::LANES);
            }
        }

        // The boundary cells follow the interior because a final partial
        // vector may overwrite the column-zero cell at the end of the
        // diagonal.
        M::write_boundaries(
            diagonal_index,
            &traceback.layout,
            &mut wavefront.current,
            &wavefront.previous,
            &wavefront.lane_scoring,
            &mut traceback.entries,
        )?;

        if let Some(interior) = interior {
            M::track_endpoint(
                &mut endpoint,
                diagonal_index,
                diagonal,
                &wavefront.current,
                interior,
            );
        }

        wavefront.rotate();
    }

    let endpoint = M::finish(endpoint, &traceback.layout, &wavefront.previous);
    Ok((traceback, endpoint))
}

/// The working state of one sweep: the padded inputs, the scoring constants,
/// and the three most recent anti-diagonals of scores.
///
/// A cell on `current` takes its left and up neighbours from `previous` and
/// its diagonal neighbour from `two_back`. Every older diagonal is no longer
/// needed, so [`Wavefront::rotate`] reuses the three buffers in turn.
struct Wavefront<K: Kernel> {
    /// The kernel computing every vector.
    kernel: K,
    /// The reference followed by `K::LANES - 1` padding bytes.
    reference: Vec<u8>,
    /// The reversed query followed by `K::LANES - 1` padding bytes.
    query_reversed: Vec<u8>,
    /// The scoring parameters in the kernel's lane type.
    lane_scoring: LaneScoring<K::Score>,
    /// The recurrence constants broadcast across every lane.
    broadcast: BroadcastScoring<K>,
    /// The scores of the diagonal being computed.
    current: DiagonalScores<K::Score>,
    /// The scores of the diagonal before `current`.
    previous: DiagonalScores<K::Score>,
    /// The scores of the diagonal before `previous`.
    two_back: DiagonalScores<K::Score>,
    /// Scratch for one vector of packed traceback bytes.
    packed: Vec<K::Score>,
}

impl<K: Kernel> Wavefront<K> {
    /// Allocates the padded inputs and the diagonal buffers for anti-diagonals
    /// of at most `max_diagonal_len` cells.
    ///
    /// Each buffer has `K::LANES - 1` slack scores after the longest diagonal
    /// and each sequence copy has `K::LANES - 1` padding bytes, so a final
    /// partial vector can load and store all `K::LANES` lanes in bounds.
    fn try_new(
        kernel: K,
        reference: &[u8],
        query: &[u8],
        lane_scoring: LaneScoring<K::Score>,
        max_diagonal_len: usize,
    ) -> Result<Self, Error> {
        if K::LANES == 0 {
            return Err(Error::MatrixSizeOverflow);
        }

        let slack = K::LANES - 1;
        let capacity = max_diagonal_len
            .checked_add(slack)
            .ok_or(Error::MatrixSizeOverflow)?;

        Ok(Self {
            kernel,
            reference: padded_copy(reference, slack, false)?,
            query_reversed: padded_copy(query, slack, true)?,
            broadcast: BroadcastScoring::new(kernel, &lane_scoring),
            lane_scoring,
            current: DiagonalScores::try_new(capacity)?,
            previous: DiagonalScores::try_new(capacity)?,
            two_back: DiagonalScores::try_new(capacity)?,
            packed: allocated_scores(K::LANES)?,
        })
    }

    /// Computes one vector of interior cells starting at `cursor`.
    ///
    /// Only the first `valid` lanes are cells of this diagonal. The other
    /// lanes read and write the diagonal buffers' slack (and possibly the
    /// column-zero boundary cell, which is written afterwards), compute
    /// arbitrary wrapping values, and never reach the traceback matrix. A
    /// valid lane reads only valid neighbours, so overhanging lanes never
    /// influence a result.
    #[cfg_attr(all(target_arch = "x86_64", target_os = "linux"), inline(always))]
    #[cfg_attr(not(all(target_arch = "x86_64", target_os = "linux")), inline)]
    fn fill_vector<M: AlignmentMode>(
        &mut self,
        cursor: Cursor,
        valid: usize,
        traceback: &mut [u8],
    ) {
        let kernel = self.kernel;
        let broadcast = &self.broadcast;

        debug_assert!((1..=K::LANES).contains(&valid));
        debug_assert!(cursor.diagonal + K::LANES <= self.two_back.aligned.len());
        debug_assert!(cursor.left + K::LANES <= self.previous.aligned.len());
        debug_assert!(cursor.reference + K::LANES <= self.reference.len());
        debug_assert!(cursor.query + K::LANES <= self.query_reversed.len());

        // SAFETY: The first lane is a valid interior cell, so its neighbour
        // indexes are below the maximum diagonal length, and its byte indexes
        // are below the sequence lengths. Every diagonal buffer has
        // `K::LANES - 1` slack scores and every sequence copy has
        // `K::LANES - 1` padding bytes, so each `K::LANES`-wide access stays
        // in bounds. `up` is `left - 1`.
        let (neighbours, substitution) = unsafe {
            (
                ByNeighbour {
                    diagonal: load_states(kernel, &self.two_back, cursor.diagonal),
                    left: load_states(kernel, &self.previous, cursor.left),
                    up: load_states(kernel, &self.previous, cursor.up),
                },
                kernel.substitution(
                    self.reference.as_ptr().add(cursor.reference),
                    self.query_reversed.as_ptr().add(cursor.query),
                    self.lane_scoring.match_score,
                    self.lane_scoring.mismatch_score,
                ),
            )
        };
        let deltas = ByDelta {
            substitution,
            gap_open: broadcast.gap_open_score,
            gap_extend: broadcast.gap_extend_score,
        };
        let (aligned, aligned_predecessor) =
            choose_vector(kernel, State::Aligned, &neighbours, &deltas, broadcast);
        let (insertion, insertion_predecessor) =
            choose_vector(kernel, State::Insertion, &neighbours, &deltas, broadcast);
        let (deletion, deletion_predecessor) =
            choose_vector(kernel, State::Deletion, &neighbours, &deltas, broadcast);
        let (aligned, aligned_predecessor) =
            M::clamp_aligned(kernel, broadcast.zero, aligned, aligned_predecessor);

        let packed_predecessors = kernel.add(
            kernel.add(aligned_predecessor, insertion_predecessor),
            deletion_predecessor,
        );

        debug_assert!(cursor.current + K::LANES <= self.current.aligned.len());
        debug_assert_eq!(self.packed.len(), K::LANES);

        // SAFETY: `cursor.current` is a valid index on the current diagonal
        // and the buffer has `K::LANES - 1` slack scores after the longest
        // diagonal. `packed` holds exactly `K::LANES` scores.
        unsafe {
            let current = &mut self.current;
            kernel.store(current.aligned.as_mut_ptr().add(cursor.current), aligned);
            kernel.store(
                current.insertion.as_mut_ptr().add(cursor.current),
                insertion,
            );
            kernel.store(current.deletion.as_mut_ptr().add(cursor.current), deletion);
            kernel.store(self.packed.as_mut_ptr(), packed_predecessors);
        }

        for (entry, lane) in traceback[cursor.traceback..cursor.traceback + valid]
            .iter_mut()
            .zip(self.packed.iter())
        {
            *entry = lane.into_code();
        }
    }

    /// Retires `two_back` and makes the finished `current` diagonal the new
    /// `previous`; the retired buffer becomes the next `current`.
    fn rotate(&mut self) {
        std::mem::swap(&mut self.two_back, &mut self.previous);
        std::mem::swap(&mut self.previous, &mut self.current);
    }
}

/// Copies `sequence` (reversed when `reverse` is set) followed by `padding`
/// zero bytes.
fn padded_copy(sequence: &[u8], padding: usize, reverse: bool) -> Result<Vec<u8>, Error> {
    let length = sequence
        .len()
        .checked_add(padding)
        .ok_or(Error::MatrixSizeOverflow)?;
    let mut copy = Vec::new();
    copy.try_reserve_exact(length)
        .map_err(|source| Error::MatrixAllocation { source })?;

    if reverse {
        copy.extend(sequence.iter().rev());
    } else {
        copy.extend_from_slice(sequence);
    }
    copy.resize(length, 0);

    Ok(copy)
}
