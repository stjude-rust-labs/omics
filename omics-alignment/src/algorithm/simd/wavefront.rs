#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "target-specific SIMD backends are added after the generic wavefront engine"
    )
)]

use std::cmp::Reverse;

use crate::algorithm::Error;
use crate::algorithm::Outcome;
use crate::algorithm::Score;
use crate::algorithm::Scoring;
use crate::algorithm::engine;
use crate::cigar::OperationKind;

/// The narrowest score lane proven safe for one alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum LaneWidth {
    /// A sixteen-bit signed lane is sufficient.
    I16,
    /// A thirty-two-bit signed lane is sufficient.
    I32,
    /// No SIMD lane proof is available.
    Scalar,
}

/// Chooses the narrowest lane width whose full-path bound fits.
fn lane_width(reference_len: usize, query_len: usize, scoring: Scoring) -> LaneWidth {
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

/// Metadata for one anti-diagonal in diagonal-major storage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Diagonal {
    /// The first linear index stored on this diagonal.
    offset: usize,
    /// The first valid row on this diagonal.
    first_row: usize,
    /// The number of cells stored on this diagonal.
    len: usize,
}

/// A checked mapping between matrix coordinates and diagonal-major indexes.
struct Layout {
    /// The matrix row count.
    rows: usize,
    /// The matrix column count.
    columns: usize,
    /// The total number of matrix cells.
    cells: usize,
    /// The metadata for each anti-diagonal.
    diagonals: Vec<Diagonal>,
}

impl Layout {
    /// Builds a checked diagonal-major layout for `rows * columns` cells.
    fn try_new(rows: usize, columns: usize) -> Result<Self, Error> {
        if rows == 0 || columns == 0 {
            return Err(Error::MatrixSizeOverflow);
        }

        let cells = rows.checked_mul(columns).ok_or(Error::MatrixSizeOverflow)?;
        let diagonal_count = rows
            .checked_add(columns)
            .and_then(|sum| sum.checked_sub(1))
            .ok_or(Error::MatrixSizeOverflow)?;
        let mut diagonals = Vec::new();
        diagonals
            .try_reserve_exact(diagonal_count)
            .map_err(|source| Error::MatrixAllocation { source })?;

        let mut offset = 0usize;
        for diagonal in 0..diagonal_count {
            let first_row = diagonal.saturating_sub(columns.saturating_sub(1));
            let last_row = diagonal.min(rows.saturating_sub(1));
            let len = last_row
                .checked_sub(first_row)
                .and_then(|span| span.checked_add(1))
                .ok_or(Error::MatrixSizeOverflow)?;
            diagonals.push(Diagonal {
                offset,
                first_row,
                len,
            });
            offset = offset.checked_add(len).ok_or(Error::MatrixSizeOverflow)?;
        }

        debug_assert_eq!(offset, cells);

        Ok(Self {
            rows,
            columns,
            cells,
            diagonals,
        })
    }

    /// Maps one matrix coordinate to its diagonal-major index.
    fn index(&self, row: usize, column: usize) -> usize {
        debug_assert!(row < self.rows);
        debug_assert!(column < self.columns);

        let diagonal = row + column;
        let diagonal = self.diagonals[diagonal];
        diagonal.offset + row - diagonal.first_row
    }

    /// Maps one diagonal-major index back to its matrix coordinate.
    fn coordinate(&self, index: usize) -> (usize, usize) {
        debug_assert!(index < self.cells);

        let diagonal = match self
            .diagonals
            .binary_search_by(|entry| entry.offset.cmp(&index))
        {
            Ok(diagonal) => diagonal,
            Err(diagonal) => {
                debug_assert!(diagonal > 0, "index precedes all diagonal offsets");
                diagonal - 1
            }
        };
        let metadata = self.diagonals[diagonal];
        let row = metadata.first_row + (index - metadata.offset);

        debug_assert!(row < metadata.first_row + metadata.len);

        (row, diagonal - row)
    }
}

/// A score lane that reserves its minimum value for unreachable states.
pub(super) trait Narrow: Copy + Eq + Ord {
    /// The sentinel used for an unreachable score.
    const MIN: Self;
    /// The greatest score that can remain reachable in this lane.
    const MAX: Self;
    /// The score used for a reachable origin or local reset.
    const ZERO: Self;

    /// Converts a scalar score when it fits without using the sentinel.
    fn try_from_score(score: Score) -> Option<Self>;

    /// Converts this lane value to the scalar score representation.
    fn into_score(self) -> Score;
}

impl Narrow for i16 {
    const MAX: Self = i16::MAX;
    const MIN: Self = i16::MIN;
    const ZERO: Self = 0;

    fn try_from_score(score: Score) -> Option<Self> {
        let score = Self::try_from(score).ok()?;
        (score != i16::MIN).then_some(score)
    }

    fn into_score(self) -> Score {
        Score::from(self)
    }
}

impl Narrow for i32 {
    const MAX: Self = i32::MAX;
    const MIN: Self = i32::MIN;
    const ZERO: Self = 0;

    fn try_from_score(score: Score) -> Option<Self> {
        let score = Self::try_from(score).ok()?;
        (score != i32::MIN).then_some(score)
    }

    fn into_score(self) -> Score {
        Score::from(self)
    }
}

/// A lane-wise backend for the anti-diagonal affine-gap recurrence.
///
/// # Safety
///
/// Implementations use `LANES` greater than zero and preserve the lane-wise
/// semantics of every method. The load and store methods access exactly
/// `LANES` scores. The substitution method reads `LANES` ascending bytes from
/// `reference` and `LANES` descending bytes immediately before `query_end`.
pub(super) unsafe trait Kernel {
    /// The narrow score representation used by each lane.
    type Score: Narrow;
    /// The vector or scalar-array representation containing every lane.
    type Vector: Copy;

    /// The number of independently processed lanes.
    const LANES: usize;

    /// Loads `LANES` contiguous score values.
    ///
    /// # Safety
    ///
    /// `pointer` points to at least `Self::LANES` readable scores.
    unsafe fn load(pointer: *const Self::Score) -> Self::Vector;

    /// Stores `LANES` contiguous score values.
    ///
    /// # Safety
    ///
    /// `pointer` points to at least `Self::LANES` writable scores.
    unsafe fn store(pointer: *mut Self::Score, value: Self::Vector);

    /// Broadcasts one score into every lane.
    fn splat(value: Self::Score) -> Self::Vector;

    /// Adds corresponding lanes with the backend's native arithmetic.
    fn add(left: Self::Vector, right: Self::Vector) -> Self::Vector;

    /// Returns a lane mask for strict greater-than comparisons.
    fn greater_than(left: Self::Vector, right: Self::Vector) -> Self::Vector;

    /// Returns a lane mask for equality comparisons.
    fn equal(left: Self::Vector, right: Self::Vector) -> Self::Vector;

    /// Selects one lane from `if_true` or `if_false` for each mask lane.
    fn select(mask: Self::Vector, if_true: Self::Vector, if_false: Self::Vector) -> Self::Vector;

    /// Computes substitutions from ascending reference and descending query
    /// bytes.
    ///
    /// # Safety
    ///
    /// `reference` points to `Self::LANES` readable bytes. `query_end` points
    /// one byte beyond `Self::LANES` readable bytes that precede it.
    unsafe fn substitution(
        reference: *const u8,
        query_end: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector;
}

/// One terminal state in the affine-gap recurrence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum State {
    /// Ends by consuming one symbol from each input.
    Aligned = 1,
    /// Ends by consuming one query symbol.
    Insertion = 2,
    /// Ends by consuming one reference symbol.
    Deletion = 3,
}

impl State {
    /// Returns this state's two-bit predecessor field offset.
    const fn field_shift(self) -> u8 {
        match self {
            Self::Aligned => 0,
            Self::Insertion => 2,
            Self::Deletion => 4,
        }
    }

    /// Decodes one two-bit predecessor code.
    const fn from_code(code: u8) -> Option<Self> {
        match code {
            1 => Some(Self::Aligned),
            2 => Some(Self::Insertion),
            3 => Some(Self::Deletion),
            _ => None,
        }
    }
}

/// Stores all three state scores for one retained anti-diagonal.
struct Scores<S> {
    /// Aligned-state scores.
    aligned: Vec<S>,
    /// Insertion-state scores.
    insertion: Vec<S>,
    /// Deletion-state scores.
    deletion: Vec<S>,
}

impl<S: Narrow> Scores<S> {
    /// Allocates unreachable scores for one maximum-width anti-diagonal.
    fn try_new(length: usize) -> Result<Self, Error> {
        Ok(Self {
            aligned: allocated_scores(length)?,
            insertion: allocated_scores(length)?,
            deletion: allocated_scores(length)?,
        })
    }

    /// Marks every retained score as unreachable before reusing the diagonal.
    fn clear(&mut self) {
        self.aligned.fill(S::MIN);
        self.insertion.fill(S::MIN);
        self.deletion.fill(S::MIN);
    }
}

/// Allocates one unreachable score array.
fn allocated_scores<S: Narrow>(length: usize) -> Result<Vec<S>, Error> {
    let mut scores = Vec::new();
    scores
        .try_reserve_exact(length)
        .map_err(|source| Error::MatrixAllocation { source })?;
    scores.resize(length, S::MIN);
    Ok(scores)
}

/// Holds state-code vectors while their traceback bytes are packed.
struct PredecessorScratch<S> {
    /// Aligned-state predecessor codes.
    aligned: Vec<S>,
    /// Insertion-state predecessor codes.
    insertion: Vec<S>,
    /// Deletion-state predecessor codes.
    deletion: Vec<S>,
}

impl<S: Narrow> PredecessorScratch<S> {
    /// Allocates storage for one complete kernel vector.
    fn try_new(lanes: usize) -> Result<Self, Error> {
        Ok(Self {
            aligned: allocated_scores(lanes)?,
            insertion: allocated_scores(lanes)?,
            deletion: allocated_scores(lanes)?,
        })
    }
}

/// Stores narrow scoring values and encoded predecessor state values.
struct Values<S> {
    /// Score for an equal byte pair.
    match_score: S,
    /// Score for an unequal byte pair.
    mismatch_score: S,
    /// Score for opening a gap.
    gap_open_score: S,
    /// Score for extending a gap.
    gap_extend_score: S,
    /// Encoded aligned-state predecessor.
    aligned_code: S,
    /// Encoded insertion-state predecessor.
    insertion_code: S,
    /// Encoded deletion-state predecessor.
    deletion_code: S,
}

impl<S: Narrow> Values<S> {
    /// Converts all scalar scoring values and state codes into narrow lanes.
    fn try_new(scoring: Scoring) -> Result<Self, Error> {
        Ok(Self {
            match_score: narrow(scoring.match_score())?,
            mismatch_score: narrow(scoring.mismatch_score())?,
            gap_open_score: narrow(scoring.gap_open_score())?,
            gap_extend_score: narrow(scoring.gap_extend_score())?,
            aligned_code: narrow(Score::from(State::Aligned as u8))?,
            insertion_code: narrow(Score::from(State::Insertion as u8))?,
            deletion_code: narrow(Score::from(State::Deletion as u8))?,
        })
    }
}

/// Converts one scalar score into a usable narrow lane.
fn narrow<S: Narrow>(score: Score) -> Result<S, Error> {
    S::try_from_score(score).ok_or(Error::ScoreOverflow)
}

/// Returns whether the full-path proof permits this narrow score type.
fn kernel_is_eligible<S: Narrow>(reference_len: usize, query_len: usize, scoring: Scoring) -> bool {
    score_bound(reference_len, query_len, scoring)
        .map(|bound| bound <= S::MAX.into_score().unsigned_abs())
        .unwrap_or(false)
}

/// Packs the three predecessor state codes for one matrix cell.
fn pack_predecessors(
    aligned: Option<State>,
    insertion: Option<State>,
    deletion: Option<State>,
) -> u8 {
    let mut packed = 0;

    for (state, predecessor) in [
        (State::Aligned, aligned),
        (State::Insertion, insertion),
        (State::Deletion, deletion),
    ] {
        let code = predecessor.map_or(0, |predecessor| predecessor as u8);
        packed |= code << state.field_shift();
    }

    packed
}

/// Decodes one packed predecessor field.
fn unpack_predecessor(packed: u8, state: State) -> Option<State> {
    State::from_code((packed >> state.field_shift()) & 0b11)
}

/// Converts one kernel predecessor code into a state.
fn state_from_lane<S: Narrow>(value: S) -> Option<State> {
    u8::try_from(value.into_score())
        .ok()
        .and_then(State::from_code)
}

/// Selects the greatest candidate and preserves candidate order on ties.
fn select_score<S: Narrow>(candidates: [(S, State); 3]) -> (S, Option<State>) {
    let mut best_score = S::MIN;
    let mut best_predecessor = None;

    for (candidate, predecessor) in candidates {
        if candidate != S::MIN && (best_score == S::MIN || candidate > best_score) {
            best_score = candidate;
            best_predecessor = Some(predecessor);
        }
    }

    (best_score, best_predecessor)
}

/// Adds a delta to a reachable scalar lane.
fn scalar_candidate<S: Narrow>(predecessor: S, delta: S) -> Result<S, Error> {
    if predecessor == S::MIN {
        return Ok(S::MIN);
    }

    let score = predecessor
        .into_score()
        .checked_add(delta.into_score())
        .ok_or(Error::ScoreOverflow)?;
    narrow(score)
}

/// Evaluates three scalar candidates in the required deterministic order.
fn choose_scalar<S: Narrow>(candidates: [(S, S, State); 3]) -> Result<(S, Option<State>), Error> {
    let candidates = candidates.map(|(predecessor, delta, state)| {
        scalar_candidate(predecessor, delta).map(|score| (score, state))
    });
    let [first, second, third] = candidates;

    Ok(select_score([first?, second?, third?]))
}

/// Forms one candidate score while retaining unreachable lanes.
fn vector_candidate<K: Kernel>(predecessor: K::Vector, delta: K::Vector) -> K::Vector {
    let unreachable = K::splat(K::Score::MIN);
    let candidate = K::add(predecessor, delta);

    K::select(K::equal(predecessor, unreachable), unreachable, candidate)
}

/// Evaluates three vector candidates in the required deterministic order.
fn choose_vector<K: Kernel>(
    candidates: [(K::Vector, K::Vector, K::Score); 3],
) -> (K::Vector, K::Vector) {
    let unreachable = K::splat(K::Score::MIN);
    let mut best_score = unreachable;
    let mut best_predecessor = K::splat(K::Score::ZERO);

    for (predecessor, delta, state) in candidates {
        let candidate = vector_candidate::<K>(predecessor, delta);
        let replace = K::greater_than(candidate, best_score);
        best_score = K::select(replace, candidate, best_score);
        best_predecessor = K::select(replace, K::splat(state), best_predecessor);
    }

    (best_score, best_predecessor)
}

/// Selects global or local recurrence initialization and endpoint behavior.
#[derive(Clone, Copy, Eq, PartialEq)]
enum Mode {
    /// Aligns both complete input slices.
    Global,
    /// Aligns the highest-scoring positive subsequences.
    Local,
}

/// Records the terminal state and coordinate selected by an alignment mode.
#[derive(Clone, Copy)]
struct Endpoint {
    /// Terminal alignment score.
    score: Score,
    /// Terminal recurrence state.
    state: State,
    /// Exclusive reference coordinate.
    row: usize,
    /// Exclusive query coordinate.
    column: usize,
}

/// Owns diagonal-major traceback fields for one computed matrix.
struct TracebackMatrix {
    /// Coordinate layout for every traceback field.
    layout: Layout,
    /// Packed two-bit predecessor fields for every matrix cell.
    entries: Vec<u8>,
}

impl TracebackMatrix {
    /// Allocates empty predecessor fields for a checked layout.
    fn try_new(layout: Layout) -> Result<Self, Error> {
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(layout.cells)
            .map_err(|source| Error::MatrixAllocation { source })?;
        entries.resize(layout.cells, 0);

        Ok(Self { layout, entries })
    }
}

/// Writes scores and predecessor fields for one cell in the current diagonal.
fn set_cell<S: Narrow>(
    scores: &mut Scores<S>,
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

/// Packs state-code lanes emitted by one complete kernel vector.
fn pack_vector_predecessors<S: Narrow>(traceback: &mut [u8], predecessors: &PredecessorScratch<S>) {
    debug_assert_eq!(traceback.len(), predecessors.aligned.len());

    for (lane, packed) in traceback.iter_mut().enumerate() {
        *packed = pack_predecessors(
            state_from_lane(predecessors.aligned[lane]),
            state_from_lane(predecessors.insertion[lane]),
            state_from_lane(predecessors.deletion[lane]),
        );
    }
}

/// Computes one complete interior vector of the recurrence.
#[expect(
    clippy::too_many_arguments,
    reason = "the recurrence needs all three diagonal states and coordinate offsets"
)]
fn fill_vector<K: Kernel>(
    current: &mut Scores<K::Score>,
    previous: &Scores<K::Score>,
    before_previous: &Scores<K::Score>,
    current_index: usize,
    left_index: usize,
    up_index: usize,
    diagonal_index: usize,
    row: usize,
    column: usize,
    reference: &[u8],
    query: &[u8],
    values: &Values<K::Score>,
    mode: Mode,
    predecessors: &mut PredecessorScratch<K::Score>,
    traceback: &mut [u8],
    traceback_index: usize,
) {
    // SAFETY: The caller passes a complete interior vector. Every predecessor
    // range holds `K::LANES` scores, and the input pointers span the matching
    // ascending reference and descending query byte ranges.
    let (
        diagonal_aligned,
        diagonal_insertion,
        diagonal_deletion,
        left_aligned,
        left_insertion,
        left_deletion,
        up_aligned,
        up_insertion,
        up_deletion,
        substitution,
    ) = unsafe {
        (
            K::load(before_previous.aligned.as_ptr().add(diagonal_index)),
            K::load(before_previous.insertion.as_ptr().add(diagonal_index)),
            K::load(before_previous.deletion.as_ptr().add(diagonal_index)),
            K::load(previous.aligned.as_ptr().add(left_index)),
            K::load(previous.insertion.as_ptr().add(left_index)),
            K::load(previous.deletion.as_ptr().add(left_index)),
            K::load(previous.aligned.as_ptr().add(up_index)),
            K::load(previous.insertion.as_ptr().add(up_index)),
            K::load(previous.deletion.as_ptr().add(up_index)),
            K::substitution(
                reference.as_ptr().add(row - 1),
                query.as_ptr().add(column),
                values.match_score,
                values.mismatch_score,
            ),
        )
    };
    let (mut aligned, mut aligned_predecessor) = choose_vector::<K>([
        (diagonal_aligned, substitution, values.aligned_code),
        (diagonal_deletion, substitution, values.deletion_code),
        (diagonal_insertion, substitution, values.insertion_code),
    ]);
    let (insertion, insertion_predecessor) = choose_vector::<K>([
        (
            left_insertion,
            K::splat(values.gap_extend_score),
            values.insertion_code,
        ),
        (
            left_aligned,
            K::splat(values.gap_open_score),
            values.aligned_code,
        ),
        (
            left_deletion,
            K::splat(values.gap_open_score),
            values.deletion_code,
        ),
    ]);
    let (deletion, deletion_predecessor) = choose_vector::<K>([
        (
            up_deletion,
            K::splat(values.gap_extend_score),
            values.deletion_code,
        ),
        (
            up_aligned,
            K::splat(values.gap_open_score),
            values.aligned_code,
        ),
        (
            up_insertion,
            K::splat(values.gap_open_score),
            values.insertion_code,
        ),
    ]);

    if mode == Mode::Local {
        let zero = K::splat(K::Score::ZERO);
        let positive = K::greater_than(aligned, zero);
        aligned = K::select(positive, aligned, zero);
        aligned_predecessor = K::select(positive, aligned_predecessor, zero);
    }

    // SAFETY: The complete interior vector supplies `K::LANES` writable
    // current scores and predecessor scratch entries.
    unsafe {
        K::store(current.aligned.as_mut_ptr().add(current_index), aligned);
        K::store(current.insertion.as_mut_ptr().add(current_index), insertion);
        K::store(current.deletion.as_mut_ptr().add(current_index), deletion);
        K::store(predecessors.aligned.as_mut_ptr(), aligned_predecessor);
        K::store(predecessors.insertion.as_mut_ptr(), insertion_predecessor);
        K::store(predecessors.deletion.as_mut_ptr(), deletion_predecessor);
    }
    pack_vector_predecessors(
        &mut traceback[traceback_index..traceback_index + K::LANES],
        predecessors,
    );
}

/// Computes one scalar tail cell of the recurrence.
#[expect(
    clippy::too_many_arguments,
    reason = "the recurrence needs all three diagonal states and coordinate offsets"
)]
fn fill_scalar<S: Narrow>(
    current: &mut Scores<S>,
    previous: &Scores<S>,
    before_previous: &Scores<S>,
    current_index: usize,
    left_index: usize,
    up_index: usize,
    diagonal_index: usize,
    row: usize,
    column: usize,
    reference: &[u8],
    query: &[u8],
    values: &Values<S>,
    mode: Mode,
    traceback: &mut [u8],
    traceback_index: usize,
) -> Result<(), Error> {
    let substitution = if reference[row - 1] == query[column - 1] {
        values.match_score
    } else {
        values.mismatch_score
    };
    let (aligned_score, aligned_predecessor) = choose_scalar([
        (
            before_previous.aligned[diagonal_index],
            substitution,
            State::Aligned,
        ),
        (
            before_previous.deletion[diagonal_index],
            substitution,
            State::Deletion,
        ),
        (
            before_previous.insertion[diagonal_index],
            substitution,
            State::Insertion,
        ),
    ])?;
    let (insertion_score, insertion_predecessor) = choose_scalar([
        (
            previous.insertion[left_index],
            values.gap_extend_score,
            State::Insertion,
        ),
        (
            previous.aligned[left_index],
            values.gap_open_score,
            State::Aligned,
        ),
        (
            previous.deletion[left_index],
            values.gap_open_score,
            State::Deletion,
        ),
    ])?;
    let (deletion_score, deletion_predecessor) = choose_scalar([
        (
            previous.deletion[up_index],
            values.gap_extend_score,
            State::Deletion,
        ),
        (
            previous.aligned[up_index],
            values.gap_open_score,
            State::Aligned,
        ),
        (
            previous.insertion[up_index],
            values.gap_open_score,
            State::Insertion,
        ),
    ])?;
    let aligned = if mode == Mode::Local && aligned_score <= S::ZERO {
        (S::ZERO, None)
    } else {
        (aligned_score, aligned_predecessor)
    };

    set_cell(
        current,
        current_index,
        aligned,
        (insertion_score, insertion_predecessor),
        (deletion_score, deletion_predecessor),
        traceback,
        traceback_index,
    );
    Ok(())
}

/// Initializes the reachable global boundary cells on one anti-diagonal.
fn initialize_global_boundary<S: Narrow>(
    diagonal_index: usize,
    diagonal: Diagonal,
    layout: &Layout,
    current: &mut Scores<S>,
    previous: &Scores<S>,
    values: &Values<S>,
    traceback: &mut [u8],
) -> Result<(), Error> {
    if diagonal_index == 0 {
        set_cell(
            current,
            0,
            (S::ZERO, None),
            (S::MIN, None),
            (S::MIN, None),
            traceback,
            diagonal.offset,
        );
        return Ok(());
    }

    let Some(previous_diagonal) = diagonal_index
        .checked_sub(1)
        .map(|index| layout.diagonals[index])
    else {
        return Ok(());
    };

    if diagonal_index < layout.columns {
        let (insertion, predecessor) = choose_scalar([
            (
                previous.insertion[0],
                values.gap_extend_score,
                State::Insertion,
            ),
            (previous.aligned[0], values.gap_open_score, State::Aligned),
            (previous.deletion[0], values.gap_open_score, State::Deletion),
        ])?;
        set_cell(
            current,
            0,
            (S::MIN, None),
            (insertion, predecessor),
            (S::MIN, None),
            traceback,
            diagonal.offset,
        );
    }

    if diagonal_index < layout.rows {
        let score_index = diagonal_index - diagonal.first_row;
        let previous_index = diagonal_index - 1 - previous_diagonal.first_row;
        let (deletion, predecessor) = choose_scalar([
            (
                previous.deletion[previous_index],
                values.gap_extend_score,
                State::Deletion,
            ),
            (
                previous.aligned[previous_index],
                values.gap_open_score,
                State::Aligned,
            ),
            (
                previous.insertion[previous_index],
                values.gap_open_score,
                State::Insertion,
            ),
        ])?;
        set_cell(
            current,
            score_index,
            (S::MIN, None),
            (S::MIN, None),
            (deletion, predecessor),
            traceback,
            diagonal.offset + score_index,
        );
    }

    Ok(())
}

/// Initializes the reachable local reset cells on one anti-diagonal.
fn initialize_local_boundary<S: Narrow>(
    diagonal_index: usize,
    diagonal: Diagonal,
    rows: usize,
    columns: usize,
    current: &mut Scores<S>,
    traceback: &mut [u8],
) {
    if diagonal_index < columns {
        set_cell(
            current,
            0,
            (S::ZERO, None),
            (S::MIN, None),
            (S::MIN, None),
            traceback,
            diagonal.offset,
        );
    }

    if diagonal_index != 0 && diagonal_index < rows {
        let score_index = diagonal_index - diagonal.first_row;
        set_cell(
            current,
            score_index,
            (S::ZERO, None),
            (S::MIN, None),
            (S::MIN, None),
            traceback,
            diagonal.offset + score_index,
        );
    }
}

/// Returns the score-index interval containing interior cells on one diagonal.
fn interior_range(diagonal_index: usize, diagonal: Diagonal) -> Option<(usize, usize)> {
    if diagonal_index < 2 {
        return None;
    }

    let first_row = diagonal.first_row.max(1);
    let last_row = (diagonal.first_row + diagonal.len - 1).min(diagonal_index - 1);

    (first_row <= last_row).then_some((
        first_row - diagonal.first_row,
        last_row + 1 - diagonal.first_row,
    ))
}

/// Updates the local endpoint using score, row, and column tie order.
fn update_local_endpoint(endpoint: &mut Option<Endpoint>, score: Score, row: usize, column: usize) {
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

/// Fills all anti-diagonals and returns their packed traceback matrix.
fn compute<K: Kernel>(
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
    mode: Mode,
) -> Result<(TracebackMatrix, Option<Endpoint>), Error> {
    engine::validate_lengths(reference.len(), query.len())?;
    let (rows, columns) = engine::matrix_dimensions(reference.len(), query.len())?;
    let values = Values::<K::Score>::try_new(scoring)?;
    let traceback = TracebackMatrix::try_new(Layout::try_new(rows, columns)?)?;
    let maximum_diagonal_length = rows.min(columns);

    if K::LANES == 0 {
        return Err(Error::MatrixSizeOverflow);
    }

    let mut current = Scores::try_new(maximum_diagonal_length)?;
    let mut previous = Scores::try_new(maximum_diagonal_length)?;
    let mut before_previous = Scores::try_new(maximum_diagonal_length)?;
    let mut predecessors = PredecessorScratch::try_new(K::LANES)?;
    let mut traceback = traceback;
    let mut endpoint = None;
    let diagonal_count = traceback.layout.diagonals.len();

    for diagonal_index in 0..diagonal_count {
        let diagonal = traceback.layout.diagonals[diagonal_index];
        current.clear();

        match mode {
            Mode::Global => initialize_global_boundary(
                diagonal_index,
                diagonal,
                &traceback.layout,
                &mut current,
                &previous,
                &values,
                &mut traceback.entries,
            )?,
            Mode::Local => initialize_local_boundary(
                diagonal_index,
                diagonal,
                rows,
                columns,
                &mut current,
                &mut traceback.entries,
            ),
        }

        if let Some((interior_start, interior_end)) = interior_range(diagonal_index, diagonal) {
            let previous_diagonal = traceback.layout.diagonals[diagonal_index - 1];
            let before_previous_diagonal = traceback.layout.diagonals[diagonal_index - 2];
            let mut current_index = interior_start;

            while K::LANES <= interior_end - current_index {
                let row = diagonal.first_row + current_index;
                let column = diagonal_index - row;
                let left_index = row - previous_diagonal.first_row;
                let up_index = row - 1 - previous_diagonal.first_row;
                let diagonal_score_index = row - 1 - before_previous_diagonal.first_row;

                fill_vector::<K>(
                    &mut current,
                    &previous,
                    &before_previous,
                    current_index,
                    left_index,
                    up_index,
                    diagonal_score_index,
                    row,
                    column,
                    reference,
                    query,
                    &values,
                    mode,
                    &mut predecessors,
                    &mut traceback.entries,
                    diagonal.offset + current_index,
                );
                current_index += K::LANES;
            }

            while current_index < interior_end {
                let row = diagonal.first_row + current_index;
                let column = diagonal_index - row;
                let left_index = row - previous_diagonal.first_row;
                let up_index = row - 1 - previous_diagonal.first_row;
                let diagonal_score_index = row - 1 - before_previous_diagonal.first_row;

                fill_scalar(
                    &mut current,
                    &previous,
                    &before_previous,
                    current_index,
                    left_index,
                    up_index,
                    diagonal_score_index,
                    row,
                    column,
                    reference,
                    query,
                    &values,
                    mode,
                    &mut traceback.entries,
                    diagonal.offset + current_index,
                )?;
                current_index += 1;
            }

            if mode == Mode::Local {
                for score_index in interior_start..interior_end {
                    let row = diagonal.first_row + score_index;
                    let column = diagonal_index - row;
                    update_local_endpoint(
                        &mut endpoint,
                        current.aligned[score_index].into_score(),
                        row,
                        column,
                    );
                }
            }
        }

        if mode == Mode::Global && diagonal_index + 1 == diagonal_count {
            let score_index = reference.len() - diagonal.first_row;
            let (score, state) = select_score([
                (current.aligned[score_index], State::Aligned),
                (current.deletion[score_index], State::Deletion),
                (current.insertion[score_index], State::Insertion),
            ]);

            if let Some(state) = state {
                endpoint = Some(Endpoint {
                    score: score.into_score(),
                    state,
                    row: reference.len(),
                    column: query.len(),
                });
            }
        }

        std::mem::swap(&mut before_previous, &mut previous);
        std::mem::swap(&mut previous, &mut current);
    }

    Ok((traceback, endpoint))
}

/// Walks packed predecessor fields backward to recover a unit-operation path.
fn traceback(
    matrix: &TracebackMatrix,
    reference: &[u8],
    query: &[u8],
    mut state: State,
    mut row: usize,
    mut column: usize,
) -> Result<(Vec<OperationKind>, usize, usize), Error> {
    let maximum_path_length = reference
        .len()
        .checked_add(query.len())
        .ok_or(Error::MatrixSizeOverflow)?;
    let mut path = Vec::new();
    path.try_reserve_exact(maximum_path_length)
        .map_err(|source| Error::MatrixAllocation { source })?;

    loop {
        let packed = matrix.entries[matrix.layout.index(row, column)];
        let Some(predecessor) = unpack_predecessor(packed, state) else {
            break;
        };

        match state {
            State::Aligned => {
                path.push(if reference[row - 1] == query[column - 1] {
                    OperationKind::SequenceMatch
                } else {
                    OperationKind::SequenceMismatch
                });
                row -= 1;
                column -= 1;
            }
            State::Insertion => {
                path.push(OperationKind::Insertion);
                column -= 1;
            }
            State::Deletion => {
                path.push(OperationKind::Deletion);
                row -= 1;
            }
        }

        state = predecessor;
    }

    path.reverse();
    Ok((path, row, column))
}

/// Computes a checked global alignment through a generic wavefront kernel.
pub(super) fn global<K: Kernel>(
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Outcome, Error> {
    if reference.is_empty() && query.is_empty() {
        return Err(Error::EmptyGlobal);
    }
    if !kernel_is_eligible::<K::Score>(reference.len(), query.len(), scoring) {
        return engine::global(reference, query, scoring);
    }

    let (matrix, endpoint) = compute::<K>(reference, query, scoring, Mode::Global)?;
    let Some(endpoint) = endpoint else {
        // SAFETY: A non-empty global matrix always has a reachable endpoint.
        unreachable!()
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
pub(super) fn local<K: Kernel>(
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    if reference.is_empty() || query.is_empty() {
        return Ok(None);
    }
    if !kernel_is_eligible::<K::Score>(reference.len(), query.len(), scoring) {
        return engine::local(reference, query, scoring);
    }

    let (matrix, endpoint) = compute::<K>(reference, query, scoring, Mode::Local)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::Score;
    use crate::algorithm::engine;

    /// A scalar-array kernel that exercises eight independent vector lanes.
    #[derive(Clone, Copy)]
    struct TestI16;

    // SAFETY: Every operation implements the documented lane-wise behavior
    // for exactly `LANES` elements.
    unsafe impl Kernel for TestI16 {
        type Score = i16;
        type Vector = [i16; 8];

        const LANES: usize = 8;

        unsafe fn load(pointer: *const Self::Score) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for (lane, value) in values.iter_mut().enumerate() {
                // SAFETY: The Kernel caller provides `LANES` readable scores.
                *value = unsafe { pointer.add(lane).read() };
            }
            values
        }

        unsafe fn store(pointer: *mut Self::Score, value: Self::Vector) {
            for (lane, value) in value.into_iter().enumerate() {
                // SAFETY: The Kernel caller provides `LANES` writable scores.
                unsafe { pointer.add(lane).write(value) };
            }
        }

        fn splat(value: Self::Score) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for element in &mut values {
                *element = value;
            }
            values
        }

        fn add(left: Self::Vector, right: Self::Vector) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for lane in 0..Self::LANES {
                values[lane] = left[lane].wrapping_add(right[lane]);
            }
            values
        }

        fn greater_than(left: Self::Vector, right: Self::Vector) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for lane in 0..Self::LANES {
                values[lane] = if left[lane] > right[lane] { -1 } else { 0 };
            }
            values
        }

        fn equal(left: Self::Vector, right: Self::Vector) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for lane in 0..Self::LANES {
                values[lane] = if left[lane] == right[lane] { -1 } else { 0 };
            }
            values
        }

        fn select(
            mask: Self::Vector,
            if_true: Self::Vector,
            if_false: Self::Vector,
        ) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for lane in 0..Self::LANES {
                values[lane] = if mask[lane] != 0 {
                    if_true[lane]
                } else {
                    if_false[lane]
                };
            }
            values
        }

        unsafe fn substitution(
            reference: *const u8,
            query_end: *const u8,
            match_score: Self::Score,
            mismatch_score: Self::Score,
        ) -> Self::Vector {
            let mut values = [0; Self::LANES];
            for (lane, value) in values.iter_mut().enumerate() {
                // SAFETY: The Kernel caller supplies `LANES` reference bytes
                // and a query end pointer with `LANES` preceding bytes.
                let reference_byte = unsafe { reference.add(lane).read() };
                // SAFETY: The Kernel caller supplies `LANES` preceding query
                // bytes ending at `query_end`.
                let query_byte = unsafe { query_end.sub(lane + 1).read() };
                *value = if reference_byte == query_byte {
                    match_score
                } else {
                    mismatch_score
                };
            }
            values
        }
    }

    /// Returns every binary byte sequence through `maximum_length`.
    fn binary_sequences(maximum_length: usize) -> Vec<Vec<u8>> {
        let mut sequences = Vec::new();

        for length in 0..=maximum_length {
            for bits in 0..(1usize << length) {
                sequences.push((0..length).map(|bit| ((bits >> bit) & 1) as u8).collect());
            }
        }

        sequences
    }

    /// Returns the scoring fixtures shared by exhaustive differential tests.
    fn exhaustive_scorings() -> Result<[Scoring; 4], Error> {
        Ok([
            Scoring::try_new(2, -3, -2, -1)?,
            Scoring::try_new(1, -5, -1, -1)?,
            Scoring::try_new(1, 0, 0, 0)?,
            Scoring::try_new(1, -1, 0, -1)?,
        ])
    }

    /// Compares the generic global kernel result with the scalar reference.
    fn differential_global<K: Kernel>(maximum_length: usize) -> Result<(), Error> {
        let sequences = binary_sequences(maximum_length);

        for scoring in exhaustive_scorings()? {
            for reference in &sequences {
                for query in &sequences {
                    let actual = global::<K>(reference, query, scoring);
                    let expected = engine::global(reference, query, scoring);

                    match (actual, expected) {
                        (Ok(actual), Ok(expected)) => assert_eq!(
                            actual, expected,
                            "global mismatch; reference={reference:?} query={query:?} \
                             scoring={scoring:?}"
                        ),
                        (Err(Error::EmptyGlobal), Err(Error::EmptyGlobal)) => {}
                        (actual, expected) => panic!(
                            "global result mismatch; reference={reference:?} query={query:?} \
                             scoring={scoring:?} actual={actual:?} expected={expected:?}"
                        ),
                    }
                }
            }
        }

        Ok(())
    }

    /// Compares the generic local kernel result with the scalar reference.
    fn differential_local<K: Kernel>(maximum_length: usize) -> Result<(), Error> {
        let sequences = binary_sequences(maximum_length);

        for scoring in exhaustive_scorings()? {
            for reference in &sequences {
                for query in &sequences {
                    let actual = local::<K>(reference, query, scoring);
                    let expected = engine::local(reference, query, scoring);

                    match (actual, expected) {
                        (Ok(actual), Ok(expected)) => assert_eq!(
                            actual, expected,
                            "local mismatch; reference={reference:?} query={query:?} \
                             scoring={scoring:?}"
                        ),
                        (actual, expected) => panic!(
                            "local result mismatch; reference={reference:?} query={query:?} \
                             scoring={scoring:?} actual={actual:?} expected={expected:?}"
                        ),
                    }
                }
            }
        }

        Ok(())
    }

    #[derive(Clone, Copy)]
    struct DifferentialCase<'a> {
        name: &'static str,
        reference: &'a [u8],
        query: &'a [u8],
        full_vectors: usize,
        tail_cells: usize,
    }

    /// Returns the deterministic `TestI16` cases used to reach vector paths.
    fn vector_differential_cases() -> [DifferentialCase<'static>; 3] {
        [
            DifferentialCase {
                name: "vector-only",
                reference: b"ACGTACGT",
                query: b"ACGTACGT",
                full_vectors: 1,
                tail_cells: 0,
            },
            DifferentialCase {
                name: "vector-plus-tail",
                reference: b"CGTACGTAC",
                query: b"CGTACGTAC",
                full_vectors: 1,
                tail_cells: 1,
            },
            DifferentialCase {
                name: "two-full-vectors",
                reference: b"TGCATGCATGCATGCA",
                query: b"TGCATGCATGCATGCA",
                full_vectors: 2,
                tail_cells: 0,
            },
        ]
    }

    /// Returns the `TestI16` full-vector and scalar-tail counts implied by the
    /// input lengths.
    fn vector_path_counts(reference_len: usize, query_len: usize) -> (usize, usize) {
        let maximum_interior = reference_len.min(query_len);
        (
            maximum_interior / TestI16::LANES,
            maximum_interior % TestI16::LANES,
        )
    }

    fn assert_vector_case_shape(case: DifferentialCase<'_>) {
        assert_eq!(
            vector_path_counts(case.reference.len(), case.query.len()),
            (case.full_vectors, case.tail_cells),
            "{} lengths should imply {} full vectors and {} tail cells",
            case.name,
            case.full_vectors,
            case.tail_cells,
        );
    }

    fn differential_global_cases<K: Kernel>(cases: &[DifferentialCase<'_>]) -> Result<(), Error> {
        for case in cases {
            assert_vector_case_shape(*case);
        }

        for scoring in exhaustive_scorings()? {
            for case in cases {
                assert_eq!(
                    global::<K>(case.reference, case.query, scoring)?,
                    engine::global(case.reference, case.query, scoring)?,
                    "global mismatch; case={} scoring={scoring:?}",
                    case.name,
                );
            }
        }

        Ok(())
    }

    fn differential_local_cases<K: Kernel>(cases: &[DifferentialCase<'_>]) -> Result<(), Error> {
        for case in cases {
            assert_vector_case_shape(*case);
        }

        for scoring in exhaustive_scorings()? {
            for case in cases {
                assert_eq!(
                    local::<K>(case.reference, case.query, scoring)?,
                    engine::local(case.reference, case.query, scoring)?,
                    "local mismatch; case={} scoring={scoring:?}",
                    case.name,
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_vector_differential_cases_cover_test_i16_paths() {
        assert_eq!(vector_path_counts(5, 5), (0, 5));

        for case in vector_differential_cases() {
            assert_vector_case_shape(case);
        }
    }

    #[test]
    fn test_kernel_matches_scalar_for_global_vector_paths() -> Result<(), Error> {
        differential_global_cases::<TestI16>(&vector_differential_cases())
    }

    #[test]
    fn test_kernel_matches_scalar_for_local_vector_paths() -> Result<(), Error> {
        differential_local_cases::<TestI16>(&vector_differential_cases())
    }

    #[test]
    fn test_kernel_falls_back_outside_the_narrow_range() -> Result<(), Error> {
        let cases = [
            (
                b"A".as_slice(),
                b"A".as_slice(),
                Scoring::try_new(32_768, -1, -1, -1)?,
            ),
            (
                b"AAAAAAAA".as_slice(),
                b"AAAAAAAA".as_slice(),
                Scoring::try_new(5_000, -1, -1, -1)?,
            ),
        ];

        for (reference, query, scoring) in cases {
            assert_eq!(
                global::<TestI16>(reference, query, scoring)?,
                engine::global(reference, query, scoring)?
            );
            assert_eq!(
                local::<TestI16>(reference, query, scoring)?,
                engine::local(reference, query, scoring)?
            );
        }

        Ok(())
    }

    #[test]
    fn test_kernel_matches_scalar_for_short_global_inputs() -> Result<(), Error> {
        differential_global::<TestI16>(5)
    }

    #[test]
    fn test_kernel_matches_scalar_for_short_local_inputs() -> Result<(), Error> {
        differential_local::<TestI16>(5)
    }

    #[test]
    fn selects_the_narrowest_proven_lane() -> Result<(), Error> {
        let scoring = Scoring::try_new(2, -3, -2, -1)?;
        assert_eq!(lane_width(100, 100, scoring), LaneWidth::I16);
        assert_eq!(lane_width(10_000, 10_000, scoring), LaneWidth::I32);

        let extreme = Scoring::try_new(1, -1, Score::MIN, -1)?;
        assert_eq!(lane_width(1, 1, extreme), LaneWidth::Scalar);
        Ok(())
    }

    #[test]
    fn diagonal_layout_round_trips_coordinates() -> Result<(), Error> {
        let layout = Layout::try_new(4, 3)?;

        for row in 0..4 {
            for column in 0..3 {
                let index = layout.index(row, column);
                assert_eq!(layout.coordinate(index), (row, column));
            }
        }

        Ok(())
    }

    #[test]
    fn diagonal_layout_reports_dimension_overflow() {
        assert!(matches!(
            Layout::try_new(usize::MAX, 2),
            Err(Error::MatrixSizeOverflow)
        ));
    }

    #[test]
    fn diagonal_layout_rejects_zero_sized_dimensions() {
        assert!(matches!(
            Layout::try_new(0, 1),
            Err(Error::MatrixSizeOverflow)
        ));
        assert!(matches!(
            Layout::try_new(1, 0),
            Err(Error::MatrixSizeOverflow)
        ));
    }

    #[test]
    fn diagonal_layout_reports_allocation_failure() {
        let error = match Layout::try_new(1, usize::MAX - 1) {
            Err(error) => error,
            Ok(_) => panic!("expected allocation failure"),
        };

        assert!(matches!(error, Error::MatrixAllocation { .. }));
        assert_eq!(error.to_string(), "failed to allocate alignment storage");
    }
}
