//! Dynamic-programming and traceback implementation.

use omics_coordinate::position::Number;

use super::Error;
use super::Outcome;
use super::Score;
use super::Scoring;
use crate::cigar::Axis;
use crate::cigar::Cigar;
use crate::cigar::Operation;
use crate::cigar::OperationKind;

/// One terminal state in the affine-gap recurrence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum State {
    /// Ends by consuming one symbol from each input.
    Aligned,
    /// Ends by consuming one query symbol.
    Insertion,
    /// Ends by consuming one reference symbol.
    Deletion,
}

/// A score and traceback predecessor for one state.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Entry {
    /// Score, or `None` when the state is unreachable.
    score: Option<Score>,
    /// State at the predecessor coordinate.
    predecessor: Option<State>,
}

impl Entry {
    /// A zero-valued origin or local reset.
    const RESET: Self = Self {
        score: Some(0),
        predecessor: None,
    };
    /// An unreachable state.
    const UNREACHABLE: Self = Self {
        score: None,
        predecessor: None,
    };

    /// Constructs a scored state with a predecessor.
    #[cfg(test)]
    const fn reachable(score: Score, predecessor: State) -> Self {
        Self {
            score: Some(score),
            predecessor: Some(predecessor),
        }
    }
}

/// The three recurrence states at one matrix coordinate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Cell {
    /// State ending with aligned symbols.
    aligned: Entry,
    /// State ending with an insertion.
    insertion: Entry,
    /// State ending with a deletion.
    deletion: Entry,
}

impl Cell {
    /// A cell with no reachable states.
    const UNREACHABLE: Self = Self {
        aligned: Entry::UNREACHABLE,
        insertion: Entry::UNREACHABLE,
        deletion: Entry::UNREACHABLE,
    };

    /// Returns one state entry.
    const fn entry(self, state: State) -> Entry {
        match state {
            State::Aligned => self.aligned,
            State::Insertion => self.insertion,
            State::Deletion => self.deletion,
        }
    }
}

/// Row-major dynamic-programming storage.
struct Matrix {
    /// Number of matrix columns.
    columns: usize,
    /// Row-major cells.
    cells: Vec<Cell>,
}

impl Matrix {
    /// Allocates an unreachable matrix with checked dimensions.
    fn try_new(rows: usize, columns: usize) -> Result<Self, Error> {
        let length = rows.checked_mul(columns).ok_or(Error::MatrixSizeOverflow)?;
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(length)
            .map_err(|source| Error::MatrixAllocation { source })?;
        cells.resize(length, Cell::UNREACHABLE);
        Ok(Self { columns, cells })
    }

    /// Returns a copied cell.
    fn get(&self, row: usize, column: usize) -> Cell {
        self.cells[row * self.columns + column]
    }

    /// Returns a mutable reference to a cell.
    fn get_mut(&mut self, row: usize, column: usize) -> &mut Cell {
        &mut self.cells[row * self.columns + column]
    }
}

/// Adds to a reachable score with overflow checking.
fn checked_add(score: Option<Score>, delta: Score) -> Result<Option<Score>, Error> {
    score
        .map(|s| s.checked_add(delta).ok_or(Error::ScoreOverflow))
        .transpose()
}

/// Chooses the greatest score and preserves candidate-array order on ties.
fn choose<const N: usize>(candidates: [(Option<Score>, State); N]) -> Entry {
    let mut best = Entry::UNREACHABLE;

    for (score, predecessor) in candidates {
        let replace = match (score, best.score) {
            (Some(candidate), Some(current)) => candidate > current,
            (Some(_), None) => true,
            (None, _) => false,
        };
        if replace {
            best = Entry {
                score,
                predecessor: Some(predecessor),
            };
        }
    }

    best
}

/// Appends one checked, coalesced operation to the accumulator.
fn push_operation(
    operations: &mut Vec<Operation>,
    kind: OperationKind,
    run_length: usize,
) -> Result<(), Error> {
    // SAFETY: input lengths are validated against Number before traceback,
    // and a coalesced run cannot exceed the input length for its axis.
    let length = Number::try_from(run_length).unwrap();
    operations.push(Operation::try_new(kind, length)?);
    Ok(())
}

/// Converts a non-empty unit-operation path into a canonical CIGAR.
pub(super) fn path_to_cigar(path: &[OperationKind]) -> Result<Cigar, Error> {
    // SAFETY: global alignment rejects two empty inputs, and local alignment
    // reaches this function only after finding a positive aligned endpoint.
    let first = path.first().copied().unwrap();
    let mut operations = Vec::new();
    operations
        .try_reserve_exact(path.len())
        .map_err(|source| Error::MatrixAllocation { source })?;
    let mut kind = first;
    let mut run_length = 1usize;

    for next in path.iter().copied().skip(1) {
        if next == kind {
            run_length += 1;
        } else {
            push_operation(&mut operations, kind, run_length)?;
            kind = next;
            run_length = 1;
        }
    }
    push_operation(&mut operations, kind, run_length)?;

    Ok(Cigar::try_new(operations)?)
}

/// Computes a checked global alignment.
pub(super) fn global<T: Eq>(
    reference: &[T],
    query: &[T],
    scoring: Scoring,
) -> Result<Outcome, Error> {
    compute_global(reference, query, scoring)
}

/// Computes a checked local alignment.
pub(super) fn local<T: Eq>(
    reference: &[T],
    query: &[T],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    compute_local(reference, query, scoring)
}

/// Validates that traceback run lengths fit the active CIGAR number type.
pub(super) fn validate_lengths(reference: usize, query: usize) -> Result<(), Error> {
    Number::try_from(reference).map_err(|_| Error::LengthOutOfRange {
        axis: Axis::Reference,
        length: reference,
    })?;
    Number::try_from(query).map_err(|_| Error::LengthOutOfRange {
        axis: Axis::Query,
        length: query,
    })?;
    Ok(())
}

/// Returns the checked dynamic-programming dimensions.
pub(super) fn matrix_dimensions(reference: usize, query: usize) -> Result<(usize, usize), Error> {
    let rows = reference.checked_add(1).ok_or(Error::MatrixSizeOverflow)?;
    let columns = query.checked_add(1).ok_or(Error::MatrixSizeOverflow)?;
    Ok((rows, columns))
}

/// Recovers the unit-operation path and start coordinates reached by walking
/// backward to an entry without a predecessor.
fn traceback<T: Eq>(
    matrix: &Matrix,
    reference: &[T],
    query: &[T],
    mut state: State,
    mut i: usize,
    mut j: usize,
) -> Result<(Vec<OperationKind>, usize, usize), Error> {
    let max_path_length = reference
        .len()
        .checked_add(query.len())
        .ok_or(Error::MatrixSizeOverflow)?;
    let mut path = Vec::new();
    path.try_reserve_exact(max_path_length)
        .map_err(|source| Error::MatrixAllocation { source })?;

    loop {
        let entry = matrix.get(i, j).entry(state);
        let Some(predecessor) = entry.predecessor else {
            break;
        };

        match state {
            State::Aligned => {
                path.push(if reference[i - 1] == query[j - 1] {
                    OperationKind::SequenceMatch
                } else {
                    OperationKind::SequenceMismatch
                });
                i -= 1;
                j -= 1;
            }
            State::Insertion => {
                path.push(OperationKind::Insertion);
                j -= 1;
            }
            State::Deletion => {
                path.push(OperationKind::Deletion);
                i -= 1;
            }
        }

        state = predecessor;
    }

    path.reverse();
    Ok((path, i, j))
}

/// Drives global matrix initialization, fill, endpoint selection, and
/// traceback.
fn compute_global<T: Eq>(reference: &[T], query: &[T], scoring: Scoring) -> Result<Outcome, Error> {
    if reference.is_empty() && query.is_empty() {
        return Err(Error::EmptyGlobal);
    }

    validate_lengths(reference.len(), query.len())?;
    let (rows, columns) = matrix_dimensions(reference.len(), query.len())?;

    let mut matrix = Matrix::try_new(rows, columns)?;

    matrix.get_mut(0, 0).aligned = Entry::RESET;

    for j in 1..columns {
        let prev = matrix.get(0, j - 1);
        let insertion = choose([
            (
                checked_add(prev.insertion.score, scoring.gap_extend_score())?,
                State::Insertion,
            ),
            (
                checked_add(prev.aligned.score, scoring.gap_open_score())?,
                State::Aligned,
            ),
            (
                checked_add(prev.deletion.score, scoring.gap_open_score())?,
                State::Deletion,
            ),
        ]);
        matrix.get_mut(0, j).insertion = insertion;
    }

    for i in 1..rows {
        let prev = matrix.get(i - 1, 0);
        let deletion = choose([
            (
                checked_add(prev.deletion.score, scoring.gap_extend_score())?,
                State::Deletion,
            ),
            (
                checked_add(prev.aligned.score, scoring.gap_open_score())?,
                State::Aligned,
            ),
            (
                checked_add(prev.insertion.score, scoring.gap_open_score())?,
                State::Insertion,
            ),
        ]);
        matrix.get_mut(i, 0).deletion = deletion;
    }

    for i in 1..rows {
        for j in 1..columns {
            let subst = scoring.substitution(&reference[i - 1], &query[j - 1]);
            let diag = matrix.get(i - 1, j - 1);
            let left = matrix.get(i, j - 1);
            let up = matrix.get(i - 1, j);

            let new_aligned = choose([
                (checked_add(diag.aligned.score, subst)?, State::Aligned),
                (checked_add(diag.deletion.score, subst)?, State::Deletion),
                (checked_add(diag.insertion.score, subst)?, State::Insertion),
            ]);
            let new_insertion = choose([
                (
                    checked_add(left.insertion.score, scoring.gap_extend_score())?,
                    State::Insertion,
                ),
                (
                    checked_add(left.aligned.score, scoring.gap_open_score())?,
                    State::Aligned,
                ),
                (
                    checked_add(left.deletion.score, scoring.gap_open_score())?,
                    State::Deletion,
                ),
            ]);
            let new_deletion = choose([
                (
                    checked_add(up.deletion.score, scoring.gap_extend_score())?,
                    State::Deletion,
                ),
                (
                    checked_add(up.aligned.score, scoring.gap_open_score())?,
                    State::Aligned,
                ),
                (
                    checked_add(up.insertion.score, scoring.gap_open_score())?,
                    State::Insertion,
                ),
            ]);

            let cell = matrix.get_mut(i, j);
            cell.aligned = new_aligned;
            cell.insertion = new_insertion;
            cell.deletion = new_deletion;
        }
    }

    let end_cell = matrix.get(reference.len(), query.len());
    let endpoint = choose([
        (end_cell.aligned.score, State::Aligned),
        (end_cell.deletion.score, State::Deletion),
        (end_cell.insertion.score, State::Insertion),
    ]);

    // SAFETY: every non-empty global alignment has at least one reachable
    // terminal state after matrix initialization and recurrence evaluation.
    let score = endpoint.score.unwrap();
    // SAFETY: choose assigns a predecessor whenever it selects a reachable
    // terminal state.
    let state = endpoint.predecessor.unwrap();
    let (path, reference_start, query_start) = traceback(
        &matrix,
        reference,
        query,
        state,
        reference.len(),
        query.len(),
    )?;
    debug_assert_eq!(reference_start, 0);
    debug_assert_eq!(query_start, 0);
    let cigar = path_to_cigar(&path)?;
    Ok(Outcome::new(
        score,
        cigar,
        0..reference.len(),
        0..query.len(),
    ))
}

/// Drives local matrix initialization, fill, endpoint selection, and traceback.
fn compute_local<T: Eq>(
    reference: &[T],
    query: &[T],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    if reference.is_empty() || query.is_empty() {
        return Ok(None);
    }

    validate_lengths(reference.len(), query.len())?;
    let (rows, columns) = matrix_dimensions(reference.len(), query.len())?;

    let mut matrix = Matrix::try_new(rows, columns)?;

    for j in 0..columns {
        matrix.get_mut(0, j).aligned = Entry::RESET;
    }
    for i in 1..rows {
        matrix.get_mut(i, 0).aligned = Entry::RESET;
    }

    for i in 1..rows {
        for j in 1..columns {
            let subst = scoring.substitution(&reference[i - 1], &query[j - 1]);
            let diag = matrix.get(i - 1, j - 1);
            let left = matrix.get(i, j - 1);
            let up = matrix.get(i - 1, j);

            let aligned = choose([
                (checked_add(diag.aligned.score, subst)?, State::Aligned),
                (checked_add(diag.deletion.score, subst)?, State::Deletion),
                (checked_add(diag.insertion.score, subst)?, State::Insertion),
            ]);
            let insertion = choose([
                (
                    checked_add(left.insertion.score, scoring.gap_extend_score())?,
                    State::Insertion,
                ),
                (
                    checked_add(left.aligned.score, scoring.gap_open_score())?,
                    State::Aligned,
                ),
                (
                    checked_add(left.deletion.score, scoring.gap_open_score())?,
                    State::Deletion,
                ),
            ]);
            let deletion = choose([
                (
                    checked_add(up.deletion.score, scoring.gap_extend_score())?,
                    State::Deletion,
                ),
                (
                    checked_add(up.aligned.score, scoring.gap_open_score())?,
                    State::Aligned,
                ),
                (
                    checked_add(up.insertion.score, scoring.gap_open_score())?,
                    State::Insertion,
                ),
            ]);

            let cell = matrix.get_mut(i, j);
            cell.aligned = match aligned.score {
                Some(score) if score > 0 => aligned,
                _ => Entry::RESET,
            };
            cell.insertion = insertion;
            cell.deletion = deletion;
        }
    }

    let mut endpoint: Option<(Score, usize, usize)> = None;
    for i in 1..rows {
        for j in 1..columns {
            let score = matrix.get(i, j).aligned.score;
            let Some(score) = score.filter(|score| *score > 0) else {
                continue;
            };

            if endpoint
                .as_ref()
                .map(|(best, ..)| score > *best)
                .unwrap_or(true)
            {
                endpoint = Some((score, i, j));
            }
        }
    }

    let Some((score, reference_end, query_end)) = endpoint else {
        return Ok(None);
    };

    let (path, reference_start, query_start) = traceback(
        &matrix,
        reference,
        query,
        State::Aligned,
        reference_end,
        query_end,
    )?;
    debug_assert!(matches!(
        path.first(),
        Some(OperationKind::SequenceMatch | OperationKind::SequenceMismatch)
    ));
    debug_assert!(matches!(
        path.last(),
        Some(OperationKind::SequenceMatch | OperationKind::SequenceMismatch)
    ));
    let cigar = path_to_cigar(&path)?;

    Ok(Some(Outcome::new(
        score,
        cigar,
        reference_start..reference_end,
        query_start..query_end,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_order_breaks_equal_score_ties() {
        let entry = choose([
            (Some(4), State::Aligned),
            (Some(4), State::Deletion),
            (Some(4), State::Insertion),
        ]);
        assert_eq!(entry, Entry::reachable(4, State::Aligned));
    }

    #[test]
    fn checked_score_addition_reports_overflow() {
        assert!(matches!(
            checked_add(Some(Score::MAX), 1),
            Err(Error::ScoreOverflow)
        ));
        assert_eq!(checked_add(None, 1).unwrap(), None);
    }

    #[test]
    fn traceback_path_coalesces_operations() -> Result<(), Error> {
        let cigar = path_to_cigar(&[
            OperationKind::SequenceMatch,
            OperationKind::SequenceMatch,
            OperationKind::Insertion,
            OperationKind::SequenceMismatch,
            OperationKind::SequenceMismatch,
        ])?;
        assert_eq!(cigar.to_string(), "2=1I2X");
        Ok(())
    }

    #[test]
    fn empty_traceback_path_panics() {
        let result = std::panic::catch_unwind(|| path_to_cigar(&[]));
        assert!(result.is_err());
    }

    #[test]
    fn matrix_dimension_overflow_is_reported() {
        assert!(matches!(
            Matrix::try_new(usize::MAX, 2),
            Err(Error::MatrixSizeOverflow)
        ));
    }

    #[test]
    fn matrix_allocation_failure_is_reported() {
        let error = match Matrix::try_new(1, usize::MAX) {
            Err(error) => error,
            Ok(_) => panic!("expected allocation failure"),
        };
        assert!(matches!(error, Error::MatrixAllocation { .. }));
        assert_eq!(error.to_string(), "failed to allocate alignment storage");
    }
}
