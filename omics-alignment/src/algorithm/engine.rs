//! Dynamic-programming and traceback implementation.

use omics_coordinate::position::Number;

use crate::cigar::Axis;
use crate::cigar::Cigar;
use crate::cigar::Operation;
use crate::cigar::OperationKind;

use super::Error;
use super::Outcome;
use super::Score;
use super::Scoring;

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
    /// An unreachable state.
    const UNREACHABLE: Self = Self {
        score: None,
        predecessor: None,
    };

    /// A zero-valued origin or local reset.
    const RESET: Self = Self {
        score: Some(0),
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
    let length = Number::try_from(run_length).map_err(|_| Error::LengthOutOfRange {
        axis: if kind == OperationKind::Insertion {
            Axis::Query
        } else {
            Axis::Reference
        },
        length: run_length,
    })?;
    operations.push(Operation::try_new(kind, length)?);
    Ok(())
}

/// Converts a non-empty unit-operation path into a canonical CIGAR.
fn path_to_cigar(path: &[OperationKind]) -> Result<Cigar, Error> {
    let first = path.first().copied().ok_or(Error::TracebackInvariant)?;
    let mut operations = Vec::new();
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

/// Drives global matrix initialization, fill, endpoint selection, and traceback.
fn compute_global<T: Eq>(reference: &[T], query: &[T], scoring: Scoring) -> Result<Outcome, Error> {
    if reference.is_empty() && query.is_empty() {
        return Err(Error::EmptyGlobal);
    }

    Number::try_from(reference.len()).map_err(|_| Error::LengthOutOfRange {
        axis: Axis::Reference,
        length: reference.len(),
    })?;
    Number::try_from(query.len()).map_err(|_| Error::LengthOutOfRange {
        axis: Axis::Query,
        length: query.len(),
    })?;

    let rows = reference
        .len()
        .checked_add(1)
        .ok_or(Error::MatrixSizeOverflow)?;
    let columns = query
        .len()
        .checked_add(1)
        .ok_or(Error::MatrixSizeOverflow)?;

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

    let score = endpoint.score.ok_or_else(|| {
        debug_assert!(false, "global endpoint is unreachable");
        Error::TracebackInvariant
    })?;
    let mut state = endpoint.predecessor.ok_or_else(|| {
        debug_assert!(false, "global endpoint has no predecessor");
        Error::TracebackInvariant
    })?;

    let mut path: Vec<OperationKind> = Vec::new();
    let mut i = reference.len();
    let mut j = query.len();

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
    let cigar = path_to_cigar(&path)?;
    Ok(Outcome::new(
        score,
        cigar,
        0..reference.len(),
        0..query.len(),
    ))
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
    fn matrix_dimension_overflow_is_reported() {
        assert!(matches!(
            Matrix::try_new(usize::MAX, 2),
            Err(Error::MatrixSizeOverflow)
        ));
    }

    #[test]
    fn matrix_allocation_failure_is_reported() {
        assert!(matches!(
            Matrix::try_new(1, usize::MAX),
            Err(Error::MatrixAllocation { .. })
        ));
    }
}
