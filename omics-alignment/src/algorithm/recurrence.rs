//! The affine-gap recurrence vocabulary shared by both alignment engines.
//!
//! Every matrix cell `(row, column)` holds three states. `row` counts consumed
//! reference symbols and `column` counts consumed query symbols:
//!
//! - [`State::Aligned`] ends by consuming one symbol from each input, so it
//!   extends a state of the diagonal neighbour `(row - 1, column - 1)`.
//! - [`State::Insertion`] ends by consuming one query symbol, so it extends a
//!   state of the left neighbour `(row, column - 1)`.
//! - [`State::Deletion`] ends by consuming one reference symbol, so it extends
//!   a state of the up neighbour `(row - 1, column)`.
//!
//! Each state is the best of three candidates, one per predecessor state.
//! [`State::candidates`] lists them in tie-break order (see "Determinism" in
//! the [`algorithm`](super) module docs): the first candidate with the
//! greatest score wins.
//!
//! The scalar engine (`engine.rs`) deliberately spells this recurrence out by
//! hand so that it remains an independent oracle for the SIMD wavefront
//! engine, which is driven by the [`State::candidates`] table. Both engines
//! share [`State`] and the [`traceback`] walker.

use super::Error;
use crate::cigar::OperationKind;

/// One state in the affine-gap recurrence.
///
/// The discriminants are the two-bit predecessor codes stored by the SIMD
/// wavefront engine; zero is reserved for "no predecessor".
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub(super) enum State {
    /// Ends by consuming one symbol from each input.
    Aligned = 1,
    /// Ends by consuming one query symbol.
    Insertion = 2,
    /// Ends by consuming one reference symbol.
    Deletion = 3,
}

/// The neighbouring cell that a state extends.
#[cfg_attr(
    not(any(
        test,
        all(target_arch = "aarch64", target_os = "macos"),
        all(target_arch = "x86_64", target_os = "linux")
    )),
    expect(
        dead_code,
        reason = "only the SIMD wavefront engine reads the recurrence table"
    )
)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Neighbour {
    /// The cell at `(row - 1, column - 1)`.
    Diagonal,
    /// The cell at `(row, column - 1)`.
    Left,
    /// The cell at `(row - 1, column)`.
    Up,
}

/// The score added to a predecessor's score to form a candidate.
#[cfg_attr(
    not(any(
        test,
        all(target_arch = "aarch64", target_os = "macos"),
        all(target_arch = "x86_64", target_os = "linux")
    )),
    expect(
        dead_code,
        reason = "only the SIMD wavefront engine reads the recurrence table"
    )
)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Delta {
    /// The match or mismatch score of the two consumed symbols.
    Substitution,
    /// The score for opening a gap after a different state.
    GapOpen,
    /// The score for extending a gap in the same state.
    GapExtend,
}

#[cfg_attr(
    not(any(
        test,
        all(target_arch = "aarch64", target_os = "macos"),
        all(target_arch = "x86_64", target_os = "linux")
    )),
    expect(
        dead_code,
        reason = "only the SIMD wavefront engine reads the recurrence table"
    )
)]
impl State {
    /// The global endpoint's candidate states in tie-break order.
    pub(super) const GLOBAL_ENDPOINT_ORDER: [Self; 3] =
        [Self::Aligned, Self::Deletion, Self::Insertion];

    /// Returns the neighbouring cell whose states this state extends.
    #[inline]
    pub(super) const fn neighbour(self) -> Neighbour {
        match self {
            Self::Aligned => Neighbour::Diagonal,
            Self::Insertion => Neighbour::Left,
            Self::Deletion => Neighbour::Up,
        }
    }

    /// Returns this state's predecessor states and their deltas in tie-break
    /// order.
    ///
    /// A gap state prefers extending itself, then opening from the aligned
    /// state, then opening from the other gap state.
    #[inline]
    pub(super) const fn candidates(self) -> [(Self, Delta); 3] {
        match self {
            Self::Aligned => [
                (Self::Aligned, Delta::Substitution),
                (Self::Deletion, Delta::Substitution),
                (Self::Insertion, Delta::Substitution),
            ],
            Self::Insertion => [
                (Self::Insertion, Delta::GapExtend),
                (Self::Aligned, Delta::GapOpen),
                (Self::Deletion, Delta::GapOpen),
            ],
            Self::Deletion => [
                (Self::Deletion, Delta::GapExtend),
                (Self::Aligned, Delta::GapOpen),
                (Self::Insertion, Delta::GapOpen),
            ],
        }
    }

    /// Decodes a two-bit predecessor code; zero means "no predecessor".
    #[inline]
    pub(super) const fn from_code(code: u8) -> Option<Self> {
        match code {
            1 => Some(Self::Aligned),
            2 => Some(Self::Insertion),
            3 => Some(Self::Deletion),
            _ => None,
        }
    }
}

/// A completed matrix that records the chosen predecessor of every state.
pub(super) trait Predecessors {
    /// Returns the predecessor chosen for `state` at `(row, column)`, or
    /// `None` when the state starts the alignment (the origin or a local
    /// reset).
    fn predecessor(&self, row: usize, column: usize, state: State) -> Option<State>;
}

/// Walks predecessors backward from `(row, column)` in `state` and returns
/// the unit-operation path together with its start coordinate.
pub(super) fn traceback<T: Eq, P: Predecessors>(
    matrix: &P,
    reference: &[T],
    query: &[T],
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

    while let Some(predecessor) = matrix.predecessor(row, column, state) {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A synthetic predecessor table keyed by `(row, column, state)`.
    struct Table(Vec<((usize, usize, State), State)>);

    impl Predecessors for Table {
        fn predecessor(&self, row: usize, column: usize, state: State) -> Option<State> {
            self.0
                .iter()
                .find(|(key, _)| *key == (row, column, state))
                .map(|(_, predecessor)| *predecessor)
        }
    }

    #[test]
    fn candidates_follow_the_documented_tie_order() {
        use Delta::*;
        use State::*;

        assert_eq!(
            Aligned.candidates(),
            [
                (Aligned, Substitution),
                (Deletion, Substitution),
                (Insertion, Substitution)
            ]
        );
        assert_eq!(
            Insertion.candidates(),
            [
                (Insertion, GapExtend),
                (Aligned, GapOpen),
                (Deletion, GapOpen)
            ]
        );
        assert_eq!(
            Deletion.candidates(),
            [
                (Deletion, GapExtend),
                (Aligned, GapOpen),
                (Insertion, GapOpen)
            ]
        );
        assert_eq!(State::GLOBAL_ENDPOINT_ORDER, [Aligned, Deletion, Insertion]);
        assert_eq!(Aligned.neighbour(), Neighbour::Diagonal);
        assert_eq!(Insertion.neighbour(), Neighbour::Left);
        assert_eq!(Deletion.neighbour(), Neighbour::Up);
    }

    #[test]
    fn codes_round_trip() {
        for state in [State::Aligned, State::Insertion, State::Deletion] {
            assert_eq!(State::from_code(state as u8), Some(state));
        }
        assert_eq!(State::from_code(0), None);
    }

    #[test]
    fn traceback_stops_without_a_predecessor() -> Result<(), Error> {
        let table = Table(Vec::new());
        let (path, row, column) = traceback(&table, b"AC", b"AG", State::Aligned, 2, 2)?;

        assert!(path.is_empty());
        assert_eq!((row, column), (2, 2));
        Ok(())
    }

    #[test]
    fn traceback_follows_every_state_transition() -> Result<(), Error> {
        // Forward path for reference "ACG" and query "ATTG":
        // (0,0) -A-> (1,1) -I-> (1,2) -I-> (1,3) -D-> (2,3) -A-> (3,4)
        let table = Table(vec![
            ((3, 4, State::Aligned), State::Deletion),
            ((2, 3, State::Deletion), State::Insertion),
            ((1, 3, State::Insertion), State::Insertion),
            ((1, 2, State::Insertion), State::Aligned),
            ((1, 1, State::Aligned), State::Aligned),
        ]);
        // The origin has no entry, so the walk stops there.
        let (path, row, column) = traceback(&table, b"ACG", b"ATTG", State::Aligned, 3, 4)?;

        assert_eq!(
            path,
            [
                OperationKind::SequenceMatch,
                OperationKind::Insertion,
                OperationKind::Insertion,
                OperationKind::Deletion,
                OperationKind::SequenceMatch,
            ]
        );
        assert_eq!((row, column), (0, 0));
        Ok(())
    }

    #[test]
    fn traceback_reports_mismatches_and_local_starts() -> Result<(), Error> {
        let table = Table(vec![
            ((3, 2, State::Aligned), State::Aligned),
            ((2, 1, State::Aligned), State::Aligned),
        ]);
        let (path, row, column) = traceback(&table, b"GAC", b"AT", State::Aligned, 3, 2)?;

        assert_eq!(
            path,
            [
                OperationKind::SequenceMatch,
                OperationKind::SequenceMismatch
            ]
        );
        assert_eq!((row, column), (1, 0));
        Ok(())
    }
}
