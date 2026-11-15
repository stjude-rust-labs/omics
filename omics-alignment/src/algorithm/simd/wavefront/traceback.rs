//! The packed traceback byte and the diagonal-major traceback matrix.
//!
//! Each cell stores one byte holding three two-bit predecessor fields:
//! aligned in bits 0-1, insertion in bits 2-3, and deletion in bits 4-5. A
//! field holds the predecessor [`State`] code, or zero when the state starts
//! the alignment. [`TracebackMatrix`] implements `Predecessors` so the
//! shared traceback walker can read it.

use super::layout::Layout;
use crate::algorithm::Error;
use crate::algorithm::recurrence::Predecessors;
use crate::algorithm::recurrence::State;

/// Returns the bit offset of `state`'s two-bit field in a packed traceback
/// byte.
#[inline]
pub(super) const fn field_shift(state: State) -> u8 {
    match state {
        State::Aligned => 0,
        State::Insertion => 2,
        State::Deletion => 4,
    }
}

/// Returns `predecessor`'s code shifted into `destination`'s field.
#[inline]
pub(super) const fn field_code(predecessor: State, destination: State) -> u8 {
    (predecessor as u8) << field_shift(destination)
}

/// Packs the three predecessor state codes for one matrix cell.
#[inline]
pub(super) fn pack_predecessors(
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
        packed |= code << field_shift(state);
    }

    packed
}

/// Decodes one packed predecessor field.
#[inline]
pub(super) fn unpack_predecessor(packed: u8, state: State) -> Option<State> {
    State::from_code((packed >> field_shift(state)) & 0b11)
}

/// Owns diagonal-major traceback fields for one computed matrix.
pub(super) struct TracebackMatrix {
    /// Coordinate layout for every traceback field.
    pub(super) layout: Layout,
    /// Packed two-bit predecessor fields for every matrix cell.
    pub(super) entries: Vec<u8>,
}

impl TracebackMatrix {
    /// Allocates empty predecessor fields for a checked layout.
    #[inline]
    pub(super) fn try_new(layout: Layout) -> Result<Self, Error> {
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(layout.cells)
            .map_err(|source| Error::MatrixAllocation { source })?;
        entries.resize(layout.cells, 0);

        Ok(Self { layout, entries })
    }
}

impl Predecessors for TracebackMatrix {
    #[inline]
    fn predecessor(&self, row: usize, column: usize, state: State) -> Option<State> {
        unpack_predecessor(self.entries[self.layout.index(row, column)], state)
    }
}
