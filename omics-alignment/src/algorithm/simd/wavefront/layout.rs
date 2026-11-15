//! The diagonal-major matrix layout and the per-vector index cursor.
//!
//! [`Layout`] records where each anti-diagonal starts and how many cells it
//! holds, and maps `(row, column)` to an index in the concatenated traceback
//! matrix. [`Cursor`] holds the indexes of one vector's cell, its three
//! neighbours, and its two sequence bytes; stepping down a diagonal advances
//! every index by the same amount, so they are computed once per diagonal.

use crate::algorithm::Error;

/// Metadata for one anti-diagonal in diagonal-major storage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Diagonal {
    /// The first linear index stored on this diagonal.
    pub(super) offset: usize,
    /// The first valid row on this diagonal.
    pub(super) first_row: usize,
    /// The number of cells stored on this diagonal.
    pub(super) len: usize,
}

/// A checked mapping between matrix coordinates and diagonal-major indexes.
pub(super) struct Layout {
    /// The matrix row count.
    pub(super) rows: usize,
    /// The matrix column count.
    pub(super) columns: usize,
    /// The total number of matrix cells.
    pub(super) cells: usize,
    /// The metadata for each anti-diagonal.
    pub(super) diagonals: Vec<Diagonal>,
}

impl Layout {
    /// Builds a checked diagonal-major layout for `rows * columns` cells.
    #[inline]
    pub(super) fn try_new(rows: usize, columns: usize) -> Result<Self, Error> {
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
    #[inline]
    pub(super) fn index(&self, row: usize, column: usize) -> usize {
        debug_assert!(row < self.rows);
        debug_assert!(column < self.columns);

        let diagonal = row + column;
        let diagonal = self.diagonals[diagonal];
        diagonal.offset + row - diagonal.first_row
    }

    /// Maps one diagonal-major index back to its matrix coordinate.
    #[cfg(test)]
    pub(super) fn coordinate(&self, index: usize) -> (usize, usize) {
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

/// The indexes of one interior cell, its three neighbours, and its bytes.
///
/// For a vector, the cursor addresses the first lane and every other lane
/// follows at consecutive indexes.
#[derive(Clone, Copy, Debug)]
pub(super) struct Cursor {
    /// The cell's index on the current diagonal.
    pub(super) current: usize,
    /// The left neighbour's index on the previous diagonal.
    pub(super) left: usize,
    /// The up neighbour's index on the previous diagonal.
    pub(super) up: usize,
    /// The diagonal neighbour's index on the diagonal two back.
    pub(super) diagonal: usize,
    /// The index of the cell's byte in the padded reference.
    pub(super) reference: usize,
    /// The index of the cell's byte in the padded reversed query.
    pub(super) query: usize,
    /// The cell's index in the diagonal-major traceback matrix.
    pub(super) traceback: usize,
}

impl Cursor {
    /// Returns the cursor of the first interior cell on one diagonal.
    ///
    /// Moving down an anti-diagonal increments the row and decrements the
    /// column, which moves every index forward by the same amount. The
    /// indexes are therefore computed once per diagonal.
    #[inline]
    pub(super) fn first(layout: &Layout, diagonal_index: usize, interior_start: usize) -> Self {
        let diagonal = layout.diagonals[diagonal_index];
        let previous = layout.diagonals[diagonal_index - 1];
        let two_back = layout.diagonals[diagonal_index - 2];
        let row = diagonal.first_row + interior_start;
        let column = diagonal_index - row;
        let left = row - previous.first_row;

        Self {
            current: interior_start,
            left,
            up: left - 1,
            diagonal: row - 1 - two_back.first_row,
            reference: row - 1,
            query: layout.columns - 1 - column,
            traceback: diagonal.offset + interior_start,
        }
    }

    /// Moves the cursor `step` cells down the diagonal.
    #[inline]
    pub(super) fn advance(self, step: usize) -> Self {
        Self {
            current: self.current + step,
            left: self.left + step,
            up: self.up + step,
            diagonal: self.diagonal + step,
            reference: self.reference + step,
            query: self.query + step,
            traceback: self.traceback + step,
        }
    }
}

/// Returns the score-index interval containing interior cells on one diagonal.
#[inline]
pub(super) fn interior_range(diagonal_index: usize, diagonal: Diagonal) -> Option<(usize, usize)> {
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
