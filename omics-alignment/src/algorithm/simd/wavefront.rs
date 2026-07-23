#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "Task 2 adds private SIMD wavefront helpers before engine integration"
    )
)]

use crate::algorithm::Error;
use crate::algorithm::Scoring;

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
    let Some(steps) = reference_len.checked_add(query_len) else {
        return LaneWidth::Scalar;
    };
    let magnitude = [
        scoring.match_score().unsigned_abs(),
        scoring.mismatch_score().unsigned_abs(),
        scoring.gap_open_score().unsigned_abs(),
        scoring.gap_extend_score().unsigned_abs(),
    ]
    .into_iter()
    .fold(0u64, u64::max);
    let Ok(steps) = u64::try_from(steps) else {
        return LaneWidth::Scalar;
    };
    let Some(bound) = steps.checked_mul(magnitude) else {
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
            Err(diagonal) => diagonal - 1,
        };
        let metadata = self.diagonals[diagonal];
        let row = metadata.first_row + (index - metadata.offset);

        debug_assert!(row < metadata.first_row + metadata.len);

        (row, diagonal - row)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::Score;

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
