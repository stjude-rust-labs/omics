//! Convenient access to a 0-based [`Coordinate`](crate::Coordinate).

use crate::Contig;
use crate::Strand;
use crate::position::zero::Position;
use crate::system::Zero;

/// A 0-based [`Coordinate`](crate::Coordinate).
pub type Coordinate = crate::Coordinate<Zero>;

impl crate::r#trait::Coordinate<Zero> for Coordinate {}

impl Coordinate {
    /// Creates a new [`Coordinate`] with a
    /// [`Value::LowerBound`](crate::position::Value::LowerBound).
    ///
    /// Note that a lower bound position can only sit on the
    /// [`Strand::Negative`], so this hardcodes the strand to
    /// [`Strand::Negative`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_coordinate::Coordinate;
    /// use omics_coordinate::Strand;
    /// use omics_coordinate::system::Zero;
    ///
    /// let coordinate = Coordinate::<Zero>::lower_bound("seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn lower_bound<C: TryInto<Contig>>(contig: C) -> Self
    where
        <C as TryInto<Contig>>::Error: std::error::Error,
    {
        // SAFETY: the lower bound is always a valid position for a 0-based position,
        // and we have hardcoded the strand to [`Strand::Negative`], so this
        // will always unwrap.
        Self::try_new(contig, Strand::Negative, Position::lower_bound()).unwrap()
    }
}
