//! Arithmetic operations.

/// Safe addition.
pub trait CheckedAdd<T>: Sized {
    /// The output type.
    type Output;

    /// Adds two items.
    ///
    /// - If the addition occurs succesfully, then [`Some<Self>`] is returned.
    /// - If the addition would overflow, [`None`] is returned.
    fn checked_add(&self, rhs: T) -> Option<Self::Output>;
}

/// Safe subtraction.
pub trait CheckedSub<T>: Sized {
    /// The output type.
    type Output;

    /// Subtracts two items.
    ///
    /// - If the subtraction occurs successfully, then [`Some<Self>`] is
    ///   returned.
    /// - If the subtraction would overflow, [`None`] is returned.
    fn checked_sub(&self, rhs: T) -> Option<Self::Output>;
}
