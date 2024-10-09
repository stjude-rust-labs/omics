//! Genomic coordinate systems.

pub mod one;
pub mod zero;

pub use one::One;
pub use zero::Zero;

/// The requirements to be a coordinate system.
pub trait System: Clone + Default + std::fmt::Debug + std::fmt::Display {}
