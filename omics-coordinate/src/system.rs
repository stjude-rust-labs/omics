//! Genomic coordinate systems.

pub mod base;
pub mod interbase;

use std::fmt::Debug;
use std::fmt::Display;

pub use base::Base;
pub use interbase::Interbase;

/// A coordinate system.
pub trait System: Clone + Default + Debug + Display + PartialOrd {}
