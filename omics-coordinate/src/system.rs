//! Genomic coordinate systems.

pub mod base;
mod corpus;
pub mod interbase;

pub use base::Base;
pub use corpus::Corpus;
pub use interbase::Interbase;

/// A trait for coordinate systems.
pub mod r#trait {
    use std::fmt::Debug;
    use std::fmt::Display;

    /// A coordinate system.
    pub trait System: Clone + Default + Debug + Display + PartialOrd {}
}

pub struct System<S: r#trait::System> {
    corpus: Corpus,

    system: S,
}
