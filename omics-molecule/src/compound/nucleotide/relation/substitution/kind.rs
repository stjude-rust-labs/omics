//! A kind of nucleotide substitution.

/// A change in the type of nucleotide compound.
#[derive(Debug, Eq, PartialEq)]
pub enum Kind {
    /// The nucleobase within the reference nucleotide changes to another
    /// nucleobase of the same [compound
    /// kind](crate::compound::nucleotide::Kind).
    Transition,

    /// The nucleobase within the reference nucleotide changes to another
    /// nucleobase of a different [compound
    /// kind](crate::compound::nucleotide::Kind).
    Transversion,
}
