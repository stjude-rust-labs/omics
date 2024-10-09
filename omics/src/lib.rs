//! Foundational representations for Rust omics ecosystem.
//!
//! See each individual module's documentation for more information.

#[cfg(feature = "coordinate")]
#[doc(inline)]
pub use omics_coordinate as coordinate;
#[cfg(feature = "core")]
#[doc(inline)]
pub use omics_core as core;
#[cfg(feature = "molecule")]
#[doc(inline)]
pub use omics_molecule as molecule;
#[cfg(feature = "variation")]
#[doc(inline)]
pub use omics_variation as variation;
