//! Foundational representations for Rust omics ecosystem.
//!
//! See each individual module's documentation for more information.
//!
//! # Optional components
//!
//! Each component is gated behind a Cargo feature of the same name and can be
//! enabled independently.
//!
//! ## Alignment
//!
//! Enable the `alignment` feature to access validated pairwise sequence
//! traversal. The feature automatically enables `coordinate`.
//!
//! ```
//! # #[cfg(feature = "alignment")]
//! # {
//! use omics::alignment::Alignment;
//! use omics::alignment::cigar::Cigar;
//! use omics::coordinate::interbase::Coordinate;
//!
//! let reference_start = "ref:+:0".parse::<Coordinate>().unwrap();
//! let query_start = "query:-:5".parse::<Coordinate>().unwrap();
//! let cigar = "3M1I1M".parse::<Cigar>().unwrap();
//! let alignment = Alignment::try_new(reference_start, query_start, cigar).unwrap();
//! assert_eq!(alignment.steps().count(), 3);
//! assert_eq!(alignment.aligned_blocks().count(), 2);
//! # }
//! ```

#[cfg(feature = "alignment")]
#[doc(inline)]
pub use omics_alignment as alignment;
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
