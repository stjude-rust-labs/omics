//! Foundational representations of sequence alignments in the Rust omics
//! ecosystem.
//!
//! This crate provides validated pairwise alignments and lossless
//! per-operation traversal. It has no dependency on SAM records, BAM, PAF, or
//! chainfile libraries.
//!
//! # Quick start
//!
//! ```
//! use omics_alignment::Alignment;
//! use omics_alignment::cigar::Cigar;
//! use omics_coordinate::interbase::Coordinate;
//!
//! let reference_start = "ref:+:0".parse::<Coordinate>()?;
//! let query_start = "query:-:5".parse::<Coordinate>()?;
//! let cigar = "3M1I1M".parse::<Cigar>()?;
//!
//! let alignment = Alignment::try_new(reference_start, query_start, cigar)?;
//! assert_eq!(alignment.steps().count(), 3);
//! assert_eq!(
//!     alignment.steps().filter(|step| step.is_aligned()).count(),
//!     2
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Operation consumption
//!
//! Each CIGAR operation advances zero, one, or both coordinate axes. Aligned
//! operations (`M`, `=`, `X`) advance both the reference and the query by the
//! same length. Insertions (`I`) and soft clips (`S`) advance only the query.
//! Deletions (`D`) and reference skips (`N`) advance only the reference. Hard
//! clips (`H`) and padding (`P`) advance neither axis; they yield a [`Step`]
//! with no reference or query interval.
//!
//! See [`cigar::OperationKind`] for the full per-variant documentation.
//!
//! # Filtering aligned steps
//!
//! [`Alignment::steps`] yields one [`Step`] per CIGAR operation in order,
//! preserving every kind including hard clips and padding (which carry no
//! interval). Callers can select the aligned operations (`M`, `=`, `X`) with
//! `alignment.steps().filter(|step| step.is_aligned())`. Filtering preserves
//! each CIGAR operation boundary rather than imposing a coalescing policy.
//!
//! # Alignments
//!
//! The [`alignment`] module provides [`Alignment`], which eagerly validates
//! that a starting reference coordinate, a starting query coordinate, and a
//! [`cigar::Cigar`] together describe a traversal that never moves out of
//! representable coordinate bounds. Once constructed, [`Alignment::steps`]
//! losslessly yields one [`Step`] per operation.
//!
//! # CIGAR model
//!
//! The [`cigar`] module provides the checked CIGAR model, including
//! [`cigar::OperationKind`], [`cigar::Operation`], and [`cigar::Cigar`].

pub use alignment::Alignment;
pub use step::Step;

/// Validated pairwise alignments with lossless per-operation traversal.
pub mod alignment;

/// The CIGAR model for representing pairwise sequence alignments.
pub mod cigar;

/// Individual lossless alignment steps.
pub mod step;
