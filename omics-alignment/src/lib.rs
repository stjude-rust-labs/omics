//! Foundational representations of sequence alignments in the Rust omics
//! ecosystem.
//!
//! This crate provides validated pairwise alignments, lossless per-operation
//! traversal, and neutral aligned blocks. It has no dependency on SAM
//! records, BAM, PAF, or chainfile libraries.
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
//! assert_eq!(alignment.aligned_blocks().count(), 2);
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
//! # Stranded starts
//!
//! Both the reference and query starting coordinates carry a strand. The
//! reference is typically positive-stranded, but the query may be on the
//! positive strand (a forward read) or the negative strand (a
//! reverse-complement read). [`Alignment::try_new`] treats both strands
//! symmetrically; every operation that consumes an axis moves the
//! corresponding coordinate "forward" for that axis, which means incrementing
//! the position on the positive strand and decrementing it on the negative
//! strand. A negative-strand query traverses from a higher position toward a
//! lower one, as in the quick-start example above where `query:-:5` ends at
//! `query:-:0` after consuming five query bases.
//!
//! # Eager validation
//!
//! [`Alignment::try_new`] replays every operation against the starting
//! coordinates before returning. If any operation would move a coordinate out
//! of representable bounds (position underflow on the negative strand, or
//! overflow on the positive strand), construction returns
//! [`alignment::Error::OutOfBounds`]. Once an [`Alignment`] exists,
//! [`Alignment::steps`] and [`Alignment::aligned_blocks`] are infallible.
//!
//! # `steps()` versus `aligned_blocks()`
//!
//! [`Alignment::steps`] yields one [`Step`] per CIGAR operation in order,
//! preserving every kind including hard clips and padding (which carry no
//! interval). [`Alignment::aligned_blocks`] skips all non-aligned operations
//! and merges consecutive aligned operations (`M`, `=`, `X`) into single
//! [`AlignedBlock`] values; a single insertion or deletion between two aligned
//! runs produces two separate blocks. Prefer [`Alignment::steps`] when you
//! need per-operation precision; prefer [`Alignment::aligned_blocks`] when you
//! need contiguous matched regions.
//!
//! # Alignments
//!
//! The [`alignment`] module provides [`Alignment`], which eagerly validates
//! that a starting reference coordinate, a starting query coordinate, and a
//! [`cigar::Cigar`] together describe a traversal that never moves out of
//! representable coordinate bounds. Once constructed, [`Alignment::steps`]
//! losslessly yields one [`Step`] per operation, and
//! [`Alignment::aligned_blocks`] coalesces consecutive aligned steps into
//! neutral [`AlignedBlock`]s.
//!
//! # Aligned blocks
//!
//! The [`block`] module provides [`block::AlignedBlock`], a neutral primitive
//! pairing equal-length reference and query interbase intervals.
//!
//! # CIGAR model
//!
//! The [`cigar`] module provides the checked CIGAR model, including
//! [`cigar::OperationKind`], [`cigar::Operation`], and [`cigar::Cigar`].

pub use alignment::Alignment;
pub use block::AlignedBlock;
pub use step::Step;

/// Validated pairwise alignments with lossless per-operation traversal.
pub mod alignment;

/// Neutral aligned blocks pairing equal-length reference and query intervals.
pub mod block;

/// The CIGAR model for representing pairwise sequence alignments.
pub mod cigar;

/// Individual lossless alignment steps.
pub mod step;
