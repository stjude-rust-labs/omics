//! Sequence alignment primitives and algorithms.
//!
//! This crate provides foundational types and algorithms for pairwise sequence
//! alignment in the `omics` ecosystem, including:
//!
//! - **Scoring**: substitution matrices ([`score::NucleotideMatrix`],
//!   [`score::Blosum62`]) and affine gap penalty models ([`score::GapPenalty`])
//! - **Encoding**: integer-encoded sequence representations
//!   ([`sequence::EncodedSequence`]) with conversions from `omics-molecule`
//!   DNA, RNA, and protein types
//! - **Results**: CIGAR string representation ([`cigar::Cigar`]) and alignment
//!   result type ([`cigar::Alignment`])
//! - **Algorithms**: scalar Needleman-Wunsch ([`needleman_wunsch`]) for global
//!   alignment and Smith-Waterman ([`smith_waterman`]) for local alignment,
//!   both with affine gap penalties and full traceback

pub mod cigar;
pub mod needleman_wunsch;
pub mod score;
pub mod sequence;
pub mod smith_waterman;
