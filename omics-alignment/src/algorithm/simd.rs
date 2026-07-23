//! SIMD-dispatched affine-gap alignment for byte slices.

use super::Error;
use super::Outcome;
use super::Scoring;
use super::engine;

/// Shared private wavefront helpers for future SIMD kernels.
mod wavefront;

/// Returns the global alignment for byte slices.
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    engine::global(reference, query, scoring)
}

/// Returns the local alignment for byte slices.
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    engine::local(reference, query, scoring)
}
