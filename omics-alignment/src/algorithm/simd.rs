//! SIMD-dispatched affine-gap alignment for byte slices.
//!
//! The public API retains the scalar fallback until a target-specific kernel
//! becomes available. The private wavefront engine supports backend
//! differential tests without changing byte-API results.

use super::Error;
use super::Outcome;
use super::Scoring;
use super::engine;

/// Apple Silicon NEON wavefront kernels and dispatch.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
mod aarch64;

/// Shared private wavefront helpers for future SIMD kernels.
mod wavefront;

/// Returns the global alignment for byte slices.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    aarch64::global(reference, query, scoring)
}

/// Returns the global alignment for byte slices.
#[cfg(not(all(target_arch = "aarch64", target_os = "macos")))]
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    engine::global(reference, query, scoring)
}

/// Returns the local alignment for byte slices.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    aarch64::local(reference, query, scoring)
}

/// Returns the local alignment for byte slices.
#[cfg(not(all(target_arch = "aarch64", target_os = "macos")))]
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    engine::local(reference, query, scoring)
}

#[cfg(all(test, target_arch = "aarch64", target_os = "macos"))]
mod tests {
    use super::*;

    #[test]
    fn neon_dispatch_matches_scalar() -> Result<(), Error> {
        let cases = [
            (
                b"ACGTACGT".as_slice(),
                b"ACGTACGT".as_slice(),
                Scoring::try_new(2, -3, -2, -1)?,
            ),
            (
                b"ACGTACGT".as_slice(),
                b"ACGTACGT".as_slice(),
                Scoring::try_new(10_000, -1, -1, -1)?,
            ),
        ];

        for (reference, query, scoring) in cases {
            assert_eq!(
                aarch64::global(reference, query, scoring)?,
                engine::global(reference, query, scoring)?
            );
            assert_eq!(
                aarch64::local(reference, query, scoring)?,
                engine::local(reference, query, scoring)?
            );
        }

        Ok(())
    }
}
