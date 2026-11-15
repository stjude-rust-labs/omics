//! SIMD-dispatched affine-gap alignment for byte slices.
//!
//! [`global`] and [`local`] run a NEON (macOS `aarch64`) or AVX2 (Linux
//! `x86_64`, detected at runtime) wavefront kernel when one is available and
//! the score bound fits a SIMD lane. Otherwise they use the scalar engine.
//! Every path returns the same result as [`super::global`] and
//! [`super::local`].

use std::sync::Once;

use super::Error;
use super::Outcome;
use super::Scoring;
use super::engine;

/// Apple Silicon NEON wavefront kernels and dispatch.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
mod aarch64;

/// Linux x86_64 AVX2 wavefront kernels and dispatch.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
mod x86_64;

mod wavefront;

/// The vector instruction set used for one alignment.
#[derive(Clone, Copy, Debug)]
enum Backend {
    /// NEON, which every macOS `aarch64` processor supports.
    #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
    Neon,
    /// AVX2, detected at runtime on Linux `x86_64`.
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    Avx2(x86_64::Avx2),
    /// No vector backend; the scalar engine computes every alignment.
    #[cfg(not(all(target_arch = "aarch64", target_os = "macos")))]
    Scalar,
}

impl Backend {
    /// Returns the backend's name for structured log fields.
    const fn name(self) -> &'static str {
        match self {
            #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
            Self::Neon => "neon",
            #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
            Self::Avx2(_) => "avx2",
            #[cfg(not(all(target_arch = "aarch64", target_os = "macos")))]
            Self::Scalar => "scalar",
        }
    }

    /// Returns the best backend for the running processor.
    ///
    /// Detection runs on every call, but only the first detection in a process
    /// is logged.
    fn detect() -> Self {
        /// Guards the once-per-process backend log event.
        static LOGGED: Once = Once::new();

        #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
        let backend = Self::Neon;
        #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
        let backend = x86_64::Avx2::detect().map_or(Self::Scalar, Self::Avx2);
        #[cfg(not(any(
            all(target_arch = "aarch64", target_os = "macos"),
            all(target_arch = "x86_64", target_os = "linux")
        )))]
        let backend = Self::Scalar;

        LOGGED.call_once(|| {
            tracing::debug!(backend = backend.name(), "selected SIMD alignment backend");
        });

        backend
    }
}

/// Returns the global alignment for byte slices.
///
/// The result always equals [`super::global`] for the same inputs.
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    match Backend::detect() {
        #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
        Backend::Neon => aarch64::global(reference, query, scoring),
        #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
        Backend::Avx2(avx2) => x86_64::global(avx2, reference, query, scoring),
        #[cfg(not(all(target_arch = "aarch64", target_os = "macos")))]
        Backend::Scalar => {
            tracing::trace!(
                mode = "global",
                reason = "no SIMD backend",
                "falling back to the scalar engine"
            );
            engine::global(reference, query, scoring)
        }
    }
}

/// Returns the local alignment for byte slices.
///
/// The result always equals [`super::local`] for the same inputs.
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    match Backend::detect() {
        #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
        Backend::Neon => aarch64::local(reference, query, scoring),
        #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
        Backend::Avx2(avx2) => x86_64::local(avx2, reference, query, scoring),
        #[cfg(not(all(target_arch = "aarch64", target_os = "macos")))]
        Backend::Scalar => {
            tracing::trace!(
                mode = "local",
                reason = "no SIMD backend",
                "falling back to the scalar engine"
            );
            engine::local(reference, query, scoring)
        }
    }
}

/// Tests that force each native backend and compare it with the scalar engine.
#[cfg(all(
    test,
    any(
        all(target_arch = "aarch64", target_os = "macos"),
        all(target_arch = "x86_64", target_os = "linux")
    )
))]
mod tests;
