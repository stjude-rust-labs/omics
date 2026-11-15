//! Linux `x86_64` AVX2 wavefront kernels and dispatch.
//!
//! `Avx2` is a capability token: only `Avx2::detect` can create one, so
//! holding it proves the processor supports AVX2. `Avx2I16` computes
//! sixteen `i16` lanes per vector and `Avx2I32` eight `i32` lanes; both
//! store the token, which makes their `Kernel` methods sound to call.
//! [`global`] and [`local`] enter an `avx2`-enabled function, pick the
//! narrowest lane width whose score bound fits, and run the shared wavefront
//! engine inlined with that kernel, or fall back to the scalar engine.

use std::arch::x86_64::__m128i;
use std::arch::x86_64::__m256i;
use std::arch::x86_64::_mm_cmpeq_epi8;
use std::arch::x86_64::_mm_loadl_epi64;
use std::arch::x86_64::_mm_loadu_si128;
use std::arch::x86_64::_mm256_add_epi16;
use std::arch::x86_64::_mm256_add_epi32;
use std::arch::x86_64::_mm256_blendv_epi8;
use std::arch::x86_64::_mm256_cmpeq_epi16;
use std::arch::x86_64::_mm256_cmpeq_epi32;
use std::arch::x86_64::_mm256_cmpgt_epi16;
use std::arch::x86_64::_mm256_cmpgt_epi32;
use std::arch::x86_64::_mm256_cvtepi8_epi16;
use std::arch::x86_64::_mm256_cvtepi8_epi32;
use std::arch::x86_64::_mm256_loadu_si256;
use std::arch::x86_64::_mm256_set1_epi16;
use std::arch::x86_64::_mm256_set1_epi32;
use std::arch::x86_64::_mm256_storeu_si256;

use super::Error;
use super::Outcome;
use super::Scoring;
use super::engine;
use super::wavefront;
use super::wavefront::Kernel;
use super::wavefront::LaneWidth;

/// Proof that the running processor supports AVX2.
///
/// The private field means [`Avx2::detect`] is the only way to obtain a value,
/// so any code holding an `Avx2` may call `#[target_feature(enable = "avx2")]`
/// functions.
#[derive(Clone, Copy, Debug)]
pub(super) struct Avx2 {
    /// Prevents construction outside [`Avx2::detect`].
    _private: (),
}

impl Avx2 {
    /// Returns a token when the running processor supports AVX2.
    ///
    /// The standard library caches the CPUID result, so repeated calls are
    /// cheap.
    pub(super) fn detect() -> Option<Self> {
        std::is_x86_feature_detected!("avx2").then_some(Self { _private: () })
    }
}

/// Returns an AVX2 token for tests, or `None` to skip AVX2 coverage.
///
/// Setting `OMICS_REQUIRE_AVX2=1` turns a missing AVX2 unit into a test
/// failure, so CI on AVX2 runners proves that the AVX2 kernels actually ran.
#[cfg(test)]
pub(super) fn avx2_for_tests() -> Option<Avx2> {
    let avx2 = Avx2::detect();
    if avx2.is_none() && std::env::var("OMICS_REQUIRE_AVX2").as_deref() == Ok("1") {
        panic!("OMICS_REQUIRE_AVX2=1 but this processor does not support AVX2");
    }
    avx2
}

/// A sixteen-lane signed sixteen-bit AVX2 wavefront kernel.
///
/// The kernel holds an [`Avx2`] token, so calling any of its methods proves
/// that the processor supports AVX2.
///
/// Each `*_avx2` helper enables `avx2` with `#[target_feature]`, so its
/// intrinsics are safe to call and only raw-pointer loads and stores need
/// `unsafe`. Safe trait methods cannot carry `#[target_feature]`, so each
/// [`Kernel`] method calls its helper in one `unsafe` block justified by the
/// token.
#[derive(Clone, Copy)]
pub(super) struct Avx2I16(pub(super) Avx2);

impl Avx2I16 {
    /// Loads sixteen signed sixteen-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2` and `pointer` addresses sixteen readable
    /// signed sixteen-bit values.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn load_avx2(pointer: *const i16) -> __m256i {
        // SAFETY: This function enables `avx2`, and its caller provides sixteen
        // readable signed sixteen-bit values at `pointer`. The unaligned load
        // accesses exactly those values.
        unsafe { _mm256_loadu_si256(pointer.cast()) }
    }

    /// Stores sixteen signed sixteen-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2` and `pointer` addresses sixteen writable
    /// signed sixteen-bit values.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn store_avx2(pointer: *mut i16, value: __m256i) {
        // SAFETY: This function enables `avx2`, and its caller provides sixteen
        // writable signed sixteen-bit values at `pointer`. The unaligned store
        // writes exactly those values.
        unsafe { _mm256_storeu_si256(pointer.cast(), value) };
    }

    /// Broadcasts a signed sixteen-bit score through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn splat_avx2(value: i16) -> __m256i {
        _mm256_set1_epi16(value)
    }

    /// Adds signed sixteen-bit lanes through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn add_avx2(left: __m256i, right: __m256i) -> __m256i {
        _mm256_add_epi16(left, right)
    }

    /// Compares signed sixteen-bit lanes through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn greater_than_avx2(left: __m256i, right: __m256i) -> __m256i {
        _mm256_cmpgt_epi16(left, right)
    }

    /// Compares signed sixteen-bit lane equality through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn equal_avx2(left: __m256i, right: __m256i) -> __m256i {
        _mm256_cmpeq_epi16(left, right)
    }

    /// Blends signed sixteen-bit lanes according to a vector mask through
    /// `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn select_avx2(mask: __m256i, if_true: __m256i, if_false: __m256i) -> __m256i {
        _mm256_blendv_epi8(if_false, if_true, mask)
    }

    /// Computes sixteen byte substitutions through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`. `reference` and `query_reversed` each
    /// address sixteen readable bytes.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn substitution_avx2(
        reference: *const u8,
        query_reversed: *const u8,
        match_score: i16,
        mismatch_score: i16,
    ) -> __m256i {
        // SAFETY: The caller provides sixteen readable bytes at `reference`.
        // The unaligned load accesses exactly those bytes.
        let reference_bytes = unsafe { _mm_loadu_si128(reference.cast::<__m128i>()) };
        // SAFETY: The caller provides sixteen readable bytes at
        // `query_reversed`. The unaligned load accesses exactly those bytes.
        let query_bytes = unsafe { _mm_loadu_si128(query_reversed.cast::<__m128i>()) };
        let equal_bytes = _mm_cmpeq_epi8(reference_bytes, query_bytes);
        let equal_masks = _mm256_cvtepi8_epi16(equal_bytes);

        _mm256_blendv_epi8(
            _mm256_set1_epi16(mismatch_score),
            _mm256_set1_epi16(match_score),
            equal_masks,
        )
    }
}

// SAFETY: Every method preserves the Kernel lane semantics for sixteen signed
// sixteen-bit lanes. A value of this kernel holds an `Avx2` token, so every
// method runs on a processor that supports `avx2`. `_mm256_add_epi16` wraps.
// The substitution method reads exactly sixteen bytes from each pointer.
unsafe impl Kernel for Avx2I16 {
    type Score = i16;
    type Vector = __m256i;

    const LANES: usize = 16;

    unsafe fn load(self, pointer: *const Self::Score) -> Self::Vector {
        // SAFETY: The Kernel contract provides sixteen readable signed
        // sixteen-bit lanes, and the `Avx2` token in `self` proves `avx2`.
        unsafe { Self::load_avx2(pointer) }
    }

    unsafe fn store(self, pointer: *mut Self::Score, value: Self::Vector) {
        // SAFETY: The Kernel contract provides sixteen writable signed
        // sixteen-bit lanes, and the `Avx2` token in `self` proves `avx2`.
        unsafe { Self::store_avx2(pointer, value) };
    }

    fn splat(self, value: Self::Score) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::splat_avx2(value) }
    }

    fn add(self, left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::add_avx2(left, right) }
    }

    fn greater_than(self, left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::greater_than_avx2(left, right) }
    }

    fn equal(self, left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::equal_avx2(left, right) }
    }

    fn select(
        self,
        mask: Self::Vector,
        if_true: Self::Vector,
        if_false: Self::Vector,
    ) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::select_avx2(mask, if_true, if_false) }
    }

    unsafe fn substitution(
        self,
        reference: *const u8,
        query_reversed: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector {
        // SAFETY: The Kernel contract supplies sixteen readable bytes at each
        // pointer. The `Avx2` token in `self` proves `avx2`.
        unsafe { Self::substitution_avx2(reference, query_reversed, match_score, mismatch_score) }
    }
}

/// An eight-lane signed thirty-two-bit AVX2 wavefront kernel.
///
/// The kernel holds an [`Avx2`] token, so calling any of its methods proves
/// that the processor supports AVX2. Its `unsafe` blocks follow the same rules
/// as [`Avx2I16`].
#[derive(Clone, Copy)]
pub(super) struct Avx2I32(pub(super) Avx2);

impl Avx2I32 {
    /// Loads eight signed thirty-two-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2` and `pointer` addresses eight readable
    /// signed thirty-two-bit values.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn load_avx2(pointer: *const i32) -> __m256i {
        // SAFETY: This function enables `avx2`, and its caller provides eight
        // readable signed thirty-two-bit values at `pointer`. The unaligned
        // load accesses exactly those values.
        unsafe { _mm256_loadu_si256(pointer.cast()) }
    }

    /// Stores eight signed thirty-two-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2` and `pointer` addresses eight writable
    /// signed thirty-two-bit values.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn store_avx2(pointer: *mut i32, value: __m256i) {
        // SAFETY: This function enables `avx2`, and its caller provides eight
        // writable signed thirty-two-bit values at `pointer`. The unaligned
        // store writes exactly those values.
        unsafe { _mm256_storeu_si256(pointer.cast(), value) };
    }

    /// Broadcasts a signed thirty-two-bit score through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn splat_avx2(value: i32) -> __m256i {
        _mm256_set1_epi32(value)
    }

    /// Adds signed thirty-two-bit lanes through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn add_avx2(left: __m256i, right: __m256i) -> __m256i {
        _mm256_add_epi32(left, right)
    }

    /// Compares signed thirty-two-bit lanes through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn greater_than_avx2(left: __m256i, right: __m256i) -> __m256i {
        _mm256_cmpgt_epi32(left, right)
    }

    /// Compares signed thirty-two-bit lane equality through `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn equal_avx2(left: __m256i, right: __m256i) -> __m256i {
        _mm256_cmpeq_epi32(left, right)
    }

    /// Blends signed thirty-two-bit lanes according to a vector mask through
    /// `avx2`.
    #[inline]
    #[target_feature(enable = "avx2")]
    fn select_avx2(mask: __m256i, if_true: __m256i, if_false: __m256i) -> __m256i {
        _mm256_blendv_epi8(if_false, if_true, mask)
    }

    /// Computes eight byte substitutions through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`. `reference` and `query_reversed` each
    /// address eight readable bytes.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn substitution_avx2(
        reference: *const u8,
        query_reversed: *const u8,
        match_score: i32,
        mismatch_score: i32,
    ) -> __m256i {
        // SAFETY: The caller provides eight readable reference bytes at
        // `reference`. The load accesses exactly those bytes.
        let reference_bytes = unsafe { _mm_loadl_epi64(reference.cast::<__m128i>()) };
        // SAFETY: The caller provides eight readable bytes at
        // `query_reversed`. The load accesses exactly those bytes.
        let query_bytes = unsafe { _mm_loadl_epi64(query_reversed.cast::<__m128i>()) };
        let reference_values = _mm256_cvtepi8_epi32(reference_bytes);
        let query_values = _mm256_cvtepi8_epi32(query_bytes);
        let equal_masks = _mm256_cmpeq_epi32(reference_values, query_values);

        _mm256_blendv_epi8(
            _mm256_set1_epi32(mismatch_score),
            _mm256_set1_epi32(match_score),
            equal_masks,
        )
    }
}

// SAFETY: Every method preserves the Kernel lane semantics for eight signed
// thirty-two-bit lanes. A value of this kernel holds an `Avx2` token, so every
// method runs on a processor that supports `avx2`. `_mm256_add_epi32` wraps.
// The substitution method reads exactly eight bytes from each pointer.
unsafe impl Kernel for Avx2I32 {
    type Score = i32;
    type Vector = __m256i;

    const LANES: usize = 8;

    unsafe fn load(self, pointer: *const Self::Score) -> Self::Vector {
        // SAFETY: The Kernel contract provides eight readable signed
        // thirty-two-bit lanes, and the `Avx2` token in `self` proves `avx2`.
        unsafe { Self::load_avx2(pointer) }
    }

    unsafe fn store(self, pointer: *mut Self::Score, value: Self::Vector) {
        // SAFETY: The Kernel contract provides eight writable signed
        // thirty-two-bit lanes, and the `Avx2` token in `self` proves `avx2`.
        unsafe { Self::store_avx2(pointer, value) };
    }

    fn splat(self, value: Self::Score) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::splat_avx2(value) }
    }

    fn add(self, left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::add_avx2(left, right) }
    }

    fn greater_than(self, left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::greater_than_avx2(left, right) }
    }

    fn equal(self, left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::equal_avx2(left, right) }
    }

    fn select(
        self,
        mask: Self::Vector,
        if_true: Self::Vector,
        if_false: Self::Vector,
    ) -> Self::Vector {
        // SAFETY: `self` holds an `Avx2` token, which proves `avx2` support.
        unsafe { Self::select_avx2(mask, if_true, if_false) }
    }

    unsafe fn substitution(
        self,
        reference: *const u8,
        query_reversed: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector {
        // SAFETY: The Kernel contract supplies eight readable bytes at each
        // pointer. The `Avx2` token in `self` proves `avx2`.
        unsafe { Self::substitution_avx2(reference, query_reversed, match_score, mismatch_score) }
    }
}

/// Returns the global alignment through the narrowest eligible AVX2 kernel.
pub(super) fn global(
    avx2: Avx2,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Outcome, Error> {
    // SAFETY: The `Avx2` token proves that the processor supports `avx2`.
    unsafe { global_avx2(avx2, reference, query, scoring) }
}

/// Returns the local alignment through the narrowest eligible AVX2 kernel.
pub(super) fn local(
    avx2: Avx2,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    // SAFETY: The `Avx2` token proves that the processor supports `avx2`.
    unsafe { local_avx2(avx2, reference, query, scoring) }
}

/// Returns the global alignment with `avx2` enabled for the inlined engine.
///
/// The wavefront engine is `#[inline(always)]` on this target, so enabling
/// `avx2` here lets the compiler inline the kernel methods into one loop.
///
/// Callers without `avx2` enabled must prove support, for example with an
/// [`Avx2`] token, before calling it in an `unsafe` block.
#[target_feature(enable = "avx2")]
fn global_avx2(
    avx2: Avx2,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Outcome, Error> {
    match wavefront::select_lane_width("global", reference.len(), query.len(), scoring) {
        LaneWidth::I16 => wavefront::global(Avx2I16(avx2), reference, query, scoring),
        LaneWidth::I32 => wavefront::global(Avx2I32(avx2), reference, query, scoring),
        LaneWidth::Scalar => engine::global(reference, query, scoring),
    }
}

/// Returns the local alignment with `avx2` enabled for the inlined engine.
///
/// See [`global_avx2`] for how the engine is inlined.
///
/// Callers without `avx2` enabled must prove support, for example with an
/// [`Avx2`] token, before calling it in an `unsafe` block.
#[target_feature(enable = "avx2")]
fn local_avx2(
    avx2: Avx2,
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    match wavefront::select_lane_width("local", reference.len(), query.len(), scoring) {
        LaneWidth::I16 => wavefront::local(Avx2I16(avx2), reference, query, scoring),
        LaneWidth::I32 => wavefront::local(Avx2I32(avx2), reference, query, scoring),
        LaneWidth::Scalar => engine::local(reference, query, scoring),
    }
}
