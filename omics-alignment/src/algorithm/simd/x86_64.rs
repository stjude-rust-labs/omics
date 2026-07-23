use std::arch::x86_64::__m128i;
use std::arch::x86_64::__m256i;
use std::arch::x86_64::_mm_cmpeq_epi8;
use std::arch::x86_64::_mm_loadu_si128;
use std::arch::x86_64::_mm_setr_epi8;
use std::arch::x86_64::_mm_shuffle_epi8;
use std::arch::x86_64::_mm256_add_epi16;
use std::arch::x86_64::_mm256_add_epi32;
use std::arch::x86_64::_mm256_blendv_epi8;
use std::arch::x86_64::_mm256_cmpeq_epi16;
use std::arch::x86_64::_mm256_cmpeq_epi32;
use std::arch::x86_64::_mm256_cmpgt_epi16;
use std::arch::x86_64::_mm256_cmpgt_epi32;
use std::arch::x86_64::_mm256_cvtepi8_epi16;
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

/// A sixteen-lane signed sixteen-bit AVX2 wavefront kernel.
#[derive(Clone, Copy)]
pub(super) struct Avx2I16;

// NOTE: Rust 1.81 requires these intrinsic calls inside unsafe blocks while
// newer toolchains mark them safe, so `#[expect]` cannot work consistently.
#[allow(unused_unsafe)]
impl Avx2I16 {
    /// Loads sixteen signed sixteen-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2` and `pointer` addresses sixteen readable
    /// signed sixteen-bit values.
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
    #[target_feature(enable = "avx2")]
    unsafe fn store_avx2(pointer: *mut i16, value: __m256i) {
        // SAFETY: This function enables `avx2`, and its caller provides sixteen
        // writable signed sixteen-bit values at `pointer`. The unaligned store
        // writes exactly those values.
        unsafe { _mm256_storeu_si256(pointer.cast(), value) };
    }

    /// Broadcasts a signed sixteen-bit score through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn splat_avx2(value: i16) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_set1_epi16(value) }
    }

    /// Adds signed sixteen-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn add_avx2(left: __m256i, right: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_add_epi16(left, right) }
    }

    /// Compares signed sixteen-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn greater_than_avx2(left: __m256i, right: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_cmpgt_epi16(left, right) }
    }

    /// Compares signed sixteen-bit lane equality through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn equal_avx2(left: __m256i, right: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_cmpeq_epi16(left, right) }
    }

    /// Blends signed sixteen-bit lanes according to a vector mask through
    /// `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn select_avx2(mask: __m256i, if_true: __m256i, if_false: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_blendv_epi8(if_false, if_true, mask) }
    }

    /// Computes sixteen byte substitutions through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`. `reference` addresses sixteen readable
    /// bytes and `query_end` follows sixteen readable bytes.
    #[target_feature(enable = "avx2")]
    unsafe fn substitution_avx2(
        reference: *const u8,
        query_end: *const u8,
        match_score: i16,
        mismatch_score: i16,
    ) -> __m256i {
        // SAFETY: This function enables `avx2`, and its caller provides sixteen
        // readable bytes at `reference`. The unaligned load accesses exactly
        // those bytes.
        let reference_bytes = unsafe { _mm_loadu_si128(reference.cast::<__m128i>()) };
        // SAFETY: The caller provides sixteen readable query bytes immediately
        // before `query_end`.
        let query_start = unsafe { query_end.sub(Self::LANES) };
        // SAFETY: This function enables `avx2`, and `query_start` begins the
        // sixteen readable query bytes established by the caller.
        let query_bytes = unsafe { _mm_loadu_si128(query_start.cast::<__m128i>()) };
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports the byte shuffle instructions.
        let query_bytes = unsafe {
            _mm_shuffle_epi8(
                query_bytes,
                _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0),
            )
        };
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports the byte comparison instruction.
        let equal_bytes = unsafe { _mm_cmpeq_epi8(reference_bytes, query_bytes) };
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports the byte-to-word widening
        // instruction.
        let equal_masks = unsafe { _mm256_cvtepi8_epi16(equal_bytes) };

        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports the blend and broadcast
        // instructions.
        unsafe {
            _mm256_blendv_epi8(
                _mm256_set1_epi16(mismatch_score),
                _mm256_set1_epi16(match_score),
                equal_masks,
            )
        }
    }
}

// SAFETY: Every method preserves the Kernel lane semantics for sixteen signed
// sixteen-bit lanes. The production wrappers and native differential test
// verify `avx2` before this kernel is used. The substitution method consumes
// sixteen ascending reference bytes and sixteen descending query bytes.
unsafe impl Kernel for Avx2I16 {
    type Score = i16;
    type Vector = __m256i;

    const LANES: usize = 16;

    unsafe fn load(pointer: *const Self::Score) -> Self::Vector {
        // SAFETY: The Kernel contract provides sixteen readable signed
        // sixteen-bit lanes, and every AVX2 caller verifies `avx2`.
        unsafe { Self::load_avx2(pointer) }
    }

    unsafe fn store(pointer: *mut Self::Score, value: Self::Vector) {
        // SAFETY: The Kernel contract provides sixteen writable signed
        // sixteen-bit lanes, and every AVX2 caller verifies `avx2`.
        unsafe { Self::store_avx2(pointer, value) };
    }

    fn splat(value: Self::Score) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::splat_avx2(value) }
    }

    fn add(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::add_avx2(left, right) }
    }

    fn greater_than(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::greater_than_avx2(left, right) }
    }

    fn equal(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::equal_avx2(left, right) }
    }

    fn select(mask: Self::Vector, if_true: Self::Vector, if_false: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::select_avx2(mask, if_true, if_false) }
    }

    unsafe fn substitution(
        reference: *const u8,
        query_end: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector {
        // SAFETY: The Kernel contract supplies sixteen readable reference
        // bytes and sixteen readable query bytes before `query_end`. Every
        // AVX2 caller verifies `avx2`.
        unsafe { Self::substitution_avx2(reference, query_end, match_score, mismatch_score) }
    }
}

/// An eight-lane signed thirty-two-bit AVX2 wavefront kernel.
#[derive(Clone, Copy)]
pub(super) struct Avx2I32;

// NOTE: Rust 1.81 requires these intrinsic calls inside unsafe blocks while
// newer toolchains mark them safe, so `#[expect]` cannot work consistently.
#[allow(unused_unsafe)]
impl Avx2I32 {
    /// Loads eight signed thirty-two-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2` and `pointer` addresses eight readable
    /// signed thirty-two-bit values.
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
    #[target_feature(enable = "avx2")]
    unsafe fn store_avx2(pointer: *mut i32, value: __m256i) {
        // SAFETY: This function enables `avx2`, and its caller provides eight
        // writable signed thirty-two-bit values at `pointer`. The unaligned
        // store writes exactly those values.
        unsafe { _mm256_storeu_si256(pointer.cast(), value) };
    }

    /// Broadcasts a signed thirty-two-bit score through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn splat_avx2(value: i32) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_set1_epi32(value) }
    }

    /// Adds signed thirty-two-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn add_avx2(left: __m256i, right: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_add_epi32(left, right) }
    }

    /// Compares signed thirty-two-bit lanes through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn greater_than_avx2(left: __m256i, right: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_cmpgt_epi32(left, right) }
    }

    /// Compares signed thirty-two-bit lane equality through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn equal_avx2(left: __m256i, right: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_cmpeq_epi32(left, right) }
    }

    /// Blends signed thirty-two-bit lanes according to a vector mask through
    /// `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`.
    #[target_feature(enable = "avx2")]
    unsafe fn select_avx2(mask: __m256i, if_true: __m256i, if_false: __m256i) -> __m256i {
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports it.
        unsafe { _mm256_blendv_epi8(if_false, if_true, mask) }
    }

    /// Computes eight byte substitutions through `avx2`.
    ///
    /// # Safety
    ///
    /// The processor supports `avx2`. `reference` addresses eight readable
    /// bytes and `query_end` follows eight readable bytes.
    #[target_feature(enable = "avx2")]
    unsafe fn substitution_avx2(
        reference: *const u8,
        query_end: *const u8,
        match_score: i32,
        mismatch_score: i32,
    ) -> __m256i {
        // SAFETY: The caller provides eight readable reference bytes at
        // `reference`.
        let reference_bytes = unsafe { std::slice::from_raw_parts(reference, Self::LANES) };
        // SAFETY: The caller provides eight readable query bytes immediately
        // before `query_end`.
        let query_start = unsafe { query_end.sub(Self::LANES) };
        // SAFETY: `query_start` begins the eight readable query bytes
        // established by the caller.
        let query_bytes = unsafe { std::slice::from_raw_parts(query_start, Self::LANES) };
        let reference_values = [
            i32::from(reference_bytes[0]),
            i32::from(reference_bytes[1]),
            i32::from(reference_bytes[2]),
            i32::from(reference_bytes[3]),
            i32::from(reference_bytes[4]),
            i32::from(reference_bytes[5]),
            i32::from(reference_bytes[6]),
            i32::from(reference_bytes[7]),
        ];
        let query_values = [
            i32::from(query_bytes[7]),
            i32::from(query_bytes[6]),
            i32::from(query_bytes[5]),
            i32::from(query_bytes[4]),
            i32::from(query_bytes[3]),
            i32::from(query_bytes[2]),
            i32::from(query_bytes[1]),
            i32::from(query_bytes[0]),
        ];
        // SAFETY: This function enables `avx2`, and the local staging array
        // contains eight readable signed thirty-two-bit values.
        let reference_values = unsafe { _mm256_loadu_si256(reference_values.as_ptr().cast()) };
        // SAFETY: This function enables `avx2`, and the local staging array
        // contains eight readable signed thirty-two-bit values.
        let query_values = unsafe { _mm256_loadu_si256(query_values.as_ptr().cast()) };
        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports the lane comparison instruction.
        let equal_masks = unsafe { _mm256_cmpeq_epi32(reference_values, query_values) };

        // SAFETY: This target-feature function uses `avx2`, and its caller
        // proves that the processor supports the blend and broadcast
        // instructions.
        unsafe {
            _mm256_blendv_epi8(
                _mm256_set1_epi32(mismatch_score),
                _mm256_set1_epi32(match_score),
                equal_masks,
            )
        }
    }
}

// SAFETY: Every method preserves the Kernel lane semantics for eight signed
// thirty-two-bit lanes. The production wrappers and native differential test
// verify `avx2` before this kernel is used. The substitution method consumes
// eight ascending reference bytes and eight descending query bytes.
unsafe impl Kernel for Avx2I32 {
    type Score = i32;
    type Vector = __m256i;

    const LANES: usize = 8;

    unsafe fn load(pointer: *const Self::Score) -> Self::Vector {
        // SAFETY: The Kernel contract provides eight readable signed
        // thirty-two-bit lanes, and every AVX2 caller verifies `avx2`.
        unsafe { Self::load_avx2(pointer) }
    }

    unsafe fn store(pointer: *mut Self::Score, value: Self::Vector) {
        // SAFETY: The Kernel contract provides eight writable signed
        // thirty-two-bit lanes, and every AVX2 caller verifies `avx2`.
        unsafe { Self::store_avx2(pointer, value) };
    }

    fn splat(value: Self::Score) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::splat_avx2(value) }
    }

    fn add(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::add_avx2(left, right) }
    }

    fn greater_than(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::greater_than_avx2(left, right) }
    }

    fn equal(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::equal_avx2(left, right) }
    }

    fn select(mask: Self::Vector, if_true: Self::Vector, if_false: Self::Vector) -> Self::Vector {
        // SAFETY: Every AVX2 caller verifies `avx2` before invoking this
        // kernel.
        unsafe { Self::select_avx2(mask, if_true, if_false) }
    }

    unsafe fn substitution(
        reference: *const u8,
        query_end: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector {
        // SAFETY: The Kernel contract supplies eight readable reference bytes
        // and eight readable query bytes before `query_end`. Every AVX2 caller
        // verifies `avx2`.
        unsafe { Self::substitution_avx2(reference, query_end, match_score, mismatch_score) }
    }
}

/// Returns the global alignment through AVX2 when the processor supports it.
pub(super) fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: The runtime feature check proves that this processor supports
        // `avx2` for the target-feature implementation.
        unsafe { global_avx2(reference, query, scoring) }
    } else {
        engine::global(reference, query, scoring)
    }
}

/// Returns the local alignment through AVX2 when the processor supports it.
pub(super) fn local(
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: The runtime feature check proves that this processor supports
        // `avx2` for the target-feature implementation.
        unsafe { local_avx2(reference, query, scoring) }
    } else {
        engine::local(reference, query, scoring)
    }
}

/// Returns the global alignment after the AVX2 feature check.
///
/// # Safety
///
/// The processor supports `avx2`.
#[target_feature(enable = "avx2")]
unsafe fn global_avx2(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    match wavefront::lane_width(reference.len(), query.len(), scoring) {
        LaneWidth::I16 => wavefront::global::<Avx2I16>(reference, query, scoring),
        LaneWidth::I32 => wavefront::global::<Avx2I32>(reference, query, scoring),
        LaneWidth::Scalar => engine::global(reference, query, scoring),
    }
}

/// Returns the local alignment after the AVX2 feature check.
///
/// # Safety
///
/// The processor supports `avx2`.
#[target_feature(enable = "avx2")]
unsafe fn local_avx2(
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    match wavefront::lane_width(reference.len(), query.len(), scoring) {
        LaneWidth::I16 => wavefront::local::<Avx2I16>(reference, query, scoring),
        LaneWidth::I32 => wavefront::local::<Avx2I32>(reference, query, scoring),
        LaneWidth::Scalar => engine::local(reference, query, scoring),
    }
}
