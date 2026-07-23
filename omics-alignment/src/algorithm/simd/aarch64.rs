use std::arch::aarch64::int16x8_t;
use std::arch::aarch64::int32x4_t;
use std::arch::aarch64::vaddq_s16;
use std::arch::aarch64::vaddq_s32;
use std::arch::aarch64::vbslq_s16;
use std::arch::aarch64::vbslq_s32;
use std::arch::aarch64::vceq_u8;
use std::arch::aarch64::vceqq_s16;
use std::arch::aarch64::vceqq_s32;
use std::arch::aarch64::vcgtq_s16;
use std::arch::aarch64::vcgtq_s32;
use std::arch::aarch64::vdupq_n_s16;
use std::arch::aarch64::vdupq_n_s32;
use std::arch::aarch64::vld1_u8;
use std::arch::aarch64::vld1q_s16;
use std::arch::aarch64::vld1q_s32;
use std::arch::aarch64::vmovl_s8;
use std::arch::aarch64::vreinterpret_s8_u8;
use std::arch::aarch64::vreinterpretq_s16_u16;
use std::arch::aarch64::vreinterpretq_s32_u32;
use std::arch::aarch64::vreinterpretq_u16_s16;
use std::arch::aarch64::vreinterpretq_u32_s32;
use std::arch::aarch64::vrev64_u8;
use std::arch::aarch64::vst1q_s16;
use std::arch::aarch64::vst1q_s32;

use super::Error;
use super::Outcome;
use super::Scoring;
use super::engine;
use super::wavefront;
use super::wavefront::Kernel;
use super::wavefront::LaneWidth;

/// An eight-lane signed sixteen-bit NEON wavefront kernel.
#[derive(Clone, Copy)]
pub(super) struct NeonI16;

// SAFETY: Every method preserves the Kernel lane semantics for eight signed
// sixteen-bit lanes. The substitution method consumes eight ascending
// reference bytes and eight descending query bytes.
// NOTE: Rust 1.81 requires these intrinsic calls inside unsafe blocks while
// newer toolchains mark them safe, so `#[expect]` cannot work consistently.
#[allow(unused_unsafe)]
unsafe impl Kernel for NeonI16 {
    type Score = i16;
    type Vector = int16x8_t;

    const LANES: usize = 8;

    #[target_feature(enable = "neon")]
    unsafe fn load(pointer: *const Self::Score) -> Self::Vector {
        // SAFETY: This method enables `neon`, and the Kernel contract provides
        // eight readable signed sixteen-bit lanes at `pointer`.
        unsafe { vld1q_s16(pointer) }
    }

    #[target_feature(enable = "neon")]
    unsafe fn store(pointer: *mut Self::Score, value: Self::Vector) {
        // SAFETY: This method enables `neon`, and the Kernel contract provides
        // eight writable signed sixteen-bit lanes at `pointer`.
        unsafe { vst1q_s16(pointer, value) };
    }

    fn splat(value: Self::Score) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vdupq_n_s16(value) }
    }

    fn add(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vaddq_s16(left, right) }
    }

    fn greater_than(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vreinterpretq_s16_u16(vcgtq_s16(left, right)) }
    }

    fn equal(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vreinterpretq_s16_u16(vceqq_s16(left, right)) }
    }

    fn select(mask: Self::Vector, if_true: Self::Vector, if_false: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vbslq_s16(vreinterpretq_u16_s16(mask), if_true, if_false) }
    }

    #[target_feature(enable = "neon")]
    unsafe fn substitution(
        reference: *const u8,
        query_end: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector {
        // SAFETY: This method enables `neon`, and the Kernel contract provides
        // eight readable reference bytes at `reference`.
        let reference_bytes = unsafe { vld1_u8(reference) };
        // SAFETY: The Kernel contract provides eight query bytes immediately
        // before `query_end`.
        let query_start = unsafe { query_end.sub(Self::LANES) };
        // SAFETY: This method enables `neon`, and `query_start` begins the
        // eight readable query bytes established by the Kernel contract.
        let query_bytes = unsafe { vld1_u8(query_start) };
        // SAFETY: This method enables `neon`, and the loaded vectors contain
        // the eight byte lanes required by the reversal, comparison, widening,
        // and score-selection intrinsics.
        unsafe {
            let query_bytes = vrev64_u8(query_bytes);
            let equal_bytes = vceq_u8(reference_bytes, query_bytes);
            let equal_masks = vreinterpretq_u16_s16(vmovl_s8(vreinterpret_s8_u8(equal_bytes)));

            vbslq_s16(
                equal_masks,
                vdupq_n_s16(match_score),
                vdupq_n_s16(mismatch_score),
            )
        }
    }
}

/// A four-lane signed thirty-two-bit NEON wavefront kernel.
#[derive(Clone, Copy)]
pub(super) struct NeonI32;

// SAFETY: Every method preserves the Kernel lane semantics for four signed
// thirty-two-bit lanes. The substitution method consumes four ascending
// reference bytes and four descending query bytes.
// NOTE: Rust 1.81 requires these intrinsic calls inside unsafe blocks while
// newer toolchains mark them safe, so `#[expect]` cannot work consistently.
#[allow(unused_unsafe)]
unsafe impl Kernel for NeonI32 {
    type Score = i32;
    type Vector = int32x4_t;

    const LANES: usize = 4;

    #[target_feature(enable = "neon")]
    unsafe fn load(pointer: *const Self::Score) -> Self::Vector {
        // SAFETY: This method enables `neon`, and the Kernel contract provides
        // four readable signed thirty-two-bit lanes at `pointer`.
        unsafe { vld1q_s32(pointer) }
    }

    #[target_feature(enable = "neon")]
    unsafe fn store(pointer: *mut Self::Score, value: Self::Vector) {
        // SAFETY: This method enables `neon`, and the Kernel contract provides
        // four writable signed thirty-two-bit lanes at `pointer`.
        unsafe { vst1q_s32(pointer, value) };
    }

    fn splat(value: Self::Score) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vdupq_n_s32(value) }
    }

    fn add(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vaddq_s32(left, right) }
    }

    fn greater_than(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vreinterpretq_s32_u32(vcgtq_s32(left, right)) }
    }

    fn equal(left: Self::Vector, right: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vreinterpretq_s32_u32(vceqq_s32(left, right)) }
    }

    fn select(mask: Self::Vector, if_true: Self::Vector, if_false: Self::Vector) -> Self::Vector {
        // SAFETY: This module only compiles for AArch64 macOS, where `neon`
        // is mandatory.
        unsafe { vbslq_s32(vreinterpretq_u32_s32(mask), if_true, if_false) }
    }

    #[target_feature(enable = "neon")]
    unsafe fn substitution(
        reference: *const u8,
        query_end: *const u8,
        match_score: Self::Score,
        mismatch_score: Self::Score,
    ) -> Self::Vector {
        // SAFETY: The Kernel contract provides four readable reference bytes
        // at `reference`.
        let reference_bytes = unsafe { std::slice::from_raw_parts(reference, Self::LANES) };
        // SAFETY: The Kernel contract provides four query bytes immediately
        // before `query_end`.
        let query_start = unsafe { query_end.sub(Self::LANES) };
        // SAFETY: The Kernel contract provides four readable query bytes at
        // `query_start`.
        let query_bytes = unsafe { std::slice::from_raw_parts(query_start, Self::LANES) };
        let reference_values = [
            i32::from(reference_bytes[0]),
            i32::from(reference_bytes[1]),
            i32::from(reference_bytes[2]),
            i32::from(reference_bytes[3]),
        ];
        let query_values = [
            i32::from(query_bytes[3]),
            i32::from(query_bytes[2]),
            i32::from(query_bytes[1]),
            i32::from(query_bytes[0]),
        ];
        // SAFETY: This method enables `neon`, and the four-element local
        // staging array provides four readable signed thirty-two-bit lanes.
        let reference_values = unsafe { vld1q_s32(reference_values.as_ptr()) };
        // SAFETY: This method enables `neon`, and the four-element local
        // staging array provides four readable signed thirty-two-bit lanes.
        let query_values = unsafe { vld1q_s32(query_values.as_ptr()) };
        // SAFETY: This method enables `neon`, and the loaded vectors contain
        // the four lanes required by the comparison and score-selection
        // intrinsics.
        unsafe {
            let equal_masks = vceqq_s32(reference_values, query_values);

            vbslq_s32(
                equal_masks,
                vdupq_n_s32(match_score),
                vdupq_n_s32(mismatch_score),
            )
        }
    }
}

/// Returns the global alignment through the narrowest eligible NEON kernel.
pub(super) fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    match wavefront::lane_width(reference.len(), query.len(), scoring) {
        LaneWidth::I16 => wavefront::global::<NeonI16>(reference, query, scoring),
        LaneWidth::I32 => wavefront::global::<NeonI32>(reference, query, scoring),
        LaneWidth::Scalar => engine::global(reference, query, scoring),
    }
}

/// Returns the local alignment through the narrowest eligible NEON kernel.
pub(super) fn local(
    reference: &[u8],
    query: &[u8],
    scoring: Scoring,
) -> Result<Option<Outcome>, Error> {
    match wavefront::lane_width(reference.len(), query.len(), scoring) {
        LaneWidth::I16 => wavefront::local::<NeonI16>(reference, query, scoring),
        LaneWidth::I32 => wavefront::local::<NeonI32>(reference, query, scoring),
        LaneWidth::Scalar => engine::local(reference, query, scoring),
    }
}
