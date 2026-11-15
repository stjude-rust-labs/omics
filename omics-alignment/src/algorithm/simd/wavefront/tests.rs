//! Differential and fixture tests for the wavefront engine.
//!
//! Portable array kernels with `i16` and `i32` lanes run the generic engine
//! without SIMD instructions (and under Miri), and are compared with the
//! scalar engine over exhaustive inputs, deterministic patterns, and
//! hand-worked fixtures. Where a native backend is available, the same
//! comparisons run on it.

use std::marker::PhantomData;
use std::ops::Range;

use super::traceback::pack_predecessors;
use super::traceback::unpack_predecessor;
use super::*;
use crate::algorithm::Score;
use crate::algorithm::engine;

/// A lane type supported by [`ArrayKernel`].
trait ArrayLane: LaneScore {
    /// A nonzero mask lane.
    const TRUE: Self;

    /// Adds two lanes with wrapping arithmetic.
    fn wrapping_add(self, other: Self) -> Self;
}

impl ArrayLane for i16 {
    const TRUE: Self = -1;

    fn wrapping_add(self, other: Self) -> Self {
        i16::wrapping_add(self, other)
    }
}

impl ArrayLane for i32 {
    const TRUE: Self = -1;

    fn wrapping_add(self, other: Self) -> Self {
        i32::wrapping_add(self, other)
    }
}

/// A portable kernel whose vectors are plain arrays of `N` lanes.
///
/// It runs the generic engine without SIMD intrinsics, so the engine's
/// vector paths are testable on every target and under Miri.
#[derive(Clone, Copy)]
struct ArrayKernel<S, const N: usize>(PhantomData<S>);

/// An eight-lane sixteen-bit test kernel (the NEON `i16` lane count).
type TestI16 = ArrayKernel<i16, 8>;
/// A four-lane thirty-two-bit test kernel (the NEON `i32` lane count).
type TestI32 = ArrayKernel<i32, 4>;

/// The sixteen-bit test kernel value.
const TEST_I16: TestI16 = ArrayKernel(PhantomData);
/// The thirty-two-bit test kernel value.
const TEST_I32: TestI32 = ArrayKernel(PhantomData);

// SAFETY: `N` is nonzero for every instantiation above, and every
// operation implements the documented lane-wise behavior for exactly `N`
// elements.
unsafe impl<S: ArrayLane, const N: usize> Kernel for ArrayKernel<S, N> {
    type Score = S;
    type Vector = [S; N];

    const LANES: usize = N;

    unsafe fn load(self, pointer: *const S) -> [S; N] {
        // SAFETY: The Kernel caller provides `N` readable scores.
        std::array::from_fn(|lane| unsafe { pointer.add(lane).read() })
    }

    unsafe fn store(self, pointer: *mut S, value: [S; N]) {
        for (lane, value) in value.into_iter().enumerate() {
            // SAFETY: The Kernel caller provides `N` writable scores.
            unsafe { pointer.add(lane).write(value) };
        }
    }

    fn splat(self, value: S) -> [S; N] {
        [value; N]
    }

    fn add(self, left: [S; N], right: [S; N]) -> [S; N] {
        std::array::from_fn(|lane| left[lane].wrapping_add(right[lane]))
    }

    fn greater_than(self, left: [S; N], right: [S; N]) -> [S; N] {
        std::array::from_fn(|lane| {
            if left[lane] > right[lane] {
                S::TRUE
            } else {
                S::ZERO
            }
        })
    }

    fn equal(self, left: [S; N], right: [S; N]) -> [S; N] {
        std::array::from_fn(|lane| {
            if left[lane] == right[lane] {
                S::TRUE
            } else {
                S::ZERO
            }
        })
    }

    fn select(self, mask: [S; N], if_true: [S; N], if_false: [S; N]) -> [S; N] {
        std::array::from_fn(|lane| {
            if mask[lane] != S::ZERO {
                if_true[lane]
            } else {
                if_false[lane]
            }
        })
    }

    unsafe fn substitution(
        self,
        reference: *const u8,
        query_reversed: *const u8,
        match_score: S,
        mismatch_score: S,
    ) -> [S; N] {
        std::array::from_fn(|lane| {
            // SAFETY: The Kernel caller supplies `N` bytes at each pointer.
            let (reference_byte, query_byte) =
                unsafe { (reference.add(lane).read(), query_reversed.add(lane).read()) };
            if reference_byte == query_byte {
                match_score
            } else {
                mismatch_score
            }
        })
    }
}

/// A hand-worked example: `ACG` (rows) against `AG` (columns) with match 2,
/// mismatch -1, gap open -2, and gap extend -1.
///
/// Each cell holds the best score over the three states, which is the global
/// score of `reference[..row]` against `query[..column]`. For example:
///
/// - `(1, 1)`: `A` matches `A`, so `0 + 2 = 2`.
/// - `(2, 1)`: deleting `C` after `(1, 1)` opens a gap, so `2 - 2 = 0`.
/// - `(2, 2)`: `C` mismatches `G` after `(1, 1)`, so `2 - 1 = 1`.
/// - `(3, 2)`: `G` matches `G` after `(2, 1)`, so `0 + 2 = 2`, the global score
///   with CIGAR `1=1D1=`.
///
/// Every cell is checked against the scalar engine and both test kernels.
/// The origin `(0, 0)` scores zero but is not a valid global alignment.
#[test]
fn test_kernels_match_a_hand_worked_example() -> Result<(), Error> {
    let reference = b"ACG";
    let query = b"AG";
    let scoring = Scoring::try_new(2, -1, -2, -1)?;
    #[rustfmt::skip]
    let expected: [[Score; 3]; 4] = [
        //  ""   A   G
        [   0, -2, -3], // ""
        [  -2,  2,  0], // A
        [  -3,  0,  1], // AC
        [  -4, -1,  2], // ACG
    ];

    for (row, scores) in expected.iter().enumerate() {
        for (column, &score) in scores.iter().enumerate() {
            let reference = &reference[..row];
            let query = &query[..column];
            let results = [
                engine::global(reference, query, scoring),
                global(TEST_I16, reference, query, scoring),
                global(TEST_I32, reference, query, scoring),
            ];

            for result in results {
                if (row, column) == (0, 0) {
                    assert!(matches!(result, Err(Error::EmptyGlobal)));
                } else {
                    assert_eq!(result?.score(), score, "cell ({row}, {column})");
                }
            }
        }
    }

    for outcome in [
        engine::global(reference, query, scoring)?,
        global(TEST_I16, reference, query, scoring)?,
        global(TEST_I32, reference, query, scoring)?,
    ] {
        assert_eq!(outcome.score(), 2);
        assert_eq!(outcome.cigar().to_string(), "1=1D1=");
        assert_eq!(outcome.reference_range(), &(0..3));
        assert_eq!(outcome.query_range(), &(0..2));
    }

    // Locally, `(1, 1)` and `(3, 2)` both score 2; the first in row-major
    // order wins.
    for outcome in [
        engine::local(reference, query, scoring)?,
        local(TEST_I16, reference, query, scoring)?,
        local(TEST_I32, reference, query, scoring)?,
    ] {
        let Some(outcome) = outcome else {
            panic!("the example has a positive local alignment");
        };
        assert_eq!(outcome.score(), 2);
        assert_eq!(outcome.cigar().to_string(), "1=");
        assert_eq!(outcome.reference_range(), &(0..1));
        assert_eq!(outcome.query_range(), &(0..1));
    }

    Ok(())
}

#[test]
fn vector_values_broadcast_recurrence_constants() -> Result<(), Error> {
    let lane_scoring = LaneScoring::<i16>::try_new(Scoring::try_new(2, -3, -2, -1)?)?;
    let broadcast = BroadcastScoring::new(TEST_I16, &lane_scoring);

    assert_eq!(broadcast.unreachable, [i16::MIN; TestI16::LANES]);
    assert_eq!(broadcast.zero, [0; TestI16::LANES]);
    assert_eq!(broadcast.gap_open_score, [-2; TestI16::LANES]);
    assert_eq!(broadcast.gap_extend_score, [-1; TestI16::LANES]);
    for (destination, shift) in [
        (State::Aligned, 0),
        (State::Insertion, 2),
        (State::Deletion, 4),
    ] {
        let codes = broadcast.predecessor_codes.get(destination);
        for (predecessor, code) in [
            (State::Aligned, 1),
            (State::Insertion, 2),
            (State::Deletion, 3),
        ] {
            assert_eq!(*codes.get(predecessor), [code << shift; TestI16::LANES]);
        }
    }
    Ok(())
}

/// The longest sequence compared exhaustively against the scalar engine.
///
/// Miri interprets every alignment, so it covers a smaller (but still
/// multi-diagonal) range of shapes.
const EXHAUSTIVE_LENGTH: usize = if cfg!(miri) { 3 } else { 5 };

/// Returns every binary byte sequence through `maximum_length`.
fn binary_sequences(maximum_length: usize) -> Vec<Vec<u8>> {
    let mut sequences = Vec::new();

    for length in 0..=maximum_length {
        for bits in 0..(1usize << length) {
            sequences.push((0..length).map(|bit| ((bits >> bit) & 1) as u8).collect());
        }
    }

    sequences
}

/// Returns the scoring fixtures shared by exhaustive differential tests.
fn exhaustive_scorings() -> Result<[Scoring; 4], Error> {
    Ok([
        Scoring::try_new(2, -3, -2, -1)?,
        Scoring::try_new(1, -5, -1, -1)?,
        Scoring::try_new(1, 0, 0, 0)?,
        Scoring::try_new(1, -1, 0, -1)?,
    ])
}

/// Compares the generic global kernel result with the scalar reference.
fn differential_global<K: Kernel>(kernel: K, maximum_length: usize) -> Result<(), Error> {
    let sequences = binary_sequences(maximum_length);

    for scoring in exhaustive_scorings()? {
        for reference in &sequences {
            for query in &sequences {
                let actual = global(kernel, reference, query, scoring);
                let expected = engine::global(reference, query, scoring);

                match (actual, expected) {
                    (Ok(actual), Ok(expected)) => assert_eq!(
                        actual, expected,
                        "global mismatch; reference={reference:?} query={query:?} \
                         scoring={scoring:?}"
                    ),
                    (Err(Error::EmptyGlobal), Err(Error::EmptyGlobal)) => {}
                    (actual, expected) => panic!(
                        "global result mismatch; reference={reference:?} query={query:?} \
                         scoring={scoring:?} actual={actual:?} expected={expected:?}"
                    ),
                }
            }
        }
    }

    Ok(())
}

/// Compares the generic local kernel result with the scalar reference.
fn differential_local<K: Kernel>(kernel: K, maximum_length: usize) -> Result<(), Error> {
    let sequences = binary_sequences(maximum_length);

    for scoring in exhaustive_scorings()? {
        for reference in &sequences {
            for query in &sequences {
                let actual = local(kernel, reference, query, scoring);
                let expected = engine::local(reference, query, scoring);

                match (actual, expected) {
                    (Ok(actual), Ok(expected)) => assert_eq!(
                        actual, expected,
                        "local mismatch; reference={reference:?} query={query:?} \
                         scoring={scoring:?}"
                    ),
                    (actual, expected) => panic!(
                        "local result mismatch; reference={reference:?} query={query:?} \
                         scoring={scoring:?} actual={actual:?} expected={expected:?}"
                    ),
                }
            }
        }
    }

    Ok(())
}

#[derive(Clone, Copy)]
struct DifferentialCase<'a> {
    name: &'static str,
    reference: &'a [u8],
    query: &'a [u8],
    full_vectors: usize,
    partial_lanes: usize,
}

/// Returns the deterministic `TestI16` cases used to reach vector paths.
fn vector_differential_cases() -> [DifferentialCase<'static>; 3] {
    [
        DifferentialCase {
            name: "vector-only",
            reference: b"ACGTACGT",
            query: b"AGCATGCT",
            full_vectors: 1,
            partial_lanes: 0,
        },
        DifferentialCase {
            name: "vector-plus-partial-vector",
            reference: b"CGTACGTAC",
            query: b"CGTACGTAC",
            full_vectors: 1,
            partial_lanes: 1,
        },
        DifferentialCase {
            name: "two-full-vectors",
            reference: b"TGCATGCATGCATGCA",
            query: b"TGCATGCATGCATGCA",
            full_vectors: 2,
            partial_lanes: 0,
        },
    ]
}

/// Returns the number of full `K` vectors and the number of valid lanes
/// in the final partial vector on the longest diagonal.
fn vector_path_counts<K: Kernel>(reference_len: usize, query_len: usize) -> (usize, usize) {
    let maximum_interior = reference_len.min(query_len);
    (maximum_interior / K::LANES, maximum_interior % K::LANES)
}

fn assert_vector_case_shape(case: DifferentialCase<'_>) {
    assert_eq!(
        vector_path_counts::<TestI16>(case.reference.len(), case.query.len()),
        (case.full_vectors, case.partial_lanes),
        "{} lengths should imply {} full vectors and {} partial-vector lanes",
        case.name,
        case.full_vectors,
        case.partial_lanes,
    );
}

fn differential_global_cases<K: Kernel>(
    kernel: K,
    cases: &[DifferentialCase<'_>],
) -> Result<(), Error> {
    for case in cases {
        assert_vector_case_shape(*case);
    }

    for scoring in exhaustive_scorings()? {
        for case in cases {
            assert_eq!(
                global(kernel, case.reference, case.query, scoring)?,
                engine::global(case.reference, case.query, scoring)?,
                "global mismatch; case={} scoring={scoring:?}",
                case.name,
            );
        }
    }

    Ok(())
}

fn differential_local_cases<K: Kernel>(
    kernel: K,
    cases: &[DifferentialCase<'_>],
) -> Result<(), Error> {
    for case in cases {
        assert_vector_case_shape(*case);
    }

    for scoring in exhaustive_scorings()? {
        for case in cases {
            assert_eq!(
                local(kernel, case.reference, case.query, scoring)?,
                engine::local(case.reference, case.query, scoring)?,
                "local mismatch; case={} scoring={scoring:?}",
                case.name,
            );
        }
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
/// Compares both modes for an input that reaches a complete kernel vector.
fn differential_full_vector_case<K: Kernel>(
    kernel: K,
    name: &str,
    reference: &[u8],
    query: &[u8],
) -> Result<(), Error> {
    assert!(
        reference.len().min(query.len()).saturating_sub(1) >= K::LANES,
        "{name} does not reach a complete {}-lane vector",
        K::LANES,
    );

    for scoring in exhaustive_scorings()? {
        assert_eq!(
            global(kernel, reference, query, scoring)?,
            engine::global(reference, query, scoring)?,
            "global mismatch; case={name} scoring={scoring:?}",
        );
        assert_eq!(
            local(kernel, reference, query, scoring)?,
            engine::local(reference, query, scoring)?,
            "local mismatch; case={name} scoring={scoring:?}",
        );
    }

    Ok(())
}

#[test]
fn test_vector_differential_cases_cover_test_i16_paths() {
    assert_eq!(vector_path_counts::<TestI16>(5, 5), (0, 5));

    for case in vector_differential_cases() {
        assert_vector_case_shape(case);
    }
}

#[test]
fn test_vector_differential_cases_include_mixed_substitutions() {
    let case = vector_differential_cases()[0];
    let mut comparisons = case.reference.iter().zip(case.query.iter().rev());

    assert!(
        comparisons
            .clone()
            .any(|(reference, query)| reference == query)
    );
    assert!(comparisons.any(|(reference, query)| reference != query));
}

#[test]
fn test_kernel_matches_scalar_for_global_vector_paths() -> Result<(), Error> {
    differential_global_cases(TEST_I16, &vector_differential_cases())
}

#[test]
fn test_kernel_matches_scalar_for_local_vector_paths() -> Result<(), Error> {
    differential_local_cases(TEST_I16, &vector_differential_cases())
}

#[test]
fn test_kernel_falls_back_outside_the_narrow_range() -> Result<(), Error> {
    let cases = [
        (
            b"A".as_slice(),
            b"A".as_slice(),
            Scoring::try_new(32_768, -1, -1, -1)?,
        ),
        (
            b"AAAAAAAA".as_slice(),
            b"AAAAAAAA".as_slice(),
            Scoring::try_new(5_000, -1, -1, -1)?,
        ),
    ];

    for (reference, query, scoring) in cases {
        assert_eq!(
            global(TEST_I16, reference, query, scoring)?,
            engine::global(reference, query, scoring)?
        );
        assert_eq!(
            local(TEST_I16, reference, query, scoring)?,
            engine::local(reference, query, scoring)?
        );
    }

    Ok(())
}

#[test]
fn test_kernel_matches_scalar_for_short_global_inputs() -> Result<(), Error> {
    differential_global(TEST_I16, EXHAUSTIVE_LENGTH)
}

#[test]
fn test_kernel_matches_scalar_for_short_local_inputs() -> Result<(), Error> {
    differential_local(TEST_I16, EXHAUSTIVE_LENGTH)
}

#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
#[test]
fn neon_matches_scalar() -> Result<(), Error> {
    use crate::algorithm::simd::aarch64::NeonI16;
    use crate::algorithm::simd::aarch64::NeonI32;

    differential_global(NeonI16, 5)?;
    differential_local(NeonI16, 5)?;
    differential_global_cases(NeonI16, &vector_differential_cases())?;
    differential_local_cases(NeonI16, &vector_differential_cases())?;
    differential_global(NeonI32, 5)?;
    differential_local(NeonI32, 5)?;

    assert_packed_codes_match_scalar(NeonI16)?;
    assert_local_fixtures(NeonI16)?;
    assert_every_remainder_matches_scalar(NeonI16)?;
    assert_packed_codes_match_scalar(NeonI32)?;
    assert_local_fixtures(NeonI32)?;
    assert_every_remainder_matches_scalar(NeonI32)
}

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
#[test]
fn avx2_matches_scalar_when_available() -> Result<(), Error> {
    use crate::algorithm::simd::x86_64::Avx2I16;
    use crate::algorithm::simd::x86_64::Avx2I32;
    use crate::algorithm::simd::x86_64::avx2_for_tests;

    let Some(avx2) = avx2_for_tests() else {
        return Ok(());
    };

    differential_global(Avx2I16(avx2), 5)?;
    differential_local(Avx2I16(avx2), 5)?;
    differential_global_cases(Avx2I16(avx2), &vector_differential_cases())?;
    differential_local_cases(Avx2I16(avx2), &vector_differential_cases())?;
    differential_global(Avx2I32(avx2), 5)?;
    differential_local(Avx2I32(avx2), 5)?;
    differential_global_cases(Avx2I32(avx2), &vector_differential_cases())?;
    differential_local_cases(Avx2I32(avx2), &vector_differential_cases())?;
    differential_full_vector_case(
        Avx2I16(avx2),
        "avx2-i16-mixed-substitutions",
        b"AAAAAAAAAAAAAAAAA",
        b"AAAAAAAABBBBBBBBB",
    )?;
    assert_packed_codes_match_scalar(Avx2I16(avx2))?;
    assert_local_fixtures(Avx2I16(avx2))?;
    assert_every_remainder_matches_scalar(Avx2I16(avx2))?;
    assert_packed_codes_match_scalar(Avx2I32(avx2))?;
    assert_local_fixtures(Avx2I32(avx2))?;
    assert_every_remainder_matches_scalar(Avx2I32(avx2))
}

/// Checks that adding one field code per destination state yields the
/// scalar packed byte for all 64 predecessor combinations, including the
/// "no predecessor" code that a local reset leaves in the aligned field.
fn assert_packed_codes_match_scalar<K: Kernel>(kernel: K) -> Result<(), Error> {
    let lane_scoring = LaneScoring::<K::Score>::try_new(Scoring::try_new(2, -3, -2, -1)?)?;
    let broadcast = BroadcastScoring::new(kernel, &lane_scoring);
    let predecessors = [
        None,
        Some(State::Aligned),
        Some(State::Insertion),
        Some(State::Deletion),
    ];
    let code = |destination: State, predecessor: Option<State>| {
        predecessor.map_or(broadcast.zero, |predecessor| {
            *broadcast
                .predecessor_codes
                .get(destination)
                .get(predecessor)
        })
    };
    let mut packed = allocated_scores::<K::Score>(K::LANES)?;

    for aligned in predecessors {
        for insertion in predecessors {
            for deletion in predecessors {
                let vector = kernel.add(
                    kernel.add(
                        code(State::Aligned, aligned),
                        code(State::Insertion, insertion),
                    ),
                    code(State::Deletion, deletion),
                );
                // SAFETY: `packed` holds exactly `K::LANES` scores.
                unsafe { kernel.store(packed.as_mut_ptr(), vector) };
                let expected = pack_predecessors(aligned, insertion, deletion);

                for lane in &packed {
                    let byte = lane.into_code();
                    assert_eq!(byte, expected);
                    assert_eq!(unpack_predecessor(byte, State::Aligned), aligned);
                    assert_eq!(unpack_predecessor(byte, State::Insertion), insertion);
                    assert_eq!(unpack_predecessor(byte, State::Deletion), deletion);
                }
            }
        }
    }

    Ok(())
}

/// Checks one local fixture against its expected outcome and the scalar
/// engine.
fn assert_local_fixture<K: Kernel>(
    kernel: K,
    (reference, query): (&[u8], &[u8]),
    scoring: Scoring,
    (score, cigar, reference_range, query_range): (Score, &str, Range<usize>, Range<usize>),
) -> Result<(), Error> {
    let actual = local(kernel, reference, query, scoring)?
        .unwrap_or_else(|| panic!("no local alignment for {reference:?} and {query:?}"));

    assert_eq!(actual.score(), score);
    assert_eq!(actual.cigar().to_string(), cigar);
    assert_eq!(actual.reference_range(), &reference_range);
    assert_eq!(actual.query_range(), &query_range);
    assert_eq!(Some(actual), engine::local(reference, query, scoring)?);
    Ok(())
}

/// Checks the local fixtures for one kernel.
fn assert_local_fixtures<K: Kernel>(kernel: K) -> Result<(), Error> {
    // The insertion's aligned predecessor is packed as `1 << 2`, and the
    // aligned start cell has no predecessor.
    assert_local_fixture(
        kernel,
        (b"AB", b"ACB"),
        Scoring::try_new(2, -3, -1, -1)?,
        (3, "1=1I1=", 0..2, 0..3),
    )?;
    // Two single matches tie; the first in row-major order of exclusive
    // `(reference_end, query_end)` wins: (1, 4) before (3, 1).
    assert_local_fixture(
        kernel,
        (b"BXA", b"AYYB"),
        Scoring::try_new(2, -3, -2, -1)?,
        (2, "1=", 0..1, 3..4),
    )
}

/// Returns deterministic DNA-like bytes with mixed matches and
/// mismatches.
fn mixed_bytes(length: usize, salt: usize) -> Vec<u8> {
    (0..length)
        .map(|index| b"ACGT"[(index * 7 + salt + index / 3) % 4])
        .collect()
}

/// Compares both modes with the scalar engine for every number of valid
/// lanes in the final vector of the longest diagonal, with the shorter
/// input on either axis.
fn assert_every_remainder_matches_scalar<K: Kernel>(kernel: K) -> Result<(), Error> {
    let full_vector_counts = if cfg!(miri) { 1..=1 } else { 1..=2 };

    for scoring in exhaustive_scorings()? {
        for remainder in 0..K::LANES {
            for full_vectors in full_vector_counts.clone() {
                let short = full_vectors * K::LANES + remainder;
                let long = short + 3;

                for (reference_len, query_len) in [(short, long), (long, short)] {
                    let reference = mixed_bytes(reference_len, 0);
                    let query = mixed_bytes(query_len, 1);
                    assert_eq!(
                        vector_path_counts::<K>(reference_len, query_len),
                        (full_vectors, remainder)
                    );
                    assert_eq!(
                        global(kernel, &reference, &query, scoring)?,
                        engine::global(&reference, &query, scoring)?,
                        "global mismatch; lengths=({reference_len}, {query_len}) \
                         scoring={scoring:?}"
                    );
                    assert_eq!(
                        local(kernel, &reference, &query, scoring)?,
                        engine::local(&reference, &query, scoring)?,
                        "local mismatch; lengths=({reference_len}, {query_len}) \
                         scoring={scoring:?}"
                    );
                }
            }
        }
    }

    Ok(())
}

#[test]
fn test_kernels_pack_every_predecessor_combination() -> Result<(), Error> {
    assert_packed_codes_match_scalar(TEST_I16)?;
    assert_packed_codes_match_scalar(TEST_I32)
}

#[test]
fn test_kernels_match_local_fixtures() -> Result<(), Error> {
    assert_local_fixtures(TEST_I16)?;
    assert_local_fixtures(TEST_I32)
}

#[test]
fn test_kernels_match_scalar_for_every_remainder() -> Result<(), Error> {
    assert_every_remainder_matches_scalar(TEST_I16)?;
    assert_every_remainder_matches_scalar(TEST_I32)
}

#[test]
fn test_i32_kernel_matches_scalar_for_short_inputs() -> Result<(), Error> {
    differential_global(TEST_I32, EXHAUSTIVE_LENGTH)?;
    differential_local(TEST_I32, EXHAUSTIVE_LENGTH)
}

#[test]
fn selects_the_narrowest_proven_lane() -> Result<(), Error> {
    let scoring = Scoring::try_new(2, -3, -2, -1)?;
    assert_eq!(lane_width(100, 100, scoring), LaneWidth::I16);
    assert_eq!(lane_width(10_000, 10_000, scoring), LaneWidth::I32);

    let extreme = Scoring::try_new(1, -1, Score::MIN, -1)?;
    assert_eq!(lane_width(1, 1, extreme), LaneWidth::Scalar);
    Ok(())
}

#[test]
fn diagonal_layout_round_trips_coordinates() -> Result<(), Error> {
    let layout = Layout::try_new(4, 3)?;

    for row in 0..4 {
        for column in 0..3 {
            let index = layout.index(row, column);
            assert_eq!(layout.coordinate(index), (row, column));
        }
    }

    Ok(())
}

#[test]
fn diagonal_layout_reports_dimension_overflow() {
    assert!(matches!(
        Layout::try_new(usize::MAX, 2),
        Err(Error::MatrixSizeOverflow)
    ));
}

#[test]
fn diagonal_layout_rejects_zero_sized_dimensions() {
    assert!(matches!(
        Layout::try_new(0, 1),
        Err(Error::MatrixSizeOverflow)
    ));
    assert!(matches!(
        Layout::try_new(1, 0),
        Err(Error::MatrixSizeOverflow)
    ));
}

#[test]
fn diagonal_layout_reports_allocation_failure() {
    let error = match Layout::try_new(1, usize::MAX - 1) {
        Err(error) => error,
        Ok(_) => panic!("expected allocation failure"),
    };

    assert!(matches!(error, Error::MatrixAllocation { .. }));
    assert_eq!(error.to_string(), "failed to allocate alignment storage");
}
