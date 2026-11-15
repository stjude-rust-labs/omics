//! Tests that force each native SIMD backend and compare its results with the
//! scalar engine.
//!
//! The regression cases cover both lane widths, full and partial vectors,
//! mixed substitutions within one vector, and inputs that fall back to the
//! scalar engine.

use std::fmt::Debug;

use super::wavefront::LaneWidth;
use super::wavefront::lane_width;
use crate::algorithm::Error;
use crate::algorithm::Outcome;
use crate::algorithm::Score;
use crate::algorithm::Scoring;
use crate::algorithm::engine;

/// A forced global backend under test.
type GlobalFn<'a> = dyn Fn(&[u8], &[u8], Scoring) -> Result<Outcome, Error> + 'a;

/// A forced local backend under test.
type LocalFn<'a> = dyn Fn(&[u8], &[u8], Scoring) -> Result<Option<Outcome>, Error> + 'a;

/// One named input pair and scoring scheme.
#[derive(Clone, Copy)]
struct BackendCase<'a> {
    /// The name reported when the case fails.
    name: &'static str,
    /// The reference sequence.
    reference: &'a [u8],
    /// The query sequence.
    query: &'a [u8],
    /// The scoring scheme.
    scoring: Scoring,
}

impl<'a> BackendCase<'a> {
    /// Creates a named case.
    fn new(name: &'static str, reference: &'a [u8], query: &'a [u8], scoring: Scoring) -> Self {
        Self {
            name,
            reference,
            query,
            scoring,
        }
    }
}

/// Returns the regression cases shared by every backend, followed by
/// `i32_case`, a backend-specific case that fills whole `i32` vectors.
fn backend_regression_cases(
    i32_case: BackendCase<'static>,
) -> Result<[BackendCase<'static>; 5], Error> {
    Ok([
        BackendCase::new("empty-both", b"", b"", Scoring::try_new(2, -3, -2, -1)?),
        BackendCase::new(
            "zero-cost-gaps",
            b"AAAA",
            b"AA",
            Scoring::try_new(1, -1, 0, 0)?,
        ),
        BackendCase::new(
            "falls-back-for-score-min",
            b"AC",
            b"A",
            Scoring::try_new(1, -1, Score::MIN, -1)?,
        ),
        BackendCase::new(
            "selects-i16-lanes",
            b"ACGTACGT",
            b"ACGTACGT",
            Scoring::try_new(2, -3, -2, -1)?,
        ),
        i32_case,
    ])
}

/// Returns the number of full vectors and the number of partial-vector lanes
/// on the longest anti-diagonal.
fn vector_path_counts(lanes: usize, reference_len: usize, query_len: usize) -> (usize, usize) {
    let maximum_interior = reference_len.min(query_len);
    (maximum_interior / lanes, maximum_interior % lanes)
}

/// Checks that `case` exercises the intended lane width and vector shape.
fn assert_backend_case_shape(
    case: BackendCase<'_>,
    lanes: usize,
    expected_lane_width: LaneWidth,
    full_vectors: usize,
    partial_lanes: usize,
) {
    assert_eq!(
        lane_width(case.reference.len(), case.query.len(), case.scoring),
        expected_lane_width,
        "{} should select {expected_lane_width:?}",
        case.name,
    );
    assert_eq!(
        vector_path_counts(lanes, case.reference.len(), case.query.len()),
        (full_vectors, partial_lanes),
        "{} should imply {full_vectors} full vectors and {partial_lanes} partial-vector lanes",
        case.name,
    );
    assert_mixed_substitutions_per_vector(case, lanes, full_vectors);
}

/// Checks that every full vector on the longest anti-diagonal compares both
/// equal and unequal byte pairs.
fn assert_mixed_substitutions_per_vector(case: BackendCase<'_>, lanes: usize, full_vectors: usize) {
    let comparisons: Vec<bool> = case
        .reference
        .iter()
        .zip(case.query.iter().rev())
        .map(|(reference, query)| reference == query)
        .collect();

    for (chunk_index, chunk) in comparisons
        .chunks_exact(lanes)
        .take(full_vectors)
        .enumerate()
    {
        assert!(
            chunk.iter().any(|&is_equal| is_equal),
            "{} vector chunk {} should include equal anti-diagonal substitutions",
            case.name,
            chunk_index,
        );
        assert!(
            chunk.iter().any(|&is_equal| !is_equal),
            "{} vector chunk {} should include unequal anti-diagonal substitutions",
            case.name,
            chunk_index,
        );
    }
}

/// Checks that both forced backends match the scalar engine on `case`.
fn assert_backend_matches_scalar_case(
    case: BackendCase<'_>,
    forced_global: &GlobalFn<'_>,
    forced_local: &LocalFn<'_>,
) {
    assert_same_result(
        "global",
        case.name,
        forced_global(case.reference, case.query, case.scoring),
        engine::global(case.reference, case.query, case.scoring),
    );
    assert_same_result(
        "local",
        case.name,
        forced_local(case.reference, case.query, case.scoring),
        engine::local(case.reference, case.query, case.scoring),
    );
}

/// Checks that both forced backends match the scalar engine on every case.
fn assert_backend_matches_scalar_cases(
    cases: &[BackendCase<'_>],
    forced_global: &GlobalFn<'_>,
    forced_local: &LocalFn<'_>,
) {
    for &case in cases {
        assert_backend_matches_scalar_case(case, forced_global, forced_local);
    }
}

/// Checks that two results are equal, including the error variant and its
/// message.
fn assert_same_result<T: Debug + PartialEq>(
    mode: &str,
    case_name: &str,
    actual: Result<T, Error>,
    expected: Result<T, Error>,
) {
    match (actual, expected) {
        (Ok(actual), Ok(expected)) => {
            assert_eq!(actual, expected, "forced {mode} mismatch; case={case_name}")
        }
        (Err(actual), Err(expected)) => {
            assert_eq!(
                actual.to_string(),
                expected.to_string(),
                "forced {mode} error mismatch; case={case_name}",
            );
            assert_eq!(
                format!("{actual:?}"),
                format!("{expected:?}"),
                "forced {mode} error variant mismatch; case={case_name}",
            );
        }
        (actual, expected) => panic!(
            "forced {mode} result mismatch; case={case_name} actual={actual:?} \
             expected={expected:?}"
        ),
    }
}

/// Returns a case whose `i32` lanes fill two full NEON vectors with both
/// equal and unequal byte pairs.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
fn neon_mixed_substitution_multivector_case() -> Result<BackendCase<'static>, Error> {
    Ok(BackendCase::new(
        "neon-i32-mixed-substitutions",
        b"AAAABBBB",
        b"DBDBCACA",
        Scoring::try_new(5_000, -1, -1, -1)?,
    ))
}

#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
#[test]
fn neon_dispatch_matches_scalar() -> Result<(), Error> {
    use super::aarch64;

    assert_backend_matches_scalar_cases(
        &backend_regression_cases(neon_mixed_substitution_multivector_case()?)?,
        &aarch64::global,
        &aarch64::local,
    );
    Ok(())
}

#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
#[test]
fn neon_i32_mixed_substitutions_reach_two_full_vectors() -> Result<(), Error> {
    use super::aarch64;

    let case = neon_mixed_substitution_multivector_case()?;
    assert_backend_case_shape(case, 4, LaneWidth::I32, 2, 0);
    assert_backend_matches_scalar_case(case, &aarch64::global, &aarch64::local);
    Ok(())
}

/// Returns a case whose `i32` lanes fill two full AVX2 vectors with both
/// equal and unequal byte pairs.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
fn avx2_mixed_substitution_multivector_case() -> Result<BackendCase<'static>, Error> {
    Ok(BackendCase::new(
        "avx2-i32-mixed-substitutions",
        b"AAAAAAAABBBBBBBB",
        b"DBDBDBDBCACACACA",
        Scoring::try_new(5_000, -1, -1, -1)?,
    ))
}

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
#[test]
fn avx2_dispatch_matches_scalar() -> Result<(), Error> {
    use super::x86_64;

    let Some(avx2) = x86_64::avx2_for_tests() else {
        return Ok(());
    };
    assert_backend_matches_scalar_cases(
        &backend_regression_cases(avx2_mixed_substitution_multivector_case()?)?,
        &|reference, query, scoring| x86_64::global(avx2, reference, query, scoring),
        &|reference, query, scoring| x86_64::local(avx2, reference, query, scoring),
    );
    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
#[test]
fn avx2_i32_mixed_substitutions_reach_two_full_vectors() -> Result<(), Error> {
    use super::x86_64;

    let Some(avx2) = x86_64::avx2_for_tests() else {
        return Ok(());
    };
    let case = avx2_mixed_substitution_multivector_case()?;
    assert_backend_case_shape(case, 8, LaneWidth::I32, 2, 0);
    assert_backend_matches_scalar_case(
        case,
        &|reference, query, scoring| x86_64::global(avx2, reference, query, scoring),
        &|reference, query, scoring| x86_64::local(avx2, reference, query, scoring),
    );
    Ok(())
}
