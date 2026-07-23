//! SIMD-dispatched affine-gap alignment for byte slices.
//!
//! The public API retains the scalar fallback until a target-specific kernel
//! becomes available. The private wavefront engine supports backend
//! differential tests without changing byte-API results.

use super::Error;
use super::Outcome;
#[cfg(test)]
use super::Score;
use super::Scoring;
use super::engine;

/// Apple Silicon NEON wavefront kernels and dispatch.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
mod aarch64;

/// Linux x86_64 AVX2 wavefront kernels and dispatch.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
mod x86_64;

/// Shared private wavefront helpers for future SIMD kernels.
mod wavefront;

/// Returns the global alignment for byte slices.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    aarch64::global(reference, query, scoring)
}

/// Returns the global alignment for byte slices.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    x86_64::global(reference, query, scoring)
}

/// Returns the global alignment for byte slices.
#[cfg(not(any(
    all(target_arch = "aarch64", target_os = "macos"),
    all(target_arch = "x86_64", target_os = "linux")
)))]
pub fn global(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Outcome, Error> {
    engine::global(reference, query, scoring)
}

/// Returns the local alignment for byte slices.
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    aarch64::local(reference, query, scoring)
}

/// Returns the local alignment for byte slices.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    x86_64::local(reference, query, scoring)
}

/// Returns the local alignment for byte slices.
#[cfg(not(any(
    all(target_arch = "aarch64", target_os = "macos"),
    all(target_arch = "x86_64", target_os = "linux")
)))]
pub fn local(reference: &[u8], query: &[u8], scoring: Scoring) -> Result<Option<Outcome>, Error> {
    engine::local(reference, query, scoring)
}

#[cfg(test)]
#[derive(Clone, Copy)]
struct BackendCase<'a> {
    name: &'static str,
    reference: &'a [u8],
    query: &'a [u8],
    scoring: Scoring,
}

#[cfg(test)]
impl<'a> BackendCase<'a> {
    fn new(name: &'static str, reference: &'a [u8], query: &'a [u8], scoring: Scoring) -> Self {
        Self {
            name,
            reference,
            query,
            scoring,
        }
    }
}

#[cfg(test)]
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

#[cfg(all(test, target_arch = "aarch64", target_os = "macos"))]
fn neon_mixed_substitution_multivector_case() -> Result<BackendCase<'static>, Error> {
    Ok(BackendCase::new(
        "neon-i32-mixed-substitutions",
        b"AAAABBBB",
        b"DBDBCACA",
        Scoring::try_new(5_000, -1, -1, -1)?,
    ))
}

#[cfg(all(test, target_arch = "x86_64", target_os = "linux"))]
fn avx2_mixed_substitution_multivector_case() -> Result<BackendCase<'static>, Error> {
    Ok(BackendCase::new(
        "avx2-i32-mixed-substitutions",
        b"AAAAAAAABBBBBBBB",
        b"DBDBDBDBCACACACA",
        Scoring::try_new(5_000, -1, -1, -1)?,
    ))
}

#[cfg(test)]
fn vector_path_counts(lanes: usize, reference_len: usize, query_len: usize) -> (usize, usize) {
    let maximum_interior = reference_len.min(query_len);
    (maximum_interior / lanes, maximum_interior % lanes)
}

#[cfg(test)]
fn assert_backend_case_shape(
    case: BackendCase<'_>,
    lanes: usize,
    lane_width: wavefront::LaneWidth,
    full_vectors: usize,
    tail_cells: usize,
) {
    assert_eq!(
        wavefront::lane_width(case.reference.len(), case.query.len(), case.scoring),
        lane_width,
        "{} should select {lane_width:?}",
        case.name,
    );
    assert_eq!(
        vector_path_counts(lanes, case.reference.len(), case.query.len()),
        (full_vectors, tail_cells),
        "{} should imply {full_vectors} full vectors and {tail_cells} tail cells",
        case.name,
    );
    assert_mixed_substitutions_per_vector(case, lanes, full_vectors);
}

#[cfg(test)]
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

#[cfg(test)]
fn assert_backend_matches_scalar_case(
    case: BackendCase<'_>,
    forced_global: fn(&[u8], &[u8], Scoring) -> Result<Outcome, Error>,
    forced_local: fn(&[u8], &[u8], Scoring) -> Result<Option<Outcome>, Error>,
) -> Result<(), Error> {
    assert_same_global_result(
        forced_global(case.reference, case.query, case.scoring),
        engine::global(case.reference, case.query, case.scoring),
        case.name,
    );
    assert_same_local_result(
        forced_local(case.reference, case.query, case.scoring),
        engine::local(case.reference, case.query, case.scoring),
        case.name,
    );

    Ok(())
}

#[cfg(test)]
fn assert_same_global_result(
    actual: Result<Outcome, Error>,
    expected: Result<Outcome, Error>,
    case_name: &str,
) {
    match (actual, expected) {
        (Ok(actual), Ok(expected)) => {
            assert_eq!(actual, expected, "forced global mismatch; case={case_name}",)
        }
        (Err(actual), Err(expected)) => {
            assert_eq!(
                actual.to_string(),
                expected.to_string(),
                "forced global error mismatch; case={case_name}",
            );
            assert_eq!(
                format!("{actual:?}"),
                format!("{expected:?}"),
                "forced global error variant mismatch; case={case_name}",
            );
        }
        (actual, expected) => panic!(
            "forced global result mismatch; case={case_name} actual={actual:?} \
             expected={expected:?}"
        ),
    }
}

#[cfg(test)]
fn assert_same_local_result(
    actual: Result<Option<Outcome>, Error>,
    expected: Result<Option<Outcome>, Error>,
    case_name: &str,
) {
    match (actual, expected) {
        (Ok(actual), Ok(expected)) => {
            assert_eq!(actual, expected, "forced local mismatch; case={case_name}",)
        }
        (Err(actual), Err(expected)) => {
            assert_eq!(
                actual.to_string(),
                expected.to_string(),
                "forced local error mismatch; case={case_name}",
            );
            assert_eq!(
                format!("{actual:?}"),
                format!("{expected:?}"),
                "forced local error variant mismatch; case={case_name}",
            );
        }
        (actual, expected) => panic!(
            "forced local result mismatch; case={case_name} actual={actual:?} \
             expected={expected:?}"
        ),
    }
}

#[cfg(test)]
fn assert_backend_matches_scalar_cases(
    cases: &[BackendCase<'_>],
    forced_global: fn(&[u8], &[u8], Scoring) -> Result<Outcome, Error>,
    forced_local: fn(&[u8], &[u8], Scoring) -> Result<Option<Outcome>, Error>,
) -> Result<(), Error> {
    for &case in cases {
        assert_backend_matches_scalar_case(case, forced_global, forced_local)?;
    }

    Ok(())
}

#[cfg(all(test, target_arch = "aarch64", target_os = "macos"))]
mod tests {
    use super::*;

    #[test]
    fn neon_dispatch_matches_scalar() -> Result<(), Error> {
        assert_backend_matches_scalar_cases(
            &backend_regression_cases(neon_mixed_substitution_multivector_case()?)?,
            aarch64::global,
            aarch64::local,
        )
    }

    #[test]
    fn neon_i32_mixed_substitutions_reach_two_full_vectors() -> Result<(), Error> {
        let case = neon_mixed_substitution_multivector_case()?;
        assert_backend_case_shape(case, 4, wavefront::LaneWidth::I32, 2, 0);
        assert_backend_matches_scalar_case(case, aarch64::global, aarch64::local)
    }
}

#[cfg(all(test, target_arch = "x86_64", target_os = "linux"))]
mod x86_64_tests {
    use super::*;

    #[test]
    fn avx2_dispatch_matches_scalar() -> Result<(), Error> {
        assert_backend_matches_scalar_cases(
            &backend_regression_cases(avx2_mixed_substitution_multivector_case()?)?,
            x86_64::global,
            x86_64::local,
        )
    }

    #[test]
    fn avx2_i32_mixed_substitutions_reach_two_full_vectors() -> Result<(), Error> {
        let case = avx2_mixed_substitution_multivector_case()?;
        assert_backend_case_shape(case, 8, wavefront::LaneWidth::I32, 2, 0);
        assert_backend_matches_scalar_case(case, x86_64::global, x86_64::local)
    }
}
