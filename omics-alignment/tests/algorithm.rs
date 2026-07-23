//! Exhaustive oracle and public compatibility tests for the alignment
//! algorithms.

use omics_alignment::algorithm::Outcome;
use omics_alignment::algorithm::Score;
use omics_alignment::algorithm::Scoring;
use omics_alignment::algorithm::global;
use omics_alignment::algorithm::local;
use omics_alignment::cigar::OperationKind;
use omics_molecule::polymer::dna::Nucleotide;

/// A move in an exhaustively enumerated alignment path.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Move {
    /// Aligns one symbol from each input.
    Aligned,
    /// Consumes one query symbol without advancing the reference.
    Insertion,
    /// Consumes one reference symbol without advancing the query.
    Deletion,
}

/// Returns the affine gap score for a run of the given length.
///
/// Caller must ensure that `length` is at least `1` and that `length - 1`
/// fits in [`Score`].
fn gap_score(scoring: Scoring, length: usize) -> Score {
    // SAFETY: exhaustive callers provide maximal runs from inputs bounded by 3.
    let extensions = Score::try_from(length - 1).unwrap();
    scoring.gap_open_score() + extensions * scoring.gap_extend_score()
}

/// Scores a complete path against the given input slices.
///
/// Groups consecutive same-kind gap moves into maximal runs and
/// applies the affine formula to each run.
fn score_path<T: Eq>(reference: &[T], query: &[T], path: &[Move], scoring: Scoring) -> Score {
    let mut ri = 0;
    let mut qi = 0;
    let mut total: Score = 0;
    let mut i = 0;
    while i < path.len() {
        match path[i] {
            Move::Aligned => {
                total += if reference[ri] == query[qi] {
                    scoring.match_score()
                } else {
                    scoring.mismatch_score()
                };
                ri += 1;
                qi += 1;
                i += 1;
            }
            Move::Insertion => {
                let start = i;
                while i < path.len() && path[i] == Move::Insertion {
                    i += 1;
                }
                let run = i - start;
                total += gap_score(scoring, run);
                qi += run;
            }
            Move::Deletion => {
                let start = i;
                while i < path.len() && path[i] == Move::Deletion {
                    i += 1;
                }
                let run = i - start;
                total += gap_score(scoring, run);
                ri += run;
            }
        }
    }
    total
}

/// Recursively enumerates all complete alignment paths, updating the running
/// best.
fn enumerate_paths<T: Eq>(
    reference: &[T],
    query: &[T],
    ri: usize,
    qi: usize,
    path: &mut Vec<Move>,
    best: &mut Score,
    scoring: Scoring,
) {
    if ri == reference.len() && qi == query.len() {
        let s = score_path(reference, query, path, scoring);
        if s > *best {
            *best = s;
        }
        return;
    }
    if ri < reference.len() && qi < query.len() {
        path.push(Move::Aligned);
        enumerate_paths(reference, query, ri + 1, qi + 1, path, best, scoring);
        path.pop();
    }
    if qi < query.len() {
        path.push(Move::Insertion);
        enumerate_paths(reference, query, ri, qi + 1, path, best, scoring);
        path.pop();
    }
    if ri < reference.len() {
        path.push(Move::Deletion);
        enumerate_paths(reference, query, ri + 1, qi, path, best, scoring);
        path.pop();
    }
}

/// Returns the maximum score over all complete global alignment paths.
fn oracle_global_score<T: Eq>(reference: &[T], query: &[T], scoring: Scoring) -> Score {
    let mut best = Score::MIN;
    let mut path = Vec::new();
    enumerate_paths(reference, query, 0, 0, &mut path, &mut best, scoring);
    best
}

/// Enumerates complete paths, retaining only those whose first and last
/// moves are both Aligned.
fn enumerate_local_paths<T: Eq>(
    reference: &[T],
    query: &[T],
    ri: usize,
    qi: usize,
    path: &mut Vec<Move>,
    best: &mut Score,
    scoring: Scoring,
) {
    if ri == reference.len() && qi == query.len() {
        if matches!(path.first(), Some(Move::Aligned)) && matches!(path.last(), Some(Move::Aligned))
        {
            let s = score_path(reference, query, path, scoring);
            if s > *best {
                *best = s;
            }
        }
        return;
    }
    if ri < reference.len() && qi < query.len() {
        path.push(Move::Aligned);
        enumerate_local_paths(reference, query, ri + 1, qi + 1, path, best, scoring);
        path.pop();
    }
    if qi < query.len() {
        path.push(Move::Insertion);
        enumerate_local_paths(reference, query, ri, qi + 1, path, best, scoring);
        path.pop();
    }
    if ri < reference.len() {
        path.push(Move::Deletion);
        enumerate_local_paths(reference, query, ri + 1, qi, path, best, scoring);
        path.pop();
    }
}

/// Returns the maximum local alignment score over every non-empty range pair.
///
/// The result is clamped to zero; returns zero when no positive alignment
/// exists with aligned first and last moves.
fn oracle_local_score<T: Eq>(reference: &[T], query: &[T], scoring: Scoring) -> Score {
    let mut best: Score = 0;
    for rs in 0..reference.len() {
        for re in (rs + 1)..=reference.len() {
            for qs in 0..query.len() {
                for qe in (qs + 1)..=query.len() {
                    let mut range_best = Score::MIN;
                    let mut path = Vec::new();
                    enumerate_local_paths(
                        &reference[rs..re],
                        &query[qs..qe],
                        0,
                        0,
                        &mut path,
                        &mut range_best,
                        scoring,
                    );
                    if range_best > best {
                        best = range_best;
                    }
                }
            }
        }
    }
    best
}

/// Returns all binary sequences of length 0 through 3.
fn all_binary_seqs() -> Vec<Vec<u8>> {
    let mut seqs = Vec::new();
    for len in 0usize..=3 {
        for mask in 0u8..(1u8 << len) {
            let seq: Vec<u8> = (0..len).map(|i| (mask >> i) & 1).collect();
            seqs.push(seq);
        }
    }
    seqs
}

/// Asserts that an outcome's CIGAR is consistent with the inputs and scoring.
///
/// Checks symbol equality for matches and mismatches, reconstructs the score
/// via the affine formula, verifies consumed lengths equal both range lengths,
/// and verifies no two adjacent operations share a kind.
fn assert_outcome_consistent<T: Eq + std::fmt::Debug>(
    reference: &[T],
    query: &[T],
    scoring: Scoring,
    outcome: &Outcome,
) {
    let ref_slice = &reference[outcome.reference_range().clone()];
    let qry_slice = &query[outcome.query_range().clone()];
    let ops = outcome.cigar().operations();

    for w in ops.windows(2) {
        assert_ne!(w[0].kind(), w[1].kind(), "adjacent operations share a kind");
    }

    let mut ri = 0;
    let mut qi = 0;
    let mut reconstructed: Score = 0;

    for op in ops {
        let kind = op.kind();
        let len = op.length() as usize;

        match kind {
            OperationKind::SequenceMatch => {
                for k in 0..len {
                    assert_eq!(ref_slice[ri + k], qry_slice[qi + k], "= but symbols differ");
                }
                reconstructed += len as Score * scoring.match_score();
                ri += len;
                qi += len;
            }
            OperationKind::SequenceMismatch => {
                for k in 0..len {
                    assert_ne!(ref_slice[ri + k], qry_slice[qi + k], "X but symbols match");
                }
                reconstructed += len as Score * scoring.mismatch_score();
                ri += len;
                qi += len;
            }
            OperationKind::Insertion => {
                reconstructed += gap_score(scoring, len);
                qi += len;
            }
            OperationKind::Deletion => {
                reconstructed += gap_score(scoring, len);
                ri += len;
            }
            other => panic!("unexpected operation kind {other:?}"),
        }
    }

    assert_eq!(ri, ref_slice.len(), "reference consumption mismatch");
    assert_eq!(qi, qry_slice.len(), "query consumption mismatch");
    assert_eq!(
        reconstructed,
        outcome.score(),
        "reconstructed score mismatch"
    );
}

/// Scoring configurations used across all exhaustive tests.
fn exhaustive_scorings() -> [Scoring; 4] {
    [
        Scoring::try_new(2, -3, -2, -1).unwrap(),
        Scoring::try_new(1, -5, -1, -1).unwrap(),
        Scoring::try_new(1, 0, 0, 0).unwrap(),
        Scoring::try_new(1, -1, 0, -1).unwrap(),
    ]
}

#[test]
fn aligns_byte_slices_through_the_public_api() -> Result<(), Box<dyn std::error::Error>> {
    let scoring = Scoring::try_new(2, -3, -2, -1)?;
    assert_eq!(global(b"AC", b"AC", scoring)?.cigar().to_string(), "2=");
    Ok(())
}

#[test]
fn simd_public_api_matches_scalar_global() -> Result<(), Box<dyn std::error::Error>> {
    let scoring = Scoring::try_new(2, -3, -2, -1)?;
    assert_eq!(
        omics_alignment::algorithm::simd::global(b"ACGT", b"AGT", scoring)?,
        global(b"ACGT", b"AGT", scoring)?
    );
    Ok(())
}

#[test]
fn simd_public_api_matches_scalar_local() -> Result<(), Box<dyn std::error::Error>> {
    let scoring = Scoring::try_new(2, -3, -2, -1)?;
    assert_eq!(
        omics_alignment::algorithm::simd::local(b"GGACGT", b"TTACGA", scoring)?,
        local(b"GGACGT", b"TTACGA", scoring)?
    );
    Ok(())
}

#[test]
fn aligns_nucleotide_slices_without_a_production_dependency()
-> Result<(), Box<dyn std::error::Error>> {
    let reference = [Nucleotide::A, Nucleotide::C];
    let query = [Nucleotide::A, Nucleotide::G];
    let scoring = Scoring::try_new(2, -1, -2, -1)?;
    assert_eq!(
        global(&reference, &query, scoring)?.cigar().to_string(),
        "1=1X"
    );
    Ok(())
}

#[test]
fn exhaustive_global_score_matches_oracle() {
    let seqs = all_binary_seqs();
    for scoring in exhaustive_scorings() {
        for reference in &seqs {
            for query in &seqs {
                if reference.is_empty() && query.is_empty() {
                    continue;
                }
                let expected = oracle_global_score(reference, query, scoring);
                let outcome = global(reference, query, scoring).expect("global failed");
                assert_eq!(
                    outcome.score(),
                    expected,
                    "global score mismatch; ref={reference:?} query={query:?} scoring={scoring:?}"
                );
                assert_outcome_consistent(reference, query, scoring, &outcome);
            }
        }
    }
}

#[test]
fn exhaustive_local_score_matches_oracle() {
    let seqs = all_binary_seqs();
    for scoring in exhaustive_scorings() {
        for reference in &seqs {
            for query in &seqs {
                let expected = oracle_local_score(reference, query, scoring);
                let result = local(reference, query, scoring).expect("local failed");
                if expected == 0 {
                    assert!(
                        result.is_none(),
                        "expected None; ref={reference:?} query={query:?} scoring={scoring:?}"
                    );
                } else {
                    let outcome = result.unwrap_or_else(|| {
                        panic!(
                            "expected Some; ref={reference:?} query={query:?} scoring={scoring:?} \
                             oracle={expected}"
                        )
                    });
                    assert_eq!(
                        outcome.score(),
                        expected,
                        "local score mismatch; ref={reference:?} query={query:?} \
                         scoring={scoring:?}"
                    );
                    let ops = outcome.cigar().operations();
                    let first = ops.first().expect("empty CIGAR").kind();
                    let last = ops.last().expect("empty CIGAR").kind();
                    let is_aligned = |k: OperationKind| {
                        k == OperationKind::SequenceMatch || k == OperationKind::SequenceMismatch
                    };
                    assert!(
                        is_aligned(first),
                        "local CIGAR does not start with an aligned op; ref={reference:?} \
                         query={query:?}"
                    );
                    assert!(
                        is_aligned(last),
                        "local CIGAR does not end with an aligned op; ref={reference:?} \
                         query={query:?}"
                    );
                    assert_outcome_consistent(reference, query, scoring, &outcome);
                }
            }
        }
    }
}
