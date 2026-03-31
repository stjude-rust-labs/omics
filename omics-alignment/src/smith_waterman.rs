//! Smith-Waterman local alignment with affine gap penalties.
//!
//! Computes the optimal local alignment between two sequences using
//! the Smith-Waterman algorithm with affine gap penalties. The local
//! alignment identifies the highest-scoring subsequence pair, ignoring
//! poorly-matching flanking regions.
//!
//! # Examples
//!
//! ```
//! use omics_alignment::score::GapPenalty;
//! use omics_alignment::score::NucleotideMatrix;
//! use omics_alignment::sequence::EncodedSequence;
//! use omics_alignment::smith_waterman;
//! use omics_molecule::polymer::dna;
//!
//! let query = "AAACGTAAAA".parse::<dna::Molecule>()?;
//! let reference = "TTCGTTTT".parse::<dna::Molecule>()?;
//!
//! let result = smith_waterman::align(
//!     &EncodedSequence::from(&query),
//!     &EncodedSequence::from(&reference),
//!     &NucleotideMatrix::new(2, -1),
//!     &GapPenalty::affine(-3, -1),
//! );
//!
//! assert!(result.score() > 0);
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::cigar::Alignment;
use crate::cigar::Cigar;
use crate::cigar::CigarOp;
use crate::score::GapPenalty;
use crate::score::ScoringMatrix;
use crate::sequence::EncodedSequence;

/// Sentinel value representing negative infinity for DP initialization.
const NEG_INF: i32 = i32::MIN / 2;

/// Traceback direction in the H matrix.
#[derive(Clone, Copy, PartialEq)]
enum TraceH {
    /// Cell is zero (local alignment boundary).
    End,
    /// Came from diagonal (match/mismatch).
    Diagonal,
    /// Came from E matrix (gap in query = deletion).
    FromE,
    /// Came from F matrix (gap in reference = insertion).
    FromF,
}

/// Traceback direction in the E or F gap matrices.
#[derive(Clone, Copy)]
enum TraceGap {
    /// Gap was opened from the H matrix.
    Open,
    /// Gap was extended from the same gap matrix.
    Extend,
}

/// Computes the optimal local alignment using Smith-Waterman with affine gap
/// penalties.
///
/// Returns an [`Alignment`] covering only the highest-scoring local region.
pub fn align<M: ScoringMatrix>(
    query: &EncodedSequence,
    reference: &EncodedSequence,
    matrix: &M,
    gap_penalty: &GapPenalty,
) -> Alignment {
    let m = query.len();
    let n = reference.len();

    if m == 0 || n == 0 {
        return Alignment::new(0, Cigar::from(vec![]), 0, 0, 0, 0);
    }

    let rows = m + 1;
    let cols = n + 1;

    let mut h = vec![vec![0i32; cols]; rows];
    let mut e = vec![vec![NEG_INF; cols]; rows];
    let mut f = vec![vec![NEG_INF; cols]; rows];

    let mut trace_h = vec![vec![TraceH::End; cols]; rows];
    let mut trace_e = vec![vec![TraceGap::Open; cols]; rows];
    let mut trace_f = vec![vec![TraceGap::Open; cols]; rows];

    let mut max_score = 0i32;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    // Fill the DP matrices (first row and column remain zero).
    for i in 1..rows {
        for j in 1..cols {
            // E: horizontal gap (gap in query = deletion).
            let e_extend = e[i][j - 1] + gap_penalty.extend();
            let e_open = h[i][j - 1] + gap_penalty.open();
            if e_extend >= e_open {
                e[i][j] = e_extend;
                trace_e[i][j] = TraceGap::Extend;
            } else {
                e[i][j] = e_open;
                trace_e[i][j] = TraceGap::Open;
            }

            // F: vertical gap (gap in reference = insertion).
            let f_extend = f[i - 1][j] + gap_penalty.extend();
            let f_open = h[i - 1][j] + gap_penalty.open();
            if f_extend >= f_open {
                f[i][j] = f_extend;
                trace_f[i][j] = TraceGap::Extend;
            } else {
                f[i][j] = f_open;
                trace_f[i][j] = TraceGap::Open;
            }

            // H: best of zero, diagonal, E, or F.
            let diag = h[i - 1][j - 1] + matrix.score(query.get(i - 1), reference.get(j - 1));
            let mut best = 0;
            let mut trace = TraceH::End;

            if diag > best {
                best = diag;
                trace = TraceH::Diagonal;
            }
            if e[i][j] > best {
                best = e[i][j];
                trace = TraceH::FromE;
            }
            if f[i][j] > best {
                best = f[i][j];
                trace = TraceH::FromF;
            }

            h[i][j] = best;
            trace_h[i][j] = trace;

            if best > max_score {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }
    }

    if max_score == 0 {
        return Alignment::new(0, Cigar::from(vec![]), 0, 0, 0, 0);
    }

    let (cigar, end_i, end_j) =
        traceback(query, reference, max_i, max_j, &trace_h, &trace_e, &trace_f);

    Alignment::new(max_score, cigar, end_i, max_i, end_j, max_j)
}

/// Performs traceback from `(start_i, start_j)` backward until a zero cell,
/// returning the CIGAR and the start coordinates.
fn traceback(
    query: &EncodedSequence,
    reference: &EncodedSequence,
    start_i: usize,
    start_j: usize,
    trace_h: &[Vec<TraceH>],
    trace_e: &[Vec<TraceGap>],
    trace_f: &[Vec<TraceGap>],
) -> (Cigar, usize, usize) {
    /// Which DP matrix the traceback is currently following.
    enum State {
        /// In the H matrix.
        H,
        /// In the E matrix (horizontal gap).
        E,
        /// In the F matrix (vertical gap).
        F,
    }

    let mut ops = Vec::new();
    let mut i = start_i;
    let mut j = start_j;
    let mut state = State::H;

    loop {
        match state {
            State::H => match trace_h[i][j] {
                TraceH::End => break,
                TraceH::Diagonal => {
                    if query.get(i - 1) == reference.get(j - 1) {
                        ops.push(CigarOp::SeqMatch(1));
                    } else {
                        ops.push(CigarOp::SeqMismatch(1));
                    }
                    i -= 1;
                    j -= 1;
                }
                TraceH::FromE => state = State::E,
                TraceH::FromF => state = State::F,
            },
            State::E => {
                ops.push(CigarOp::Deletion(1));
                match trace_e[i][j] {
                    TraceGap::Extend => j -= 1,
                    TraceGap::Open => {
                        j -= 1;
                        state = State::H;
                    }
                }
            }
            State::F => {
                ops.push(CigarOp::Insertion(1));
                match trace_f[i][j] {
                    TraceGap::Extend => i -= 1,
                    TraceGap::Open => {
                        i -= 1;
                        state = State::H;
                    }
                }
            }
        }
    }

    ops.reverse();
    (Cigar::from_ops_merged(ops), i, j)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::score::Blosum62;
    use crate::score::NucleotideMatrix;

    /// Helper to locally align two DNA strings.
    fn align_dna(q: &str, r: &str, m: i32, mm: i32, go: i32, ge: i32) -> Alignment {
        let qm = q.parse::<omics_molecule::polymer::dna::Molecule>().unwrap();
        let rm = r.parse::<omics_molecule::polymer::dna::Molecule>().unwrap();
        align(
            &EncodedSequence::from(&qm),
            &EncodedSequence::from(&rm),
            &NucleotideMatrix::new(m, mm),
            &GapPenalty::affine(go, ge),
        )
    }

    #[test]
    fn it_aligns_identical_sequences() {
        let result = align_dna("ACGT", "ACGT", 2, -1, -3, -1);
        assert_eq!(result.score(), 8);
        assert_eq!(result.cigar().to_string(), "4=");
    }

    #[test]
    fn it_finds_local_alignment_in_flanking_noise() {
        let result = align_dna("AAACGTAAA", "TTTCGTTT", 2, -1, -3, -1);
        assert!(result.score() >= 6); // at least the "CGT" match
        assert!(result.cigar().to_string().contains('='));
    }

    #[test]
    fn it_returns_zero_for_completely_dissimilar_sequences() {
        // With match=1, mismatch=-2, these are always negative when paired.
        let result = align_dna("AAAA", "CCCC", 1, -2, -3, -1);
        assert_eq!(result.score(), 0);
    }

    #[test]
    fn it_handles_empty_sequences() {
        let result = align_dna("", "ACGT", 1, -1, -2, -1);
        assert_eq!(result.score(), 0);

        let result = align_dna("ACGT", "", 1, -1, -2, -1);
        assert_eq!(result.score(), 0);
    }

    #[test]
    fn it_aligns_proteins_with_blosum62() {
        let q = "MKWVTFISLLFLFSSAYS"
            .parse::<omics_molecule::polymer::protein::Molecule>()
            .unwrap();
        let r = "MKWVTFISLLFSSAYS"
            .parse::<omics_molecule::polymer::protein::Molecule>()
            .unwrap();
        let result = align(
            &EncodedSequence::from(&q),
            &EncodedSequence::from(&r),
            &Blosum62,
            &GapPenalty::affine(-11, -1),
        );
        assert!(result.score() > 0);
        assert!(!result.cigar().ops().is_empty());
    }

    #[test]
    fn it_reports_correct_alignment_coordinates() {
        let result = align_dna("AAACGTAAA", "TTTCGTTTT", 2, -1, -3, -1);
        assert!(result.query_start() <= result.query_end());
        assert!(result.reference_start() <= result.reference_end());
        assert!(result.query_end() <= 9);
        assert!(result.reference_end() <= 9);
    }
}
