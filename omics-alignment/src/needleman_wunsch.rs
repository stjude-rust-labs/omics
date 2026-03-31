//! Needleman-Wunsch global alignment with affine gap penalties.
//!
//! Computes the optimal global alignment between two sequences using
//! the Gotoh variant of the Needleman-Wunsch algorithm with separate
//! gap-open and gap-extend costs.
//!
//! # Examples
//!
//! ```
//! use omics_alignment::needleman_wunsch;
//! use omics_alignment::score::GapPenalty;
//! use omics_alignment::score::NucleotideMatrix;
//! use omics_alignment::sequence::EncodedSequence;
//! use omics_molecule::polymer::dna;
//!
//! let query = "ACGTACGT".parse::<dna::Molecule>()?;
//! let reference = "ACGTACGT".parse::<dna::Molecule>()?;
//!
//! let result = needleman_wunsch::align(
//!     &EncodedSequence::from(&query),
//!     &EncodedSequence::from(&reference),
//!     &NucleotideMatrix::new(1, -1),
//!     &GapPenalty::affine(-2, -1),
//! );
//!
//! assert_eq!(result.score(), 8);
//! assert_eq!(result.cigar().to_string(), "8=");
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
#[derive(Clone, Copy)]
enum TraceH {
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

/// Computes the optimal global alignment using Needleman-Wunsch with affine
/// gap penalties.
///
/// Returns an [`Alignment`] spanning the full length of both sequences.
pub fn align<M: ScoringMatrix>(
    query: &EncodedSequence,
    reference: &EncodedSequence,
    matrix: &M,
    gap_penalty: &GapPenalty,
) -> Alignment {
    let m = query.len();
    let n = reference.len();

    if m == 0 && n == 0 {
        return Alignment::new(0, Cigar::from(vec![]), 0, 0, 0, 0);
    }
    if m == 0 {
        return Alignment::new(
            gap_penalty.cost(n),
            Cigar::from(vec![CigarOp::Deletion(n)]),
            0,
            0,
            0,
            n,
        );
    }
    if n == 0 {
        return Alignment::new(
            gap_penalty.cost(m),
            Cigar::from(vec![CigarOp::Insertion(m)]),
            0,
            m,
            0,
            0,
        );
    }

    let rows = m + 1;
    let cols = n + 1;

    let mut h = vec![vec![0i32; cols]; rows];
    let mut e = vec![vec![NEG_INF; cols]; rows];
    let mut f = vec![vec![NEG_INF; cols]; rows];

    let mut trace_h = vec![vec![TraceH::Diagonal; cols]; rows];
    let mut trace_e = vec![vec![TraceGap::Open; cols]; rows];
    let mut trace_f = vec![vec![TraceGap::Open; cols]; rows];

    // Initialize first column (gaps in reference).
    for i in 1..rows {
        h[i][0] = gap_penalty.cost(i);
        f[i][0] = h[i][0];
        trace_h[i][0] = TraceH::FromF;
    }

    // Initialize first row (gaps in query).
    for j in 1..cols {
        h[0][j] = gap_penalty.cost(j);
        e[0][j] = h[0][j];
        trace_h[0][j] = TraceH::FromE;
    }

    // Fill the DP matrices.
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

            // H: best of diagonal, E, or F.
            let diag = h[i - 1][j - 1] + matrix.score(query.get(i - 1), reference.get(j - 1));
            h[i][j] = diag;
            trace_h[i][j] = TraceH::Diagonal;

            if e[i][j] > h[i][j] {
                h[i][j] = e[i][j];
                trace_h[i][j] = TraceH::FromE;
            }
            if f[i][j] > h[i][j] {
                h[i][j] = f[i][j];
                trace_h[i][j] = TraceH::FromF;
            }
        }
    }

    let score = h[m][n];
    let cigar = traceback(query, reference, m, n, &trace_h, &trace_e, &trace_f);

    Alignment::new(score, cigar, 0, m, 0, n)
}

/// Performs traceback from `(i, j)` to `(0, 0)` to reconstruct the CIGAR.
fn traceback(
    query: &EncodedSequence,
    reference: &EncodedSequence,
    start_i: usize,
    start_j: usize,
    trace_h: &[Vec<TraceH>],
    trace_e: &[Vec<TraceGap>],
    trace_f: &[Vec<TraceGap>],
) -> Cigar {
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

    while i > 0 || j > 0 {
        match state {
            State::H => {
                if i == 0 {
                    ops.push(CigarOp::Deletion(j));
                    break;
                }
                if j == 0 {
                    ops.push(CigarOp::Insertion(i));
                    break;
                }
                match trace_h[i][j] {
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
                }
            }
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
    Cigar::from_ops_merged(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::score::NucleotideMatrix;

    /// Helper to align two DNA strings.
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
        let result = align_dna("ACGT", "ACGT", 1, -1, -2, -1);
        assert_eq!(result.score(), 4);
        assert_eq!(result.cigar().to_string(), "4=");
    }

    #[test]
    fn it_aligns_sequences_with_a_mismatch() {
        let result = align_dna("ACGT", "ACGA", 1, -1, -2, -1);
        assert_eq!(result.score(), 2);
        assert_eq!(result.cigar().to_string(), "3=1X");
    }

    #[test]
    fn it_aligns_sequences_with_a_deletion() {
        let result = align_dna("ACGT", "ACGGT", 2, -1, -3, -1);
        let cigar = result.cigar().to_string();
        assert!(cigar.contains('D') || cigar.contains('I'));
    }

    #[test]
    fn it_aligns_sequences_with_an_insertion() {
        let result = align_dna("ACGGT", "ACGT", 2, -1, -3, -1);
        let cigar = result.cigar().to_string();
        assert!(cigar.contains('I') || cigar.contains('D'));
    }

    #[test]
    fn it_handles_empty_query() {
        let result = align_dna("", "ACGT", 1, -1, -2, -1);
        assert_eq!(result.cigar().to_string(), "4D");
        assert_eq!(result.score(), -5);
    }

    #[test]
    fn it_handles_empty_reference() {
        let result = align_dna("ACGT", "", 1, -1, -2, -1);
        assert_eq!(result.cigar().to_string(), "4I");
        assert_eq!(result.score(), -5);
    }

    #[test]
    fn it_spans_full_sequences() {
        let result = align_dna("ACGT", "ACGT", 1, -1, -2, -1);
        assert_eq!(result.query_start(), 0);
        assert_eq!(result.query_end(), 4);
        assert_eq!(result.reference_start(), 0);
        assert_eq!(result.reference_end(), 4);
    }
}
