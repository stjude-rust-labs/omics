//! Scoring matrices and gap penalty models.
//!
//! Provides the [`ScoringMatrix`] trait for abstracting over substitution
//! matrices, concrete implementations for nucleotide and protein alignment, and
//! the [`GapPenalty`] type for affine gap penalty models.

/// A substitution scoring matrix for pairwise alignment.
///
/// Implementations map pairs of encoded residue indices to integer scores.
/// Residue indices are assigned by
/// [`EncodedSequence`](crate::sequence::EncodedSequence) during encoding (e.g.,
/// A=0, C=1, G=2, T=3 for DNA).
pub trait ScoringMatrix {
    /// Returns the substitution score for aligning residue `a` against residue
    /// `b`, where both are encoded residue indices.
    fn score(&self, a: u8, b: u8) -> i32;
}

/// An affine gap penalty model.
///
/// A gap of length *n* costs `open + (n - 1) * extend`, where
/// [`open`](Self::open) is the cost of the first gap position and
/// [`extend`](Self::extend) is the cost of each subsequent position. Both
/// values should be negative (or zero).
///
/// # Examples
///
/// ```
/// use omics_alignment::score::GapPenalty;
///
/// // Affine gap: opening costs -11, extension costs -1.
/// let affine = GapPenalty::affine(-11, -1);
/// assert_eq!(affine.cost(1), -11);
/// assert_eq!(affine.cost(3), -13);
///
/// // Linear gap: every position costs -2.
/// let linear = GapPenalty::linear(-2);
/// assert_eq!(linear.cost(1), -2);
/// assert_eq!(linear.cost(3), -6);
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct GapPenalty {
    /// The cost of opening a gap (first position).
    open: i32,

    /// The cost of extending an existing gap (each additional position).
    extend: i32,
}

impl GapPenalty {
    /// Creates an affine gap penalty with the given open and extend costs.
    pub fn affine(open: i32, extend: i32) -> Self {
        Self { open, extend }
    }

    /// Creates a linear gap penalty where every gap position has the same cost.
    pub fn linear(penalty: i32) -> Self {
        Self {
            open: penalty,
            extend: penalty,
        }
    }

    /// Returns the open cost.
    pub fn open(&self) -> i32 {
        self.open
    }

    /// Returns the extend cost.
    pub fn extend(&self) -> i32 {
        self.extend
    }

    /// Returns the total cost of a gap of length `n`.
    pub fn cost(&self, n: usize) -> i32 {
        if n == 0 {
            return 0;
        }
        self.open + (n as i32 - 1) * self.extend
    }
}

/// A simple nucleotide scoring matrix using match/mismatch scores.
///
/// # Examples
///
/// ```
/// use omics_alignment::score::NucleotideMatrix;
/// use omics_alignment::score::ScoringMatrix;
///
/// let matrix = NucleotideMatrix::new(1, -1);
/// assert_eq!(matrix.score(0, 0), 1); // match
/// assert_eq!(matrix.score(0, 1), -1); // mismatch
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct NucleotideMatrix {
    /// Score awarded for matching residues.
    match_score: i32,

    /// Score penalized for mismatching residues.
    mismatch_score: i32,
}

impl NucleotideMatrix {
    /// Creates a new nucleotide scoring matrix.
    pub fn new(match_score: i32, mismatch_score: i32) -> Self {
        Self {
            match_score,
            mismatch_score,
        }
    }
}

impl ScoringMatrix for NucleotideMatrix {
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The BLOSUM62 substitution matrix for protein sequence alignment.
///
/// Residues are indexed in standard BLOSUM order:
/// A=0, R=1, N=2, D=3, C=4, Q=5, E=6, G=7, H=8, I=9,
/// L=10, K=11, M=12, F=13, P=14, S=15, T=16, W=17, Y=18, V=19.
///
/// # Examples
///
/// ```
/// use omics_alignment::score::Blosum62;
/// use omics_alignment::score::ScoringMatrix;
///
/// let matrix = Blosum62;
/// assert_eq!(matrix.score(0, 0), 4); // A-A
/// assert_eq!(matrix.score(4, 4), 9); // C-C
/// assert_eq!(matrix.score(17, 17), 11); // W-W
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Blosum62;

impl ScoringMatrix for Blosum62 {
    fn score(&self, a: u8, b: u8) -> i32 {
        BLOSUM62_DATA[a as usize][b as usize]
    }
}

/// Raw BLOSUM62 data indexed as [row][col] in standard amino acid order.
//
// Source: Henikoff & Henikoff (1992), PNAS 89(22):10915-10919.
// Order: A R N D C Q E G H I L K M F P S T W Y V
#[rustfmt::skip]
const BLOSUM62_DATA: [[i32; 20]; 20] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], // A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], // R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], // N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], // D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], // C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], // Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], // E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], // G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], // H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], // I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], // L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], // K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], // M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], // F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], // P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], // S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], // T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], // W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], // Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4], // V
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_calculates_affine_gap_costs() {
        let gap = GapPenalty::affine(-11, -1);
        assert_eq!(gap.cost(0), 0);
        assert_eq!(gap.cost(1), -11);
        assert_eq!(gap.cost(2), -12);
        assert_eq!(gap.cost(5), -15);
    }

    #[test]
    fn it_calculates_linear_gap_costs() {
        let gap = GapPenalty::linear(-2);
        assert_eq!(gap.cost(0), 0);
        assert_eq!(gap.cost(1), -2);
        assert_eq!(gap.cost(3), -6);
    }

    #[test]
    fn it_scores_nucleotide_matches_and_mismatches() {
        let m = NucleotideMatrix::new(2, -1);
        assert_eq!(m.score(0, 0), 2);
        assert_eq!(m.score(0, 1), -1);
        assert_eq!(m.score(3, 3), 2);
    }

    #[test]
    fn it_returns_correct_blosum62_diagonal_values() {
        let m = Blosum62;
        assert_eq!(m.score(0, 0), 4); // A-A
        assert_eq!(m.score(4, 4), 9); // C-C
        assert_eq!(m.score(8, 8), 8); // H-H
        assert_eq!(m.score(17, 17), 11); // W-W
    }

    #[test]
    fn it_returns_symmetric_blosum62_values() {
        let m = Blosum62;
        for a in 0..20u8 {
            for b in 0..20u8 {
                assert_eq!(m.score(a, b), m.score(b, a));
            }
        }
    }
}
