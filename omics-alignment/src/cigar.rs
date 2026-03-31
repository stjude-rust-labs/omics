//! CIGAR string representation and alignment results.
//!
//! Provides [`CigarOp`] for individual alignment operations, [`Cigar`] as a
//! sequence of operations, and [`Alignment`] as the complete result of a
//! pairwise alignment including score, CIGAR, and coordinate ranges.

use std::fmt;
use std::str::FromStr;

use thiserror::Error;

/// An error related to CIGAR parsing.
#[derive(Error, Debug)]
pub enum Error {
    /// An invalid CIGAR operation character was encountered.
    #[error("invalid CIGAR operation `{0}`")]
    InvalidOperation(char),

    /// An invalid CIGAR string format was encountered.
    #[error("invalid CIGAR format `{0}`")]
    InvalidFormat(String),
}

/// A single operation in a CIGAR string.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum CigarOp {
    /// Sequence match (=): residues are identical.
    SeqMatch(usize),

    /// Sequence mismatch (X): residues differ.
    SeqMismatch(usize),

    /// Insertion to reference (I): present in query but not reference.
    Insertion(usize),

    /// Deletion from reference (D): present in reference but not query.
    Deletion(usize),
}

impl CigarOp {
    /// Returns the length of this operation.
    pub fn len(&self) -> usize {
        match self {
            CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::Deletion(n) => *n,
        }
    }

    /// Returns `true` if this operation has zero length.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the SAM-style character code for this operation.
    fn code(&self) -> char {
        match self {
            CigarOp::SeqMatch(_) => '=',
            CigarOp::SeqMismatch(_) => 'X',
            CigarOp::Insertion(_) => 'I',
            CigarOp::Deletion(_) => 'D',
        }
    }
}

/// A CIGAR string representing the alignment path between two sequences.
///
/// # Examples
///
/// ```
/// use omics_alignment::cigar::Cigar;
/// use omics_alignment::cigar::CigarOp;
///
/// let cigar = Cigar::from(vec![
///     CigarOp::SeqMatch(3),
///     CigarOp::Insertion(1),
///     CigarOp::SeqMatch(4),
/// ]);
/// assert_eq!(cigar.to_string(), "3=1I4=");
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Cigar(Vec<CigarOp>);

impl Cigar {
    /// Returns the operations in this CIGAR string.
    pub fn ops(&self) -> &[CigarOp] {
        &self.0
    }

    /// Constructs a CIGAR from a list of single-residue operations, merging
    /// consecutive operations of the same type.
    pub(crate) fn from_ops_merged(ops: Vec<CigarOp>) -> Self {
        let mut merged: Vec<CigarOp> = Vec::new();
        for op in ops {
            if op.is_empty() {
                continue;
            }
            if let Some(last) = merged.last_mut() {
                let merged_op = match (*last, op) {
                    (CigarOp::SeqMatch(a), CigarOp::SeqMatch(b)) => Some(CigarOp::SeqMatch(a + b)),
                    (CigarOp::SeqMismatch(a), CigarOp::SeqMismatch(b)) => {
                        Some(CigarOp::SeqMismatch(a + b))
                    }
                    (CigarOp::Insertion(a), CigarOp::Insertion(b)) => {
                        Some(CigarOp::Insertion(a + b))
                    }
                    (CigarOp::Deletion(a), CigarOp::Deletion(b)) => Some(CigarOp::Deletion(a + b)),
                    _ => None,
                };
                if let Some(new_op) = merged_op {
                    *last = new_op;
                    continue;
                }
            }
            merged.push(op);
        }
        Self(merged)
    }
}

impl From<Vec<CigarOp>> for Cigar {
    fn from(ops: Vec<CigarOp>) -> Self {
        Self(ops)
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for op in &self.0 {
            write!(f, "{}{}", op.len(), op.code())?;
        }
        Ok(())
    }
}

impl FromStr for Cigar {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ops = Vec::new();
        let mut num_start = 0;

        for (i, c) in s.char_indices() {
            if c.is_ascii_digit() {
                continue;
            }
            let len_str = &s[num_start..i];
            let len: usize = len_str
                .parse()
                .map_err(|_| Error::InvalidFormat(s.to_string()))?;
            let op = match c {
                '=' => CigarOp::SeqMatch(len),
                'X' => CigarOp::SeqMismatch(len),
                'I' => CigarOp::Insertion(len),
                'D' => CigarOp::Deletion(len),
                _ => return Err(Error::InvalidOperation(c)),
            };
            ops.push(op);
            num_start = i + 1;
        }

        if num_start != s.len() {
            return Err(Error::InvalidFormat(s.to_string()));
        }

        Ok(Self(ops))
    }
}

/// The result of a pairwise sequence alignment.
///
/// Contains the alignment score, CIGAR string, and the coordinate ranges in
/// both the query and reference sequences.
///
/// # Examples
///
/// ```
/// use omics_alignment::cigar::Alignment;
/// use omics_alignment::cigar::Cigar;
/// use omics_alignment::cigar::CigarOp;
///
/// let alignment = Alignment::new(10, Cigar::from(vec![CigarOp::SeqMatch(5)]), 0, 5, 0, 5);
/// assert_eq!(alignment.score(), 10);
/// assert_eq!(alignment.cigar().to_string(), "5=");
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Alignment {
    /// The alignment score.
    score: i32,

    /// The CIGAR string.
    cigar: Cigar,

    /// The start position in the query (0-based, inclusive).
    query_start: usize,

    /// The end position in the query (0-based, exclusive).
    query_end: usize,

    /// The start position in the reference (0-based, inclusive).
    reference_start: usize,

    /// The end position in the reference (0-based, exclusive).
    reference_end: usize,
}

impl Alignment {
    /// Creates a new alignment result.
    pub fn new(
        score: i32,
        cigar: Cigar,
        query_start: usize,
        query_end: usize,
        reference_start: usize,
        reference_end: usize,
    ) -> Self {
        Self {
            score,
            cigar,
            query_start,
            query_end,
            reference_start,
            reference_end,
        }
    }

    /// Returns the alignment score.
    pub fn score(&self) -> i32 {
        self.score
    }

    /// Returns the CIGAR string.
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns the start position in the query (0-based, inclusive).
    pub fn query_start(&self) -> usize {
        self.query_start
    }

    /// Returns the end position in the query (0-based, exclusive).
    pub fn query_end(&self) -> usize {
        self.query_end
    }

    /// Returns the start position in the reference (0-based, inclusive).
    pub fn reference_start(&self) -> usize {
        self.reference_start
    }

    /// Returns the end position in the reference (0-based, exclusive).
    pub fn reference_end(&self) -> usize {
        self.reference_end
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_displays_cigar_strings() {
        let cigar = Cigar::from(vec![
            CigarOp::SeqMatch(3),
            CigarOp::Insertion(2),
            CigarOp::SeqMatch(4),
            CigarOp::Deletion(1),
        ]);
        assert_eq!(cigar.to_string(), "3=2I4=1D");
    }

    #[test]
    fn it_parses_cigar_strings() -> Result<(), Box<dyn std::error::Error>> {
        let cigar = "3=2I4=1D".parse::<Cigar>()?;
        assert_eq!(
            cigar.ops(),
            &[
                CigarOp::SeqMatch(3),
                CigarOp::Insertion(2),
                CigarOp::SeqMatch(4),
                CigarOp::Deletion(1),
            ]
        );
        Ok(())
    }

    #[test]
    fn it_round_trips_cigar_strings() -> Result<(), Box<dyn std::error::Error>> {
        let original = "5=3X2I1D4=";
        let cigar = original.parse::<Cigar>()?;
        assert_eq!(cigar.to_string(), original);
        Ok(())
    }

    #[test]
    fn it_merges_consecutive_operations() {
        let ops = vec![
            CigarOp::SeqMatch(1),
            CigarOp::SeqMatch(1),
            CigarOp::SeqMatch(1),
            CigarOp::Insertion(1),
            CigarOp::Insertion(1),
            CigarOp::SeqMismatch(1),
        ];
        let cigar = Cigar::from_ops_merged(ops);
        assert_eq!(cigar.to_string(), "3=2I1X");
    }

    #[test]
    fn it_rejects_invalid_cigar_operations() {
        let err = "3Z".parse::<Cigar>().unwrap_err();
        assert!(matches!(err, Error::InvalidOperation('Z')));
    }
}
