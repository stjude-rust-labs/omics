//! CIGAR strings for representing pairwise sequence alignments.
//!
//! A CIGAR string encodes the operations needed to transform a query sequence
//! into alignment with a reference sequence. Each operation has a character
//! code and a positive length.
//!
//! # Quick start
//!
//! ```
//! use omics_alignment::cigar::Cigar;
//!
//! let cigar: Cigar = "10M2I3D".parse().unwrap();
//! assert_eq!(cigar.reference_length(), 13);
//! assert_eq!(cigar.query_length(), 12);
//! assert_eq!(cigar.to_string(), "10M2I3D");
//! ```

use std::fmt;
use std::num::ParseIntError;
use std::str::FromStr;

use omics_coordinate::position::Number;
use thiserror::Error;

/// Which axis of an alignment a length overflow occurred on.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Axis {
    /// The reference (genomic) axis.
    Reference,
    /// The query (read) axis.
    Query,
}

impl fmt::Display for Axis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Axis::Reference => write!(f, "reference"),
            Axis::Query => write!(f, "query"),
        }
    }
}

/// An error when constructing an [`Operation`].
#[derive(Clone, Copy, Debug, Eq, Error, Hash, PartialEq)]
pub enum OperationError {
    /// The length supplied to [`Operation::try_new`] was zero.
    #[error("operation length must be non-zero")]
    ZeroLength,
}

/// An error when constructing or parsing a [`Cigar`].
#[derive(Debug, Error)]
pub enum Error {
    /// The CIGAR string or operation list was empty.
    #[error("empty cigar")]
    Empty,

    /// The CIGAR field contained the unavailable marker (`*`).
    #[error("cigar is unavailable")]
    Unavailable,

    /// An operation was encountered without a preceding length.
    #[error("operation {operation_index} at offset {offset} has no length")]
    MissingLength {
        /// Zero-based index of the operation in the string.
        operation_index: usize,
        /// Byte offset within the input string where the operation code
        /// appeared.
        offset: usize,
    },

    /// The length text for an operation could not be parsed as a
    /// [`Number`].
    #[error("operation {operation_index} has invalid length `{value}`")]
    InvalidLength {
        /// Zero-based index of the operation in the string.
        operation_index: usize,
        /// The text that failed to parse.
        value: String,
        /// The underlying parse error.
        #[source]
        source: ParseIntError,
    },

    /// An unrecognised operation code was encountered.
    #[error("operation {operation_index} at offset {offset} has unknown kind `{value}`")]
    InvalidKind {
        /// Zero-based index of the operation in the string.
        operation_index: usize,
        /// Byte offset within the input string where the unknown code appeared.
        offset: usize,
        /// The unrecognised character.
        value: char,
    },

    /// An [`Operation`] could not be constructed from a parsed length.
    #[error("operation {operation_index}: {source}")]
    Operation {
        /// Zero-based index of the operation in the string.
        operation_index: usize,
        /// The underlying construction error.
        #[source]
        source: OperationError,
    },

    /// A clipping operation appeared at a position not allowed by the SAM
    /// specification.
    #[error("operation {operation_index} has invalid clipping placement for {kind:?}")]
    InvalidClipping {
        /// Zero-based index of the offending operation.
        operation_index: usize,
        /// The kind of clipping operation that was misplaced.
        kind: OperationKind,
    },

    /// Digits were present at the end of the string with no following
    /// operation code.
    #[error("trailing length `{value}` with no operation code")]
    TrailingLength {
        /// The trailing digit text.
        value: String,
    },

    /// The accumulated length for one axis overflowed [`Number`].
    #[error("{axis} length overflows")]
    LengthOverflow {
        /// The axis on which the overflow occurred.
        axis: Axis,
    },
}

/// A single CIGAR operation kind.
///
/// Each variant corresponds to one of the nine operation codes defined in the
/// SAM specification (M, I, D, N, S, H, P, =, X).
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum OperationKind {
    /// Alignment match (M). Consumes both reference and query.
    Match,
    /// Insertion to the reference (I). Consumes query only.
    Insertion,
    /// Deletion from the reference (D). Consumes reference only.
    Deletion,
    /// Reference skip (N). Consumes reference only.
    ReferenceSkip,
    /// Soft clipping (S). Consumes query only.
    SoftClip,
    /// Hard clipping (H). Consumes neither reference nor query.
    HardClip,
    /// Padding (P). Consumes neither reference nor query.
    Padding,
    /// Sequence match (=). Consumes both reference and query.
    SequenceMatch,
    /// Sequence mismatch (X). Consumes both reference and query.
    SequenceMismatch,
}

impl OperationKind {
    /// Returns `true` if this operation advances the reference position.
    pub const fn consumes_reference(self) -> bool {
        matches!(
            self,
            OperationKind::Match
                | OperationKind::Deletion
                | OperationKind::ReferenceSkip
                | OperationKind::SequenceMatch
                | OperationKind::SequenceMismatch
        )
    }

    /// Returns `true` if this operation advances the query position.
    pub const fn consumes_query(self) -> bool {
        matches!(
            self,
            OperationKind::Match
                | OperationKind::Insertion
                | OperationKind::SoftClip
                | OperationKind::SequenceMatch
                | OperationKind::SequenceMismatch
        )
    }

    /// Returns `true` if this operation contributes to the pairwise base
    /// alignment (M, =, or X).
    pub const fn is_aligned(self) -> bool {
        matches!(
            self,
            OperationKind::Match | OperationKind::SequenceMatch | OperationKind::SequenceMismatch
        )
    }
}

impl fmt::Display for OperationKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            OperationKind::Match => 'M',
            OperationKind::Insertion => 'I',
            OperationKind::Deletion => 'D',
            OperationKind::ReferenceSkip => 'N',
            OperationKind::SoftClip => 'S',
            OperationKind::HardClip => 'H',
            OperationKind::Padding => 'P',
            OperationKind::SequenceMatch => '=',
            OperationKind::SequenceMismatch => 'X',
        };
        write!(f, "{c}")
    }
}

impl TryFrom<char> for OperationKind {
    type Error = ();

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'M' => Ok(OperationKind::Match),
            'I' => Ok(OperationKind::Insertion),
            'D' => Ok(OperationKind::Deletion),
            'N' => Ok(OperationKind::ReferenceSkip),
            'S' => Ok(OperationKind::SoftClip),
            'H' => Ok(OperationKind::HardClip),
            'P' => Ok(OperationKind::Padding),
            '=' => Ok(OperationKind::SequenceMatch),
            'X' => Ok(OperationKind::SequenceMismatch),
            _ => Err(()),
        }
    }
}

/// A single CIGAR operation with a kind and a non-zero length.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Operation {
    /// The operation kind.
    kind: OperationKind,
    /// The operation length. Always non-zero.
    length: Number,
}

impl Operation {
    /// Constructs a new [`Operation`].
    ///
    /// Returns [`OperationError::ZeroLength`] if `length` is zero.
    pub fn try_new(kind: OperationKind, length: Number) -> Result<Self, OperationError> {
        if length == 0 {
            return Err(OperationError::ZeroLength);
        }
        Ok(Self { kind, length })
    }

    /// Returns the operation kind.
    pub const fn kind(&self) -> OperationKind {
        self.kind
    }

    /// Returns the operation length.
    pub const fn length(&self) -> Number {
        self.length
    }
}

impl fmt::Display for Operation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.length, self.kind)
    }
}

/// A validated CIGAR string.
///
/// A [`Cigar`] wraps an ordered sequence of [`Operation`]s and caches the
/// total consumed lengths for the reference and query axes. All clipping
/// placement rules are enforced on construction.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Cigar {
    /// The ordered operations.
    operations: Vec<Operation>,
    /// Total bases consumed on the reference axis.
    reference_length: Number,
    /// Total bases consumed on the query axis.
    query_length: Number,
}

impl Cigar {
    /// Constructs a [`Cigar`] from a vector of operations.
    ///
    /// Returns an error if the vector is empty, if any clipping operation is
    /// misplaced, or if the accumulated reference or query length overflows
    /// [`Number`].
    pub fn try_new(operations: Vec<Operation>) -> Result<Self, Error> {
        if operations.is_empty() {
            return Err(Error::Empty);
        }

        for (index, op) in operations.iter().enumerate() {
            match op.kind() {
                OperationKind::HardClip => {
                    let terminal = index == 0 || index + 1 == operations.len();
                    if !terminal {
                        return Err(Error::InvalidClipping {
                            operation_index: index,
                            kind: OperationKind::HardClip,
                        });
                    }
                }
                OperationKind::SoftClip => {
                    let terminal = index == 0 || index + 1 == operations.len();
                    let soft_terminal = terminal
                        || (index == 1 && operations[0].kind() == OperationKind::HardClip)
                        || (index + 2 == operations.len()
                            && operations[index + 1].kind() == OperationKind::HardClip);
                    if !soft_terminal {
                        return Err(Error::InvalidClipping {
                            operation_index: index,
                            kind: OperationKind::SoftClip,
                        });
                    }
                }
                _ => {}
            }
        }

        let mut reference_length: Number = 0;
        let mut query_length: Number = 0;

        for op in &operations {
            let len = op.length();
            if op.kind().consumes_reference() {
                reference_length =
                    reference_length
                        .checked_add(len)
                        .ok_or(Error::LengthOverflow {
                            axis: Axis::Reference,
                        })?;
            }
            if op.kind().consumes_query() {
                query_length = query_length
                    .checked_add(len)
                    .ok_or(Error::LengthOverflow { axis: Axis::Query })?;
            }
        }

        Ok(Self {
            operations,
            reference_length,
            query_length,
        })
    }

    /// Returns the operations as a slice.
    pub fn operations(&self) -> &[Operation] {
        &self.operations
    }

    /// Returns an iterator over the operations.
    pub fn iter(&self) -> std::slice::Iter<'_, Operation> {
        self.operations.iter()
    }

    /// Returns the total number of reference bases consumed by this CIGAR.
    pub const fn reference_length(&self) -> Number {
        self.reference_length
    }

    /// Returns the total number of query bases consumed by this CIGAR.
    pub const fn query_length(&self) -> Number {
        self.query_length
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for op in &self.operations {
            write!(f, "{op}")?;
        }
        Ok(())
    }
}

impl FromStr for Cigar {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(Error::Empty);
        }
        if s == "*" {
            return Err(Error::Unavailable);
        }

        let mut operations: Vec<Operation> = Vec::new();
        let mut length_text = String::new();

        for (offset, ch) in s.char_indices() {
            if ch.is_ascii_digit() {
                length_text.push(ch);
            } else {
                let operation_index = operations.len();

                if length_text.is_empty() {
                    return Err(Error::MissingLength {
                        operation_index,
                        offset,
                    });
                }

                let length: Number =
                    length_text.parse().map_err(|source| Error::InvalidLength {
                        operation_index,
                        value: length_text.clone(),
                        source,
                    })?;
                length_text.clear();

                let kind = OperationKind::try_from(ch).map_err(|()| Error::InvalidKind {
                    operation_index,
                    offset,
                    value: ch,
                })?;

                let op = Operation::try_new(kind, length).map_err(|source| Error::Operation {
                    operation_index,
                    source,
                })?;

                operations.push(op);
            }
        }

        if !length_text.is_empty() {
            return Err(Error::TrailingLength { value: length_text });
        }

        Cigar::try_new(operations)
    }
}

impl<'a> IntoIterator for &'a Cigar {
    type IntoIter = std::slice::Iter<'a, Operation>;
    type Item = &'a Operation;

    fn into_iter(self) -> Self::IntoIter {
        self.operations.iter()
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::position::Number;

    use super::*;

    #[test]
    fn consumption_table() {
        let cases = [
            (OperationKind::Match, true, true, true),
            (OperationKind::Insertion, false, true, false),
            (OperationKind::Deletion, true, false, false),
            (OperationKind::ReferenceSkip, true, false, false),
            (OperationKind::SoftClip, false, true, false),
            (OperationKind::HardClip, false, false, false),
            (OperationKind::Padding, false, false, false),
            (OperationKind::SequenceMatch, true, true, true),
            (OperationKind::SequenceMismatch, true, true, true),
        ];

        for (kind, ref_, query, aligned) in cases {
            assert_eq!(kind.consumes_reference(), ref_, "{kind:?}");
            assert_eq!(kind.consumes_query(), query, "{kind:?}");
            assert_eq!(kind.is_aligned(), aligned, "{kind:?}");
        }
    }

    #[test]
    fn zero_length_rejected() {
        assert!(matches!(
            Operation::try_new(OperationKind::Match, 0),
            Err(OperationError::ZeroLength)
        ));
    }

    #[test]
    fn roundtrip() -> Result<(), Box<dyn std::error::Error>> {
        // A valid 9-operation CIGAR with H and S placed at terminal positions.
        let s = "6H5S10M2I3D4N7P8=9X5S6H";
        assert_eq!(s.parse::<Cigar>()?.to_string(), s);
        Ok(())
    }

    #[test]
    fn adjacent_preserved() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!("1M2M".parse::<Cigar>()?.operations().len(), 2);
        Ok(())
    }

    #[test]
    fn invalid_clipping_hard_in_middle() {
        assert!(matches!(
            "1M1H1M".parse::<Cigar>(),
            Err(Error::InvalidClipping { .. })
        ));
    }

    #[test]
    fn invalid_clipping_soft_in_middle() {
        assert!(matches!(
            "1M1S1M".parse::<Cigar>(),
            Err(Error::InvalidClipping { .. })
        ));
    }

    #[test]
    fn valid_clipping_hard_at_ends() -> Result<(), Box<dyn std::error::Error>> {
        let c = "1H10M1H".parse::<Cigar>()?;
        assert_eq!(c.operations().len(), 3);
        Ok(())
    }

    #[test]
    fn valid_clipping_soft_adjacent_to_hard() -> Result<(), Box<dyn std::error::Error>> {
        let c = "1H1S10M1S1H".parse::<Cigar>()?;
        assert_eq!(c.operations().len(), 5);
        Ok(())
    }

    #[test]
    fn empty_string() {
        assert!(matches!("".parse::<Cigar>(), Err(Error::Empty)));
    }

    #[test]
    fn unavailable_marker() {
        assert!(matches!("*".parse::<Cigar>(), Err(Error::Unavailable)));
    }

    #[test]
    fn missing_length() {
        assert!(matches!(
            "M".parse::<Cigar>(),
            Err(Error::MissingLength { .. })
        ));
    }

    #[test]
    fn invalid_kind() {
        assert!(matches!(
            "1Z".parse::<Cigar>(),
            Err(Error::InvalidKind { .. })
        ));
    }

    #[test]
    fn trailing_length() {
        assert!(matches!(
            "1M2".parse::<Cigar>(),
            Err(Error::TrailingLength { .. })
        ));
    }

    #[test]
    fn per_token_overflow() {
        let mut big = Number::MAX.to_string();
        big.push('0');
        let cigar_str = format!("{big}M");
        assert!(matches!(
            cigar_str.parse::<Cigar>(),
            Err(Error::InvalidLength { .. })
        ));
    }

    #[test]
    fn aggregate_reference_overflow() {
        let cigar_str = format!("{}D1D", Number::MAX);
        assert!(matches!(
            cigar_str.parse::<Cigar>(),
            Err(Error::LengthOverflow { .. })
        ));
    }

    #[test]
    fn aggregate_query_overflow() {
        let cigar_str = format!("{}I1I", Number::MAX);
        assert!(matches!(
            cigar_str.parse::<Cigar>(),
            Err(Error::LengthOverflow { .. })
        ));
    }

    #[test]
    fn reference_length() -> Result<(), Box<dyn std::error::Error>> {
        let c = "3M2I4D".parse::<Cigar>()?;
        assert_eq!(c.reference_length(), 7);
        Ok(())
    }

    #[test]
    fn query_length() -> Result<(), Box<dyn std::error::Error>> {
        let c = "3M2I4D".parse::<Cigar>()?;
        assert_eq!(c.query_length(), 5);
        Ok(())
    }

    #[test]
    fn iter_count() -> Result<(), Box<dyn std::error::Error>> {
        let c = "3M2I4D".parse::<Cigar>()?;
        assert_eq!(c.iter().count(), 3);
        Ok(())
    }

    #[test]
    fn into_iter_for_ref() -> Result<(), Box<dyn std::error::Error>> {
        let c = "3M2I4D".parse::<Cigar>()?;
        assert_eq!((&c).into_iter().count(), 3);
        Ok(())
    }
}
