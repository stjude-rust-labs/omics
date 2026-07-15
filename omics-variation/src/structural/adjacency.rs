//! Novel adjacencies between breakends.

use std::fmt;
use std::str::FromStr;

use omics_molecule::compound::Complement;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence;
use omics_molecule::sequence::Sequence;
use thiserror::Error;

use crate::structural::breakend;
use crate::structural::breakend::Breakend;

/// The separator between an adjacency's breakend and insertion fields.
///
/// It is distinct from the breakend-internal separator so a whole breakend
/// token never contains it, keeping the split unambiguous.
const ADJACENCY_SEPARATOR: &str = "::";

/// The token marking the open side of a single-ended adjacency.
///
/// It sits in the middle field. A paired breakend token is never this value,
/// so the middle field alone disambiguates single-ended from paired.
const OPEN_SIDE: &str = ".";

/// An error related to constructing an [`Adjacency`].
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// The two breakends sit at the same locus with opposite orientations and
    /// an empty insertion, which encodes no change.
    #[error("adjacency is an identity join and encodes no change")]
    IdentityJoin,

    /// The two breakends are fully identical, meaning the same contig,
    /// position, and orientation, so the adjacency joins a breakend to itself.
    #[error("adjacency cannot join a breakend to itself")]
    SelfJunction,
}

/// A parse error related to an [`Adjacency`].
#[derive(Error, Debug)]
pub enum ParseError {
    /// The adjacency did not split into exactly three `::`-separated fields.
    #[error("invalid adjacency format: `{0}`")]
    Format(String),

    /// A breakend field failed to parse.
    #[error(transparent)]
    Breakend(#[from] breakend::ParseError),

    /// The insertion sequence failed to parse.
    #[error(transparent)]
    Insertion(#[from] sequence::ParseError),

    /// The paired adjacency could not be constructed.
    #[error(transparent)]
    Construct(#[from] Error),
}

/// The private inner representation of an [`Adjacency`].
///
/// Keeping this enum private makes [`Adjacency`] opaque, so callers cannot
/// build a non-canonical or identity adjacency that would break classification.
/// Adjacencies are built only through [`Adjacency::try_new_paired`] and
/// [`Adjacency::new_single`], and are read back through
/// [`Adjacency::paired`] and [`Adjacency::single`].
#[derive(Debug, Clone, PartialEq, Eq)]
enum Inner<N: Nucleotide> {
    /// Two loci joined, with optional non-templated bases between them.
    Paired {
        /// The canonical lower-locus breakend.
        a: Breakend,

        /// The canonical higher-locus breakend.
        b: Breakend,

        /// The non-templated insertion, in the reading frame of `a`.
        insertion: Sequence<N>,
    },

    /// One locus joined to novel or unassembled sequence, open on the other
    /// side.
    Single {
        /// The single known breakend.
        breakend: Breakend,

        /// The non-templated insertion at the open side.
        insertion: Sequence<N>,
    },
}

/// A single novel junction between one or two oriented breakends.
///
/// An adjacency is opaque. A paired adjacency is always held in canonical form,
/// with its breakends ordered by their canonical key and its insertion
/// carried in the reading frame of the canonical lower breakend, so that
/// equality and later classification do not depend on the input order.
///
/// `Hash` is intentionally not derived because `Sequence` does not implement
/// it.
///
/// An adjacency serializes as three `::`-separated fields. A paired adjacency
/// renders both breakends and its insertion, while a single-ended one renders
/// its one breakend, a literal `.` marking the open side, and its insertion.
///
/// # Examples
///
/// ```
/// use omics_molecule::polymer::dna;
/// use omics_variation::structural::adjacency::Adjacency;
///
/// let rendered = "seq0:>:100(i)::seq1:<:200(i)::AT";
/// let adjacency = rendered.parse::<Adjacency<dna::Nucleotide>>()?;
/// assert_eq!(adjacency.to_string(), rendered);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adjacency<N: Nucleotide> {
    /// The private, canonical inner representation.
    inner: Inner<N>,
}

impl<N: Nucleotide + Complement> Adjacency<N> {
    /// Attempts to create a paired [`Adjacency`] in canonical form.
    ///
    /// The two breakends are ordered by their canonical key. When
    /// ordering swaps the supplied pair, the insertion is reverse-complemented
    /// into the reading frame of the canonical lower breakend. The identity
    /// join, meaning two breakends at the same locus with opposite orientations
    /// and an empty insertion, is rejected because it encodes no change.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::try_new_paired(a, b, ".".parse()?)?;
    /// assert!(adjacency.paired().is_some());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new_paired(x: Breakend, y: Breakend, insertion: Sequence<N>) -> Result<Self, Error> {
        let (a, b, insertion) = if x.canonical_key() <= y.canonical_key() {
            (x, y, insertion)
        } else {
            (y, x, insertion.reverse_complement())
        };

        let same_locus = a.contig() == b.contig() && a.position().get() == b.position().get();

        // Reject a breakend joined to itself. Fully identical breakends, meaning
        // the same contig, position, and orientation, encode a meaningless
        // zero-length self-loop that would otherwise misclassify downstream.
        if same_locus && a.orientation() == b.orientation() {
            return Err(Error::SelfJunction);
        }

        let opposite = a.orientation().is_opposite(b.orientation());
        if same_locus && opposite && insertion.is_empty() {
            return Err(Error::IdentityJoin);
        }

        Ok(Self {
            inner: Inner::Paired { a, b, insertion },
        })
    }
}

impl<N: Nucleotide> Adjacency<N> {
    /// Creates a single-ended [`Adjacency`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::new_single(breakend, "AT".parse()?);
    /// assert!(adjacency.single().is_some());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new_single(breakend: Breakend, insertion: Sequence<N>) -> Self {
        Self {
            inner: Inner::Single {
                breakend,
                insertion,
            },
        }
    }

    /// Gets the two canonical breakends and insertion of a paired adjacency, or
    /// `None` for a single-ended one.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::try_new_paired(a, b, ".".parse()?)?;
    /// let (lower, higher, _) = adjacency.paired().expect("a paired adjacency");
    /// assert_eq!(lower.position().get(), 100);
    /// assert_eq!(higher.position().get(), 200);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn paired(&self) -> Option<(&Breakend, &Breakend, &Sequence<N>)> {
        match &self.inner {
            Inner::Paired { a, b, insertion } => Some((a, b, insertion)),
            Inner::Single { .. } => None,
        }
    }

    /// Gets the breakend and insertion of a single-ended adjacency, or `None`
    /// for a paired one.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::new_single(breakend, "AT".parse()?);
    /// let (breakend, insertion) = adjacency.single().expect("a single adjacency");
    /// assert_eq!(breakend.position().get(), 100);
    /// assert_eq!(insertion.to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn single(&self) -> Option<(&Breakend, &Sequence<N>)> {
        match &self.inner {
            Inner::Single {
                breakend,
                insertion,
            } => Some((breakend, insertion)),
            Inner::Paired { .. } => None,
        }
    }
}

impl<N: Nucleotide> fmt::Display for Adjacency<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some((a, b, insertion)) = self.paired() {
            write!(
                f,
                "{a}{ADJACENCY_SEPARATOR}{b}{ADJACENCY_SEPARATOR}{insertion}"
            )
        } else if let Some((breakend, insertion)) = self.single() {
            write!(
                f,
                "{breakend}{ADJACENCY_SEPARATOR}{OPEN_SIDE}{ADJACENCY_SEPARATOR}{insertion}"
            )
        } else {
            unreachable!("an adjacency is always paired or single")
        }
    }
}

impl<N: Nucleotide + Complement> FromStr for Adjacency<N> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(ADJACENCY_SEPARATOR).collect::<Vec<_>>();
        let [first, middle, insertion] = parts.as_slice() else {
            return Err(ParseError::Format(s.to_string()));
        };

        let insertion = insertion.parse::<Sequence<N>>()?;

        if *middle == OPEN_SIDE {
            let breakend = first.parse::<Breakend>()?;
            Ok(Adjacency::new_single(breakend, insertion))
        } else {
            let a = first.parse::<Breakend>()?;
            let b = middle.parse::<Breakend>()?;
            Ok(Adjacency::try_new_paired(a, b, insertion)?)
        }
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::position::Number;
    use omics_molecule::polymer::dna;

    use super::*;
    use crate::structural::orientation::Orientation;

    fn bnd(contig: &str, orientation: Orientation, position: Number) -> Breakend {
        Breakend::try_new(contig, orientation, position).unwrap()
    }

    fn seq(value: &str) -> Sequence<dna::Nucleotide> {
        value.parse().unwrap()
    }

    #[test]
    fn it_canonicalizes_supplied_order() {
        let lo = bnd("seq0", Orientation::LowerFlank, 100);
        let hi = bnd("seq0", Orientation::HigherFlank, 200);

        let forward = Adjacency::try_new_paired(lo.clone(), hi.clone(), seq(".")).unwrap();
        let backward = Adjacency::try_new_paired(hi, lo, seq(".")).unwrap();
        assert_eq!(forward, backward);
    }

    #[test]
    fn it_reverse_complements_the_insertion_on_swap() {
        let lo = bnd("seq0", Orientation::LowerFlank, 100);
        let hi = bnd("seq0", Orientation::HigherFlank, 200);

        // Supplied high-to-low, so the insertion is reverse-complemented into
        // the canonical low-to-high frame. `AAAC` reverse-complements to `GTTT`.
        let adjacency = Adjacency::try_new_paired(hi, lo, seq("AAAC")).unwrap();
        let (_, _, insertion) = adjacency.paired().expect("a paired adjacency");
        assert_eq!(insertion.to_string(), "GTTT");
    }

    #[test]
    fn it_exposes_a_single_adjacency_through_the_accessor() {
        let breakend = bnd("seq0", Orientation::LowerFlank, 100);
        let adjacency = Adjacency::new_single(breakend, seq("AT"));
        assert!(adjacency.paired().is_none());
        let (_, insertion) = adjacency.single().expect("a single adjacency");
        assert_eq!(insertion.to_string(), "AT");
    }

    #[test]
    fn it_rejects_the_identity_join() {
        let a = bnd("seq0", Orientation::LowerFlank, 100);
        let b = bnd("seq0", Orientation::HigherFlank, 100);
        let err = Adjacency::try_new_paired(a, b, seq(".")).unwrap_err();
        assert!(matches!(err, Error::IdentityJoin));
    }

    #[test]
    fn it_rejects_a_self_junction() {
        // Fully identical breakends, meaning the same contig, position, and
        // orientation, join a breakend to itself and are rejected, both with an
        // empty and with a non-empty insertion.
        let empty = Adjacency::try_new_paired(
            bnd("seq0", Orientation::LowerFlank, 100),
            bnd("seq0", Orientation::LowerFlank, 100),
            seq("."),
        )
        .unwrap_err();
        assert!(matches!(empty, Error::SelfJunction));

        let inserted = Adjacency::try_new_paired(
            bnd("seq0", Orientation::HigherFlank, 100),
            bnd("seq0", Orientation::HigherFlank, 100),
            seq("AAAC"),
        )
        .unwrap_err();
        assert!(matches!(inserted, Error::SelfJunction));
    }

    fn round_trip(adjacency: Adjacency<dna::Nucleotide>) {
        let rendered = adjacency.to_string();
        let parsed = rendered.parse::<Adjacency<dna::Nucleotide>>().unwrap();
        assert_eq!(parsed, adjacency);
        assert_eq!(parsed.to_string(), rendered);
    }

    #[test]
    fn it_round_trips_a_paired_adjacency_with_an_insertion() {
        let adjacency = Adjacency::try_new_paired(
            bnd("seq0", Orientation::LowerFlank, 100),
            bnd("seq1", Orientation::HigherFlank, 200),
            seq("AT"),
        )
        .unwrap();
        assert_eq!(adjacency.to_string(), "seq0:>:100(i)::seq1:<:200(i)::AT");
        round_trip(adjacency);
    }

    #[test]
    fn it_round_trips_a_paired_adjacency_with_an_empty_insertion() {
        let adjacency = Adjacency::try_new_paired(
            bnd("seq0", Orientation::LowerFlank, 100),
            bnd("seq0", Orientation::HigherFlank, 200),
            seq("."),
        )
        .unwrap();
        assert_eq!(adjacency.to_string(), "seq0:>:100(i)::seq0:<:200(i)::.");
        round_trip(adjacency);
    }

    #[test]
    fn it_round_trips_a_single_ended_adjacency_with_an_insertion() {
        let adjacency = Adjacency::new_single(bnd("seq0", Orientation::LowerFlank, 100), seq("AT"));
        assert_eq!(adjacency.to_string(), "seq0:>:100(i)::.::AT");
        round_trip(adjacency);
    }

    #[test]
    fn it_round_trips_a_single_ended_adjacency_with_an_empty_insertion() {
        let adjacency = Adjacency::new_single(bnd("seq0", Orientation::LowerFlank, 100), seq("."));
        assert_eq!(adjacency.to_string(), "seq0:>:100(i)::.::.");
        round_trip(adjacency);
    }

    #[test]
    fn it_rejects_too_few_fields() {
        let err = "seq0:>:100(i)::AT"
            .parse::<Adjacency<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, ParseError::Format(_)));
    }

    #[test]
    fn it_rejects_too_many_fields() {
        let err = "seq0:>:100(i)::seq1:<:200(i)::AT::extra"
            .parse::<Adjacency<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, ParseError::Format(_)));
    }

    #[test]
    fn it_rejects_an_identity_join_from_a_string() {
        // Same locus, opposite orientations, and an empty insertion, so the
        // parsed paired join encodes no change and must be rejected.
        let err = "seq0:>:100(i)::seq0:<:100(i)::."
            .parse::<Adjacency<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, ParseError::Construct(Error::IdentityJoin)));
    }

    #[test]
    fn it_rejects_a_bad_breakend() {
        let err = "bad::seq1:<:200(i)::AT"
            .parse::<Adjacency<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, ParseError::Breakend(_)));
    }
}
