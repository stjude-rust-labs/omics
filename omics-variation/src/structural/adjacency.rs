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

/// An error related to constructing a [`PairedAdjacency`].
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

/// A junction joining two oriented breakends, held in canonical form.
///
/// The two breakends are ordered by their canonical key, and the non-templated
/// insertion is carried in the reading frame of the canonical lower breakend,
/// so that equality and later classification do not depend on the input order.
/// A value of this type can only be built through [`PairedAdjacency::try_new`],
/// so it always upholds those invariants.
///
/// `Hash` is intentionally not derived because `Sequence` does not implement
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairedAdjacency<N: Nucleotide> {
    /// The canonical lower-locus breakend.
    a: Breakend,

    /// The canonical higher-locus breakend.
    b: Breakend,

    /// The non-templated insertion, in the reading frame of `a`.
    insertion: Sequence<N>,
}

impl<N: Nucleotide + Complement> PairedAdjacency<N> {
    /// Attempts to create a [`PairedAdjacency`] in canonical form.
    ///
    /// The two breakends are ordered by their canonical key, which sorts by
    /// contig, then interbase position, then orientation. The `insertion` is
    /// supplied as it reads across the junction from the first argument `x`
    /// toward the second argument `y`. If `x` already sorts at or before `y`,
    /// both are kept as given. If `x` sorts after `y`, the pair is swapped into
    /// canonical order and the insertion is reverse-complemented so that it
    /// still reads from the canonical lower breakend toward the higher one. The
    /// identity join, meaning two breakends at the same locus with opposite
    /// orientations and an empty insertion, is rejected because it encodes no
    /// change. Two fully identical breakends are rejected because they join a
    /// breakend to itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::adjacency::PairedAdjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let paired = PairedAdjacency::<dna::Nucleotide>::try_new(a, b, ".".parse()?)?;
    /// assert_eq!(paired.a().position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(x: Breakend, y: Breakend, insertion: Sequence<N>) -> Result<Self, Error> {
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

        Ok(Self { a, b, insertion })
    }
}

impl<N: Nucleotide> PairedAdjacency<N> {
    /// Gets the canonical lower-locus breakend.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna;
    /// # use omics_variation::structural::adjacency::PairedAdjacency;
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let paired = PairedAdjacency::<dna::Nucleotide>::try_new(a, b, ".".parse()?)?;
    /// assert_eq!(paired.a().position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn a(&self) -> &Breakend {
        &self.a
    }

    /// Gets the canonical higher-locus breakend.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna;
    /// # use omics_variation::structural::adjacency::PairedAdjacency;
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let paired = PairedAdjacency::<dna::Nucleotide>::try_new(a, b, ".".parse()?)?;
    /// assert_eq!(paired.b().position().get(), 200);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn b(&self) -> &Breakend {
        &self.b
    }

    /// Gets the non-templated insertion, in the reading frame of `a`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna;
    /// # use omics_variation::structural::adjacency::PairedAdjacency;
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq1", Orientation::HigherFlank, 200)?;
    /// let paired = PairedAdjacency::<dna::Nucleotide>::try_new(a, b, "AT".parse()?)?;
    /// assert_eq!(paired.insertion().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn insertion(&self) -> &Sequence<N> {
        &self.insertion
    }
}

impl<N: Nucleotide> fmt::Display for PairedAdjacency<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{ADJACENCY_SEPARATOR}{}{ADJACENCY_SEPARATOR}{}",
            self.a, self.b, self.insertion
        )
    }
}

/// A junction joining one breakend to novel or unassembled sequence.
///
/// One side is a known breakend and the other is open, so the adjacency carries
/// a single locus and the non-templated insertion at its open side.
///
/// `Hash` is intentionally not derived because `Sequence` does not implement
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SingleAdjacency<N: Nucleotide> {
    /// The single known breakend.
    breakend: Breakend,

    /// The non-templated insertion at the open side.
    insertion: Sequence<N>,
}

impl<N: Nucleotide> SingleAdjacency<N> {
    /// Creates a [`SingleAdjacency`].
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::adjacency::SingleAdjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let single = SingleAdjacency::<dna::Nucleotide>::new(breakend, "AT".parse()?);
    /// assert_eq!(single.insertion().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(breakend: Breakend, insertion: Sequence<N>) -> Self {
        Self {
            breakend,
            insertion,
        }
    }

    /// Gets the single known breakend.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna;
    /// # use omics_variation::structural::adjacency::SingleAdjacency;
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let single = SingleAdjacency::<dna::Nucleotide>::new(breakend, "AT".parse()?);
    /// assert_eq!(single.breakend().position().get(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn breakend(&self) -> &Breakend {
        &self.breakend
    }

    /// Gets the non-templated insertion at the open side.
    ///
    /// # Examples
    ///
    /// ```
    /// # use omics_molecule::polymer::dna;
    /// # use omics_variation::structural::adjacency::SingleAdjacency;
    /// # use omics_variation::structural::breakend::Breakend;
    /// # use omics_variation::structural::orientation::Orientation;
    /// let breakend = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let single = SingleAdjacency::<dna::Nucleotide>::new(breakend, "AT".parse()?);
    /// assert_eq!(single.insertion().to_string(), "AT");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn insertion(&self) -> &Sequence<N> {
        &self.insertion
    }
}

impl<N: Nucleotide> fmt::Display for SingleAdjacency<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{ADJACENCY_SEPARATOR}{OPEN_SIDE}{ADJACENCY_SEPARATOR}{}",
            self.breakend, self.insertion
        )
    }
}

/// A single novel junction between one or two oriented breakends.
///
/// A [`Paired`](Adjacency::Paired) junction joins two known loci, while a
/// [`Single`](Adjacency::Single) junction joins one known locus to novel or
/// unassembled sequence. Each variant wraps a strong type that owns its
/// invariants, so a paired junction is always canonical.
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
pub enum Adjacency<N: Nucleotide> {
    /// Two loci joined, with optional non-templated bases between them.
    Paired(PairedAdjacency<N>),

    /// One locus joined to novel or unassembled sequence, open on the other
    /// side.
    Single(SingleAdjacency<N>),
}

impl<N: Nucleotide + Complement> Adjacency<N> {
    /// Attempts to create a paired [`Adjacency`] in canonical form.
    ///
    /// This is a convenience wrapper over [`PairedAdjacency::try_new`].
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
    /// assert!(matches!(adjacency, Adjacency::Paired(_)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new_paired(x: Breakend, y: Breakend, insertion: Sequence<N>) -> Result<Self, Error> {
        Ok(Adjacency::Paired(PairedAdjacency::try_new(
            x, y, insertion,
        )?))
    }
}

impl<N: Nucleotide> Adjacency<N> {
    /// Creates a single-ended [`Adjacency`].
    ///
    /// This is a convenience wrapper over [`SingleAdjacency::new`].
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
    /// assert!(matches!(adjacency, Adjacency::Single(_)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new_single(breakend: Breakend, insertion: Sequence<N>) -> Self {
        Adjacency::Single(SingleAdjacency::new(breakend, insertion))
    }
}

impl<N: Nucleotide> fmt::Display for Adjacency<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Adjacency::Paired(paired) => write!(f, "{paired}"),
            Adjacency::Single(single) => write!(f, "{single}"),
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
        let Adjacency::Paired(paired) = adjacency else {
            panic!("expected a paired adjacency");
        };
        assert_eq!(paired.insertion().to_string(), "GTTT");
    }

    #[test]
    fn it_builds_a_single_adjacency() {
        let breakend = bnd("seq0", Orientation::LowerFlank, 100);
        let adjacency = Adjacency::new_single(breakend, seq("AT"));
        let Adjacency::Single(single) = adjacency else {
            panic!("expected a single adjacency");
        };
        assert_eq!(single.insertion().to_string(), "AT");
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
