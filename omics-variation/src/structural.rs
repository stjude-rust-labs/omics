//! Structural variants.
//!
//! Structural variants are modeled as novel adjacencies between breakends,
//! following the breakend representation used by the VCF breakend spec and by
//! GA4GH VRS. A [`Breakend`] is one endpoint of a junction,
//! an [`Adjacency`] is a single junction joining two
//! breakends with an optional non-templated insertion, and a
//! [`StructuralVariant`] is an event built from one or more adjacencies. The
//! event [`Kind`] is derived from the geometry of those adjacencies rather than
//! stored.
//!
//! # Coordinate-anchored orientation
//!
//! A breakend carries a contig, an interbase position, and an
//! [`Orientation`]. It deliberately does not carry a
//! strand. Every position is a coordinate on the reference, so
//! [`LowerFlank`](orientation::Orientation::LowerFlank) and
//! [`HigherFlank`](orientation::Orientation::HigherFlank) name the retained
//! side purely by reference coordinate, and "lower" can only ever mean a
//! smaller reference position.
//!
//! This deliberately splits two facts that other breakend representations
//! bundle into a single per-endpoint value. The first fact is which side of the
//! cut is retained, which lives on the breakend as its
//! [`Orientation`]. The second fact is the reading
//! direction of the retained piece in the derivative molecule, or equivalently
//! whether the join reverse-complements one side. That second fact is stored on
//! no breakend. It emerges from the pair. Two breakends with opposite
//! orientations join co-linearly with no flip, while two breakends with
//! matching orientations form a fold-back that reverse-complements one side,
//! which is how an inversion arises. The [`Adjacency`]
//! carries the residual bookkeeping by storing the non-templated insertion in
//! the reading frame of its canonical first breakend and reverse-complementing
//! it when the supplied breakends are swapped into canonical order.
//!
//! Naming the variants after reference coordinates rather than after a strand
//! keeps the meaning unambiguous. A single fold-back breakend cannot say on its
//! own which way its retained piece reads, so a lone fold-back classifies as
//! [`Kind::Complex`] and only a matched pair resolves an inversion.
//!
//! # Double-stranded reference
//!
//! A breakend breaks double-stranded DNA, and a junction fuses one
//! double-stranded end to another, so both strands join at every adjacency.
//! Only the reference forward strand is written down. Its antiparallel partner
//! is left implicit and fuses along with it by Watson-Crick complementarity.
//!
//! This is why an orientation pair is enough to place a junction. Once each
//! side states which flank it retains and whether the join is co-linear or a
//! fold-back, the way both strands connect is fixed, so the partner-strand
//! connection is forced rather than separately recorded. A fold-back is exactly
//! the case where the forward strand of one side continues onto the partner
//! strand of the other, which is why that side is read reverse-complemented.
//!
//! The model therefore assumes a clean fusion of double-stranded ends,
//! optionally with a non-templated insertion between them. Single-stranded
//! junctions, and events whose two strands join different partners, are out of
//! scope. Both strands fusing at one junction is also a separate matter from
//! whether a whole event is balanced. A reciprocal event that yields more than
//! one derivative molecule is several adjacencies grouped in a
//! [`StructuralVariant`].
//!
//! # Reading the kinds
//!
//! The figures below show how each [`Kind`] looks in the data model. In them
//! `>` marks a lower-flank breakend and `<` marks a higher-flank one, the same
//! glyphs the serialized form uses. Positions are interbase boundaries, `::`
//! separates the fields of one adjacency, and `;` separates the adjacencies of
//! one event. A `.` is an empty insertion or, in the middle field, the open
//! side of a single-ended breakend.
//!
//! A deletion is one co-linear junction whose lower breakend keeps its lower
//! flank, so the two outer pieces join and the middle is dropped.
//!
//! ```text
//! ref   ==A==|100 .......... 5000|==B==
//! adj   seq0:>:100(i)::seq0:<:5000(i)::.
//! =>    ==A====B==                 Deletion  (101..5000 removed)
//! ```
//!
//! A large insertion is one co-linear junction at a single boundary carrying
//! the novel bases.
//!
//! ```text
//! ref   ==A==|100|==B==
//! adj   seq0:>:100(i)::seq0:<:100(i)::ACGT
//! =>    ==A==[ACGT]==B==            Insertion
//! ```
//!
//! A tandem duplication is one co-linear junction whose lower breakend keeps
//! its higher flank, so the join runs backward and the segment repeats.
//!
//! ```text
//! ref   ==| 100 ==seg== 160 |==
//! adj   seq0:<:100(i)::seq0:>:160(i)::GAT
//! =>    ==|100 ==seg== 160||100 ==seg== 160|==   TandemDuplication
//! ```
//!
//! An inversion is two fold-back junctions over the same two boundaries, one
//! joining the lower flanks and one joining the higher flanks.
//!
//! ```text
//! ref   ==| 100 ==seg== 300 |==
//! adj   seq0:>:100(i)::seq0:>:300(i)::. ; seq0:<:100(i)::seq0:<:300(i)::.
//! =>    ==| revcomp(seg) |==            Inversion
//! ```
//!
//! An interchromosomal translocation is one junction across two contigs. With
//! opposite orientations the piece brought over from the other contig keeps its
//! direction. Below, `[R S]` is the segment translocated from `seq1` into the
//! `seq0`-derived molecule.
//!
//! ```text
//! seq0   ==P==Q==|100 . . .            keep P Q up to boundary 100
//! seq1        . . . 900|==R==S==       keep R S from boundary 900
//! adj    seq0:>:100(i)::seq1:<:900(i)::.
//! =>     ==P==Q==[R==S]==              [R S] is spliced in from seq1
//!        Translocation { Interchromosomal, CoLinear }
//! ```
//!
//! With matching orientations the segment brought over is reverse-complemented,
//! a fold-back. Below, `[S' R']` is that same `seq1` segment, now inverted.
//!
//! ```text
//! seq0   ==P==Q==|100 . . .            keep P Q up to boundary 100
//! seq1        . . . |200 R==S==        matching orientation flips this piece
//! adj    seq0:>:100(i)::seq1:>:200(i)::.
//! =>     ==P==Q==[S'==R']==            [S' R'] is [R S] from seq1, reverse-complemented
//!        Translocation { Interchromosomal, FoldBack }
//! ```
//!
//! An intrachromosomal translocation is the balanced forward relocation of a
//! segment, three co-linear junctions over three boundaries where each boundary
//! appears once as a lower flank and once as a higher flank.
//!
//! ```text
//! adj   seq0:>:100(i)::seq0:<:200(i)::.     origin flanks rejoin
//!       seq0:<:100(i)::seq0:>:400(i)::.     segment start meets the target
//!       seq0:>:200(i)::seq0:<:400(i)::.     segment end meets the target
//! =>    Translocation { Intrachromosomal, CoLinear }
//! ```
//!
//! A single-ended breakend joins one locus to novel or unassembled sequence,
//! with the other side open.
//!
//! ```text
//! ref   ==A==|100  ->  novel sequence
//! adj   seq0:>:100(i)::.::AT
//! =>    Breakend
//! ```
//!
//! Any well-formed event whose shape matches none of these, such as a lone
//! fold-back or a reciprocal interchromosomal pair, is [`Kind::Complex`].

pub mod adjacency;
pub mod breakend;
pub mod orientation;

use std::fmt;
use std::str::FromStr;

use omics_coordinate::position::Number;
use omics_molecule::compound::Complement;
use omics_molecule::compound::Nucleotide;
use thiserror::Error;

use crate::structural::adjacency::Adjacency;
use crate::structural::breakend::Breakend;
use crate::structural::orientation::Orientation;

/// The separator between the adjacency tokens of a [`StructuralVariant`].
const EVENT_SEPARATOR: &str = ";";

/// Whether a translocation stays on one contig or crosses contigs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Locality {
    /// The translocation crosses two contigs.
    Interchromosomal,

    /// The translocation stays on one contig.
    Intrachromosomal,
}

/// How the two sides of a translocation join with respect to reading direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Join {
    /// The two sides read in the same direction, joined without a flip. Its
    /// breakends have opposite orientations.
    CoLinear,

    /// One side is reverse-complemented, so the join folds back. Its breakends
    /// have matching orientations, meaning the fused piece is inverted.
    FoldBack,
}

/// The derived class of a structural variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// A large deletion.
    Deletion,

    /// A large insertion.
    Insertion,

    /// An internal tandem duplication.
    TandemDuplication,

    /// An inversion.
    Inversion,

    /// A translocation, with its locality and how its sides join.
    Translocation {
        /// Whether the translocation crosses contigs or stays on one.
        locality: Locality,

        /// Whether the sides join co-linearly or fold back.
        join: Join,
    },

    /// A single-ended breakend joined to novel sequence.
    Breakend,

    /// A well-formed event whose shape matches no named pattern.
    Complex,
}

/// An error related to constructing a [`StructuralVariant`].
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// The event held no adjacencies.
    #[error("a structural variant must hold at least one adjacency")]
    Empty,
}

/// A parse error related to a [`StructuralVariant`].
#[derive(Error, Debug)]
pub enum ParseError {
    /// The structural-variant string was empty.
    #[error("a structural variant string must not be empty")]
    Empty,

    /// An adjacency token failed to parse.
    #[error(transparent)]
    Adjacency(#[from] adjacency::ParseError),

    /// The structural variant could not be constructed.
    #[error(transparent)]
    Construct(#[from] Error),
}

/// A structural variant, built from one or more novel adjacencies.
///
/// The adjacencies are held in a normalized form, sorted by a deterministic
/// canonical key and deduplicated, so that equality and classification do not
/// depend on the input order and are insensitive to duplicates.
///
/// `Hash` is intentionally not derived because `Adjacency` does not implement
/// it.
///
/// A structural variant serializes as its adjacency tokens joined with `;`.
///
/// # Examples
///
/// ```
/// use omics_molecule::polymer::dna;
/// use omics_variation::structural::StructuralVariant;
///
/// let rendered = "seq0:>:100(i)::seq0:<:200(i)::.";
/// let variant = rendered.parse::<StructuralVariant<dna::Nucleotide>>()?;
/// assert_eq!(variant.to_string(), rendered);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuralVariant<N: Nucleotide> {
    /// The adjacencies making up the event, held in normalized order.
    adjacencies: Vec<Adjacency<N>>,
}

impl<N: Nucleotide> StructuralVariant<N> {
    /// Attempts to create a new [`StructuralVariant`].
    ///
    /// The supplied adjacencies are normalized, meaning they are sorted by a
    /// deterministic canonical key and deduplicated. As a result, a permuted or
    /// duplicate-bearing input compares equal to its canonical form and
    /// classifies identically.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::StructuralVariant;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::try_new_paired(a, b, ".".parse()?)?;
    /// let variant = StructuralVariant::try_new(vec![adjacency])?;
    /// assert_eq!(variant.adjacencies().len(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(mut adjacencies: Vec<Adjacency<N>>) -> Result<Self, Error> {
        if adjacencies.is_empty() {
            return Err(Error::Empty);
        }

        adjacencies.sort_by_key(|adjacency| adjacency_key(adjacency));
        adjacencies.dedup();

        Ok(Self { adjacencies })
    }

    /// Gets the adjacencies of this structural variant.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::StructuralVariant;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::try_new_paired(a, b, ".".parse()?)?;
    /// let variant = StructuralVariant::try_new(vec![adjacency])?;
    /// assert_eq!(variant.adjacencies().len(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn adjacencies(&self) -> &[Adjacency<N>] {
        &self.adjacencies
    }

    /// Derives the [`Kind`] of this structural variant from its geometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::structural::Kind;
    /// use omics_variation::structural::StructuralVariant;
    /// use omics_variation::structural::adjacency::Adjacency;
    /// use omics_variation::structural::breakend::Breakend;
    /// use omics_variation::structural::orientation::Orientation;
    ///
    /// let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
    /// let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
    /// let adjacency = Adjacency::<dna::Nucleotide>::try_new_paired(a, b, ".".parse()?)?;
    /// let variant = StructuralVariant::try_new(vec![adjacency])?;
    /// assert_eq!(variant.kind(), Kind::Deletion);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn kind(&self) -> Kind {
        match self.adjacencies.as_slice() {
            [single] => classify_single(single),
            [first, second] => classify_pair(first, second),
            [first, second, third] => classify_triple(first, second, third),
            _ => Kind::Complex,
        }
    }
}

impl<N: Nucleotide> fmt::Display for StructuralVariant<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let tokens = self
            .adjacencies
            .iter()
            .map(|adjacency| adjacency.to_string())
            .collect::<Vec<_>>();
        write!(f, "{}", tokens.join(EVENT_SEPARATOR))
    }
}

impl<N: Nucleotide + Complement> FromStr for StructuralVariant<N> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let adjacencies = s
            .split(EVENT_SEPARATOR)
            .map(|token| token.parse::<Adjacency<N>>())
            .collect::<Result<Vec<_>, _>>()?;

        Ok(StructuralVariant::try_new(adjacencies)?)
    }
}

/// An owned canonical ordering key for a single breakend.
///
/// It mirrors [`Breakend::canonical_key`] but owns its contig so it can back a
/// sortable key for a whole adjacency.
type BreakendKey = (String, Number, u8);

/// Gets the owned canonical ordering key for a breakend.
fn breakend_key(breakend: &Breakend) -> BreakendKey {
    let (contig, position, rank) = breakend.canonical_key();
    (contig.to_string(), position, rank)
}

/// Gets the deterministic canonical ordering key for an adjacency.
///
/// The key tags paired adjacencies ahead of single ones, then orders by the
/// canonical breakend keys and finally by the insertion sequence. It backs the
/// normalized ordering of the adjacencies in a [`StructuralVariant`].
fn adjacency_key<N: Nucleotide>(
    adjacency: &Adjacency<N>,
) -> (u8, BreakendKey, BreakendKey, String) {
    match adjacency {
        Adjacency::Paired(paired) => (
            0,
            breakend_key(paired.a()),
            breakend_key(paired.b()),
            paired.insertion().to_string(),
        ),
        Adjacency::Single(single) => (
            1,
            breakend_key(single.breakend()),
            (String::new(), 0, 0),
            single.insertion().to_string(),
        ),
    }
}

/// Gets the two breakends of a paired adjacency, or `None` for a single one.
fn paired_breakends<N: Nucleotide>(adjacency: &Adjacency<N>) -> Option<(&Breakend, &Breakend)> {
    match adjacency {
        Adjacency::Paired(paired) => Some((paired.a(), paired.b())),
        Adjacency::Single(_) => None,
    }
}

/// Reports whether every breakend across the adjacencies is on one contig.
///
/// A single-ended adjacency has no second known locus, so its presence makes
/// the answer `false`.
fn one_contig<N: Nucleotide>(adjacencies: &[&Adjacency<N>]) -> bool {
    let mut contig: Option<&str> = None;
    for adjacency in adjacencies {
        let Some((a, b)) = paired_breakends(adjacency) else {
            return false;
        };
        for breakend in [a, b] {
            match contig {
                None => contig = Some(breakend.contig().as_str()),
                Some(seen) if seen != breakend.contig().as_str() => return false,
                Some(_) => {}
            }
        }
    }
    true
}

/// Gets the shared orientation of a fold-back adjacency, or `None` otherwise.
///
/// A fold-back is a paired adjacency whose two breakends share an orientation.
fn fold_back_orientation<N: Nucleotide>(adjacency: &Adjacency<N>) -> Option<Orientation> {
    match paired_breakends(adjacency) {
        Some((a, b)) if a.orientation() == b.orientation() => Some(a.orientation()),
        _ => None,
    }
}

/// Gets the sorted pair of interbase positions of a paired adjacency.
fn position_pair<N: Nucleotide>(adjacency: &Adjacency<N>) -> Option<(Number, Number)> {
    let (a, b) = paired_breakends(adjacency)?;
    let (lo, hi) = (a.position().get(), b.position().get());
    Some(if lo <= hi { (lo, hi) } else { (hi, lo) })
}

/// Classifies a structural variant made of exactly one adjacency.
fn classify_single<N: Nucleotide>(adjacency: &Adjacency<N>) -> Kind {
    let Adjacency::Paired(paired) = adjacency else {
        // A single-ended breakend.
        return Kind::Breakend;
    };
    let (a, b) = (paired.a(), paired.b());

    if a.contig() != b.contig() {
        let join = if a.orientation().is_opposite(b.orientation()) {
            Join::CoLinear
        } else {
            Join::FoldBack
        };
        return Kind::Translocation {
            locality: Locality::Interchromosomal,
            join,
        };
    }

    if !a.orientation().is_opposite(b.orientation()) {
        // A lone fold-back is an incomplete inversion.
        return Kind::Complex;
    }

    match a.orientation() {
        Orientation::HigherFlank => Kind::TandemDuplication,
        Orientation::LowerFlank => {
            if a.position().get() == b.position().get() {
                Kind::Insertion
            } else {
                Kind::Deletion
            }
        }
    }
}

/// Classifies a two-adjacency event.
///
/// An inversion is exactly two paired fold-back adjacencies on one contig over
/// the same two boundary positions, one both-`LowerFlank` and one
/// both-`HigherFlank`. Anything else is [`Kind::Complex`].
fn classify_pair<N: Nucleotide>(first: &Adjacency<N>, second: &Adjacency<N>) -> Kind {
    if !one_contig(&[first, second]) {
        return Kind::Complex;
    }

    let orientations = (fold_back_orientation(first), fold_back_orientation(second));
    let positions = (position_pair(first), position_pair(second));
    if let ((Some(one), Some(two)), (Some(p1), Some(p2))) = (orientations, positions) {
        // The shared boundary pair must span two distinct positions. A pair
        // with `lo == hi` is a zero-length fold-back that cannot form a real
        // inversion, so it falls through to `Kind::Complex`. With the
        // self-junction rejection in `Adjacency::try_new_paired`, such a pair is
        // already unconstructable, so this guard is defense-in-depth.
        let (lo, hi) = p1;
        if one.is_opposite(two) && p1 == p2 && lo != hi {
            return Kind::Inversion;
        }
    }

    Kind::Complex
}

/// Classifies a three-adjacency event.
///
/// The only recognized signature is the balanced forward intrachromosomal
/// relocation, namely three paired co-linear (opposite-orientation) junctions
/// on one contig forming a triangle over exactly three distinct boundaries,
/// where each boundary appears exactly once as a `LowerFlank` and exactly once
/// as a `HigherFlank` across the three junctions. Anything else, including a
/// triple that reuses a flank or one with a fold-back, is [`Kind::Complex`].
fn classify_triple<N: Nucleotide>(
    first: &Adjacency<N>,
    second: &Adjacency<N>,
    third: &Adjacency<N>,
) -> Kind {
    let adjacencies = [first, second, third];
    if !one_contig(&adjacencies) {
        return Kind::Complex;
    }

    // Every junction must be a co-linear (opposite-orientation) paired join
    // that is not a self-loop. A self-loop (both breakends at the same
    // boundary) is a degenerate join that cannot participate in a connected
    // triangle, so any such junction demotes the triple to complex even when
    // the per-boundary flank tally would otherwise be satisfied.
    for adjacency in adjacencies {
        match paired_breakends(adjacency) {
            Some((a, b))
                if a.orientation().is_opposite(b.orientation())
                    && a.position().get() != b.position().get() => {}
            _ => return Kind::Complex,
        }
    }

    // Tally, per boundary position, how often it appears as a `LowerFlank` and
    // as a `HigherFlank`. The triangle requires exactly three distinct
    // boundaries, each appearing once as each flank.
    let mut counts: Vec<(Number, usize, usize)> = Vec::new();
    for adjacency in adjacencies {
        if let Some((a, b)) = paired_breakends(adjacency) {
            for breakend in [a, b] {
                let position = breakend.position().get();
                if !counts.iter().any(|(p, ..)| *p == position) {
                    counts.push((position, 0, 0));
                }
                let entry = counts
                    .iter_mut()
                    .find(|(p, ..)| *p == position)
                    .expect("the position was just inserted");
                match breakend.orientation() {
                    Orientation::LowerFlank => entry.1 += 1,
                    Orientation::HigherFlank => entry.2 += 1,
                }
            }
        }
    }

    let triangle = counts.len() == 3
        && counts
            .iter()
            .all(|(_, lower, higher)| *lower == 1 && *higher == 1);

    if triangle {
        // The recognized relocation is three co-linear junctions, so the sides
        // join co-linearly.
        Kind::Translocation {
            locality: Locality::Intrachromosomal,
            join: Join::CoLinear,
        }
    } else {
        Kind::Complex
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::position::Number;
    use omics_molecule::polymer::dna;

    use super::*;
    use crate::structural::adjacency::Adjacency;
    use crate::structural::breakend::Breakend;
    use crate::structural::orientation::Orientation;

    type Sv = StructuralVariant<dna::Nucleotide>;

    fn bnd(orientation: Orientation, position: Number) -> Breakend {
        Breakend::try_new("seq0", orientation, position).unwrap()
    }

    fn paired(x: Breakend, y: Breakend, insertion: &str) -> Adjacency<dna::Nucleotide> {
        Adjacency::try_new_paired(x, y, insertion.parse().unwrap()).unwrap()
    }

    #[test]
    fn it_rejects_an_empty_event() {
        assert!(matches!(Sv::try_new(Vec::new()), Err(Error::Empty)));
    }

    #[test]
    fn it_classifies_a_deletion() {
        let adjacency = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        assert_eq!(Sv::try_new(vec![adjacency]).unwrap().kind(), Kind::Deletion);
    }

    #[test]
    fn it_classifies_a_tandem_duplication() {
        let adjacency = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        assert_eq!(
            Sv::try_new(vec![adjacency]).unwrap().kind(),
            Kind::TandemDuplication
        );
    }

    #[test]
    fn it_classifies_a_large_insertion() {
        let adjacency = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 100),
            "AAAA",
        );
        assert_eq!(
            Sv::try_new(vec![adjacency]).unwrap().kind(),
            Kind::Insertion
        );
    }

    #[test]
    fn it_classifies_an_interchromosomal_translocation() {
        let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100).unwrap();
        let b = Breakend::try_new("seq1", Orientation::HigherFlank, 200).unwrap();
        let adjacency = Adjacency::try_new_paired(a, b, ".".parse().unwrap()).unwrap();
        assert_eq!(
            Sv::try_new(vec![adjacency]).unwrap().kind(),
            Kind::Translocation {
                locality: Locality::Interchromosomal,
                join: Join::CoLinear
            }
        );
    }

    #[test]
    fn it_classifies_an_inverted_interchromosomal_translocation() {
        // Matching orientations across contigs fuse the other contig's piece
        // reverse-complemented, so the join folds back.
        let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100).unwrap();
        let b = Breakend::try_new("seq1", Orientation::LowerFlank, 200).unwrap();
        let adjacency = Adjacency::try_new_paired(a, b, ".".parse().unwrap()).unwrap();
        assert_eq!(
            Sv::try_new(vec![adjacency]).unwrap().kind(),
            Kind::Translocation {
                locality: Locality::Interchromosomal,
                join: Join::FoldBack
            }
        );
    }

    #[test]
    fn it_classifies_a_single_ended_breakend() {
        let breakend = bnd(Orientation::LowerFlank, 100);
        let adjacency = Adjacency::new_single(breakend, "AT".parse().unwrap());
        assert_eq!(Sv::try_new(vec![adjacency]).unwrap().kind(), Kind::Breakend);
    }

    #[test]
    fn it_classifies_a_lone_fold_back_as_complex() {
        let adjacency = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        assert_eq!(Sv::try_new(vec![adjacency]).unwrap().kind(), Kind::Complex);
    }

    #[test]
    fn it_classifies_an_inversion() {
        // Two fold-back junctions over the same two boundaries {100, 200}, one
        // both-LowerFlank and one both-HigherFlank.
        let lower = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        let higher = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        assert_eq!(
            Sv::try_new(vec![lower, higher]).unwrap().kind(),
            Kind::Inversion
        );
    }

    #[test]
    fn it_classifies_a_genuine_inversion_over_distinct_boundaries() {
        // Two fold-back junctions over two distinct boundaries {100, 300}, one
        // both-LowerFlank and one both-HigherFlank. The distinct-position guard
        // in `classify_pair` is satisfied, so this classifies as an inversion.
        let lower = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 300),
            ".",
        );
        let higher = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::HigherFlank, 300),
            ".",
        );
        assert_eq!(
            Sv::try_new(vec![lower, higher]).unwrap().kind(),
            Kind::Inversion
        );
    }

    #[test]
    fn it_classifies_an_intrachromosomal_translocation() {
        // Three co-linear junctions over boundaries {100, 200, 400}, forming a
        // triangle where each boundary appears once as each flank.
        let origin = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        let target_left = paired(
            bnd(Orientation::LowerFlank, 400),
            bnd(Orientation::HigherFlank, 100),
            ".",
        );
        let target_right = paired(
            bnd(Orientation::LowerFlank, 200),
            bnd(Orientation::HigherFlank, 400),
            ".",
        );
        assert_eq!(
            Sv::try_new(vec![origin, target_left, target_right])
                .unwrap()
                .kind(),
            Kind::Translocation {
                locality: Locality::Intrachromosomal,
                join: Join::CoLinear
            }
        );
    }

    #[test]
    fn it_classifies_two_unrelated_fold_backs_as_complex() {
        let one = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        let two = paired(
            bnd(Orientation::LowerFlank, 300),
            bnd(Orientation::LowerFlank, 400),
            ".",
        );
        assert_eq!(Sv::try_new(vec![one, two]).unwrap().kind(), Kind::Complex);
    }

    #[test]
    fn it_rejects_a_triple_that_reuses_a_flank_as_complex() {
        // Three co-linear junctions over three distinct boundaries {100, 200,
        // 300}, each appearing twice overall, but 100 appears twice as a
        // LowerFlank and never as a HigherFlank, so it is not a proper triangle.
        let one = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        let two = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 300),
            ".",
        );
        let three = paired(
            bnd(Orientation::LowerFlank, 200),
            bnd(Orientation::HigherFlank, 300),
            ".",
        );
        assert_eq!(
            Sv::try_new(vec![one, two, three]).unwrap().kind(),
            Kind::Complex
        );
    }

    #[test]
    fn it_rejects_three_self_loop_junctions_as_complex() {
        // Three independent self-loop insertions. Each junction's two breakends
        // sit at the same position, so the graph is three disconnected
        // self-loops rather than a connected triangle. Although the per-boundary
        // flank tally is satisfied, this must classify as Complex, not an
        // intrachromosomal translocation.
        let j1 = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 100),
            "A",
        );
        let j2 = paired(
            bnd(Orientation::LowerFlank, 200),
            bnd(Orientation::HigherFlank, 200),
            "A",
        );
        let j3 = paired(
            bnd(Orientation::LowerFlank, 300),
            bnd(Orientation::HigherFlank, 300),
            "A",
        );
        assert_eq!(Sv::try_new(vec![j1, j2, j3]).unwrap().kind(), Kind::Complex);
    }

    #[test]
    fn it_rejects_an_inverted_reinsertion_as_complex() {
        // A three-junction set where one junction is a fold-back (both breakends
        // share the same orientation). The co-linear check rejects any triple
        // containing a fold-back, so this classifies as Complex.
        let fold_back = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        let two = paired(
            bnd(Orientation::LowerFlank, 400),
            bnd(Orientation::HigherFlank, 100),
            ".",
        );
        let three = paired(
            bnd(Orientation::LowerFlank, 200),
            bnd(Orientation::HigherFlank, 400),
            ".",
        );
        assert_eq!(
            Sv::try_new(vec![fold_back, two, three]).unwrap().kind(),
            Kind::Complex
        );
    }

    #[test]
    fn it_normalizes_a_permuted_adjacency_list() {
        let lower = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        let higher = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );

        let canonical = Sv::try_new(vec![lower.clone(), higher.clone()]).unwrap();
        let permuted = Sv::try_new(vec![higher, lower]).unwrap();

        assert_eq!(canonical, permuted);
        assert_eq!(canonical.kind(), permuted.kind());
        assert_eq!(permuted.kind(), Kind::Inversion);
    }

    fn round_trip(variant: Sv, expected: Kind) {
        assert_eq!(variant.kind(), expected);
        let rendered = variant.to_string();
        let parsed = rendered.parse::<Sv>().unwrap();
        assert_eq!(parsed, variant);
        assert_eq!(parsed.kind(), variant.kind());
        assert_eq!(parsed.to_string(), rendered);
    }

    #[test]
    fn it_round_trips_a_large_deletion() {
        let adjacency = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        round_trip(Sv::try_new(vec![adjacency]).unwrap(), Kind::Deletion);
    }

    #[test]
    fn it_round_trips_a_large_insertion() {
        let adjacency = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 100),
            "AAAA",
        );
        round_trip(Sv::try_new(vec![adjacency]).unwrap(), Kind::Insertion);
    }

    #[test]
    fn it_round_trips_an_interchromosomal_translocation() {
        let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100).unwrap();
        let b = Breakend::try_new("seq1", Orientation::HigherFlank, 200).unwrap();
        let adjacency = Adjacency::try_new_paired(a, b, ".".parse().unwrap()).unwrap();
        round_trip(
            Sv::try_new(vec![adjacency]).unwrap(),
            Kind::Translocation {
                locality: Locality::Interchromosomal,
                join: Join::CoLinear,
            },
        );
    }

    #[test]
    fn it_round_trips_an_intrachromosomal_translocation() {
        let origin = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        let target_left = paired(
            bnd(Orientation::LowerFlank, 400),
            bnd(Orientation::HigherFlank, 100),
            ".",
        );
        let target_right = paired(
            bnd(Orientation::LowerFlank, 200),
            bnd(Orientation::HigherFlank, 400),
            ".",
        );
        round_trip(
            Sv::try_new(vec![origin, target_left, target_right]).unwrap(),
            Kind::Translocation {
                locality: Locality::Intrachromosomal,
                join: Join::CoLinear,
            },
        );
    }

    #[test]
    fn it_round_trips_an_inversion() {
        let lower = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        let higher = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );
        round_trip(Sv::try_new(vec![lower, higher]).unwrap(), Kind::Inversion);
    }

    #[test]
    fn it_round_trips_an_internal_tandem_duplication() {
        let adjacency = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        round_trip(
            Sv::try_new(vec![adjacency]).unwrap(),
            Kind::TandemDuplication,
        );
    }

    #[test]
    fn it_round_trips_a_single_ended_breakend() {
        let adjacency =
            Adjacency::new_single(bnd(Orientation::LowerFlank, 100), "AT".parse().unwrap());
        round_trip(Sv::try_new(vec![adjacency]).unwrap(), Kind::Breakend);
    }

    #[test]
    fn it_rejects_an_empty_string() {
        let err = "".parse::<Sv>().unwrap_err();
        assert!(matches!(err, ParseError::Empty));
    }

    #[test]
    fn it_rejects_a_bad_adjacency_token() {
        let err = "seq0:>:100(i)::AT".parse::<Sv>().unwrap_err();
        assert!(matches!(err, ParseError::Adjacency(_)));
    }

    #[test]
    fn it_normalizes_a_duplicate_adjacency_list() {
        let lower = paired(
            bnd(Orientation::LowerFlank, 100),
            bnd(Orientation::LowerFlank, 200),
            ".",
        );
        let higher = paired(
            bnd(Orientation::HigherFlank, 100),
            bnd(Orientation::HigherFlank, 200),
            ".",
        );

        let canonical = Sv::try_new(vec![lower.clone(), higher.clone()]).unwrap();
        let with_duplicate = Sv::try_new(vec![lower, higher.clone(), higher]).unwrap();

        assert_eq!(canonical, with_duplicate);
        assert_eq!(canonical.adjacencies().len(), 2);
        assert_eq!(with_duplicate.adjacencies().len(), 2);
        assert_eq!(canonical.kind(), with_duplicate.kind());
        assert_eq!(with_duplicate.kind(), Kind::Inversion);
    }
}
