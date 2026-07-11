//! The structural-variant tier.

pub mod adjacency;
pub mod breakend;
pub mod orientation;

use omics_coordinate::position::Number;
use omics_molecule::compound::Nucleotide;
use thiserror::Error;

use crate::structural::adjacency::Adjacency;
use crate::structural::breakend::Breakend;
use crate::structural::orientation::Orientation;

/// Whether a translocation stays on one contig or crosses contigs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Locality {
    /// The translocation crosses two contigs.
    Interchromosomal,

    /// The translocation stays on one contig.
    Intrachromosomal,
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

    /// A translocation, with its locality.
    Translocation(Locality),

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

/// A structural variant, built from one or more novel adjacencies.
///
/// The adjacencies are held in a normalized form, sorted by a deterministic
/// canonical key and deduplicated, so that equality and classification do not
/// depend on the input order and are insensitive to duplicates.
///
/// `Hash` is intentionally not derived because `Adjacency` does not implement
/// it.
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
fn adjacency_key<N: Nucleotide>(adjacency: &Adjacency<N>) -> (u8, BreakendKey, BreakendKey, String) {
    if let Some((a, b, insertion)) = adjacency.paired() {
        (0, breakend_key(a), breakend_key(b), insertion.to_string())
    } else if let Some((breakend, insertion)) = adjacency.single() {
        (
            1,
            breakend_key(breakend),
            (String::new(), 0, 0),
            insertion.to_string(),
        )
    } else {
        unreachable!("an adjacency is always paired or single")
    }
}

/// Gets the two breakends of a paired adjacency, or `None` for a single one.
fn paired_breakends<N: Nucleotide>(adjacency: &Adjacency<N>) -> Option<(&Breakend, &Breakend)> {
    adjacency.paired().map(|(a, b, _)| (a, b))
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
    let Some((a, b, _)) = adjacency.paired() else {
        // A single-ended breakend.
        return Kind::Breakend;
    };

    if a.contig() != b.contig() {
        return Kind::Translocation(Locality::Interchromosomal);
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
        if one.is_opposite(two) && p1 == p2 {
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

    // Every junction must be a co-linear (opposite-orientation) paired join.
    for adjacency in adjacencies {
        match paired_breakends(adjacency) {
            Some((a, b)) if a.orientation().is_opposite(b.orientation()) => {}
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
                if !counts.iter().any(|(p, _, _)| *p == position) {
                    counts.push((position, 0, 0));
                }
                let entry = counts
                    .iter_mut()
                    .find(|(p, _, _)| *p == position)
                    .expect("the position was just inserted");
                match breakend.orientation() {
                    Orientation::LowerFlank => entry.1 += 1,
                    Orientation::HigherFlank => entry.2 += 1,
                }
            }
        }
    }

    let triangle =
        counts.len() == 3 && counts.iter().all(|(_, lower, higher)| *lower == 1 && *higher == 1);

    if triangle {
        Kind::Translocation(Locality::Intrachromosomal)
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
        assert_eq!(Sv::try_new(vec![adjacency]).unwrap().kind(), Kind::Insertion);
    }

    #[test]
    fn it_classifies_an_interchromosomal_translocation() {
        let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100).unwrap();
        let b = Breakend::try_new("seq1", Orientation::HigherFlank, 200).unwrap();
        let adjacency = Adjacency::try_new_paired(a, b, ".".parse().unwrap()).unwrap();
        assert_eq!(
            Sv::try_new(vec![adjacency]).unwrap().kind(),
            Kind::Translocation(Locality::Interchromosomal)
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
            Kind::Translocation(Locality::Intrachromosomal)
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
