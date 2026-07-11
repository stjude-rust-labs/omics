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
            // Multi-adjacency classification is added in the next task.
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
}
