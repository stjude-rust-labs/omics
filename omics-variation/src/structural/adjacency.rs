//! Novel adjacencies between breakends.

use omics_molecule::compound::Complement;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;
use thiserror::Error;

use crate::structural::breakend::Breakend;

/// An error related to constructing an [`Adjacency`].
#[derive(Error, Debug, PartialEq, Eq)]
pub enum Error {
    /// The two breakends sit at the same locus with opposite orientations and
    /// an empty insertion, which encodes no change.
    #[error("adjacency is an identity join and encodes no change")]
    IdentityJoin,
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
/// with its breakends ordered by [`Breakend::canonical_key`] and its insertion
/// carried in the reading frame of the canonical lower breakend, so that
/// equality and later classification do not depend on the input order.
///
/// `Hash` is intentionally not derived because `Sequence` does not implement
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adjacency<N: Nucleotide> {
    /// The private, canonical inner representation.
    inner: Inner<N>,
}

impl<N: Nucleotide + Complement> Adjacency<N> {
    /// Attempts to create a paired [`Adjacency`] in canonical form.
    ///
    /// The two breakends are ordered by [`Breakend::canonical_key`]. When
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
}
