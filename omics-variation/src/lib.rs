//! Genomic variation.

use std::str::FromStr;

use omics_coordinate::Coordinate;
use omics_coordinate::Strand;
use omics_coordinate::system::Base;
use omics_coordinate::system::Interbase;
use omics_core::VARIANT_SEPARATOR;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;
use thiserror::Error;

pub mod variant;

use variant::Alteration;
use variant::Kind;
pub use variant::deletion;
pub use variant::delins;
pub use variant::insertion;
pub use variant::mnv;
pub use variant::snv;

/// An error related to a top-level [`Variant`].
#[derive(Error, Debug)]
pub enum Error {
    /// The variant string had the wrong number of `:`-separated parts.
    #[error("invalid variant format: `{0}`")]
    InvalidFormat(String),

    /// The coordinate portion failed to parse.
    #[error(transparent)]
    Coordinate(#[from] omics_coordinate::coordinate::Error),

    /// The reference allele failed to parse.
    #[error("reference allele parse error: {0}")]
    ReferenceSequence(omics_molecule::sequence::ParseError),

    /// The alternate allele failed to parse.
    #[error("alternate allele parse error: {0}")]
    AlternateSequence(omics_molecule::sequence::ParseError),

    /// The alleles did not form a valid alteration (both-empty or identical).
    #[error(transparent)]
    Alteration(#[from] variant::Error),

    /// A typed kind could not be constructed from a classified alteration.
    #[error(transparent)]
    Kind(#[from] variant::KindError),

    /// A coordinate shift during normalization overflowed the position bounds.
    #[error("normalization overflowed the coordinate position")]
    NormalizeOverflow,
}

/// A genomic variant.
///
/// `PartialEq` compares the stored (as-given) form literally: two spellings of
/// the same variant (for example `seq0:+:100:AT:AG` and `seq0:+:101:T:G`) are
/// **not** equal unless they are first put through
/// [`normalize`](Self::normalize).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Variant<N: Nucleotide> {
    /// A single nucleotide variant.
    Snv(snv::Variant<N>),
    /// A multi-nucleotide variant.
    Mnv(mnv::Variant<N>),
    /// An insertion.
    Insertion(insertion::Variant<N>),
    /// A deletion.
    Deletion(deletion::Variant<N>),
    /// A combined deletion-insertion.
    Delins(delins::Variant<N>),
}

impl<N: Nucleotide> Variant<N> {
    /// Gets the [`Kind`] of this variant.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::Variant;
    /// use omics_variation::variant::Kind;
    ///
    /// let variant = "seq0:+:100:A:C".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.kind(), Kind::Snv);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn kind(&self) -> Kind {
        match self {
            Variant::Snv(_) => Kind::Snv,
            Variant::Mnv(_) => Kind::Mnv,
            Variant::Insertion(_) => Kind::Insertion,
            Variant::Deletion(_) => Kind::Deletion,
            Variant::Delins(_) => Kind::Delins,
        }
    }

    /// Returns a normalized copy: shared prefix/suffix bases trimmed and the
    /// coordinate adjusted. May reclassify (for example an MNV collapsing to an
    /// SNV, or a delins to an insertion or deletion).
    ///
    /// Trimming is parsimony only; left-alignment through repeats requires a
    /// reference sequence and is not performed. Insertions have an empty
    /// reference, so nothing can trim and they are returned unchanged.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::Variant;
    ///
    /// let variant = "seq0:+:100:AT:AG".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.normalize()?.to_string(), "seq0:+:101:T:G");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn normalize(&self) -> Result<Variant<N>, Error> {
        // Only the four base-coordinate kinds can trim; insertions cannot.
        let (start, alteration) = match self {
            Variant::Insertion(_) => return Ok(self.clone()),
            Variant::Snv(v) => (v.coordinate().clone(), v.alteration().clone()),
            Variant::Mnv(v) => (v.coordinate().clone(), v.alteration().clone()),
            Variant::Deletion(v) => (v.coordinate().clone(), v.alteration().clone()),
            Variant::Delins(v) => (v.coordinate().clone(), v.alteration().clone()),
        };

        let (prefix, trimmed) = alteration.trimmed();
        let prefix = omics_coordinate::position::Number::try_from(prefix)
            .map_err(|_| Error::NormalizeOverflow)?;

        Ok(match trimmed.kind() {
            Kind::Insertion => {
                // Collapse-to-insertion: the boundary is one interbase step back
                // from the base start, then forward by the trimmed prefix.
                let boundary = start
                    .nudge_backward()
                    .ok_or(Error::NormalizeOverflow)?
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Insertion(insertion::Variant::try_new(boundary, trimmed)?)
            }
            Kind::Snv => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Snv(snv::Variant::from_alteration(c, trimmed)?)
            }
            Kind::Mnv => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Mnv(mnv::Variant::try_new(c, trimmed)?)
            }
            Kind::Deletion => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Deletion(deletion::Variant::try_new(c, trimmed)?)
            }
            Kind::Delins => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Delins(delins::Variant::try_new(c, trimmed)?)
            }
        })
    }
}

impl<N: Nucleotide> FromStr for Variant<N> {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(VARIANT_SEPARATOR).collect::<Vec<_>>();
        let (contig, strand, position, reference, alternate) = match parts.as_slice() {
            [c, p, r, a] => (*c, Strand::Positive.to_string(), *p, *r, *a),
            [c, st, p, r, a] => (*c, (*st).to_string(), *p, *r, *a),
            _ => return Err(Error::InvalidFormat(s.to_owned())),
        };

        let reference = reference
            .parse::<Sequence<N>>()
            .map_err(Error::ReferenceSequence)?;
        let alternate = alternate
            .parse::<Sequence<N>>()
            .map_err(Error::AlternateSequence)?;
        let alteration = Alteration::try_new(reference, alternate)?;

        let coord = [contig, strand.as_str(), position].join(VARIANT_SEPARATOR);

        Ok(match alteration.kind() {
            Kind::Insertion => {
                let coordinate = coord.parse::<Coordinate<Interbase>>()?;
                Variant::Insertion(insertion::Variant::try_new(coordinate, alteration)?)
            }
            Kind::Snv => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Snv(snv::Variant::from_alteration(coordinate, alteration)?)
            }
            Kind::Mnv => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Mnv(mnv::Variant::try_new(coordinate, alteration)?)
            }
            Kind::Deletion => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Deletion(deletion::Variant::try_new(coordinate, alteration)?)
            }
            Kind::Delins => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Delins(delins::Variant::try_new(coordinate, alteration)?)
            }
        })
    }
}

impl<N: Nucleotide> std::fmt::Display for Variant<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // The coordinate types differ (insertion is interbase), so read the
        // display primitives per arm rather than unifying the coordinate type.
        let (contig, strand, position, alteration) = match self {
            Variant::Snv(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                v.alteration(),
            ),
            Variant::Mnv(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                v.alteration(),
            ),
            Variant::Deletion(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                v.alteration(),
            ),
            Variant::Delins(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                v.alteration(),
            ),
            Variant::Insertion(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                v.alteration(),
            ),
        };

        write!(
            f,
            "{contig}{sep}{strand}{sep}{position}{sep}{}{sep}{}",
            alteration.reference(),
            alteration.alternate(),
            sep = VARIANT_SEPARATOR,
        )
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::Strand;
    use omics_molecule::polymer::dna;

    use super::*;
    use crate::variant::Kind;

    #[test]
    fn it_parses_and_classifies_all_small_variants() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(
            "seq0:+:100:A:C".parse::<Variant<dna::Nucleotide>>()?.kind(),
            Kind::Snv
        );
        assert_eq!(
            "seq0:+:100:AT:GC"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Mnv
        );
        assert_eq!(
            "seq0:+:100:.:AT"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Insertion
        );
        assert_eq!(
            "seq0:+:100:AT:."
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Deletion
        );
        assert_eq!(
            "seq0:+:100:AT:G"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Delins
        );
        Ok(())
    }

    #[test]
    fn it_defaults_missing_strand_to_positive() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:100:A:C".parse::<Variant<dna::Nucleotide>>()?;
        match variant {
            Variant::Snv(snv) => assert_eq!(snv.coordinate().strand(), Strand::Positive),
            _ => panic!("expected SNV"),
        }
        Ok(())
    }

    #[test]
    fn display_round_trips_with_explicit_strand() -> Result<(), Box<dyn std::error::Error>> {
        for input in [
            "seq0:+:100:A:C",
            "seq0:+:100:AT:GC",
            "seq0:+:100:.:AT",
            "seq0:+:100:AT:.",
            "seq0:+:100:AT:G",
        ] {
            let variant = input.parse::<Variant<dna::Nucleotide>>()?;
            assert_eq!(variant.to_string(), input);
        }
        // Missing strand canonicalizes to explicit `+`.
        let variant = "seq0:100:A:C".parse::<Variant<dna::Nucleotide>>()?;
        assert_eq!(variant.to_string(), "seq0:+:100:A:C");
        Ok(())
    }

    #[test]
    fn it_rejects_both_empty_alleles() {
        let err = "seq0:+:100:.:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, Error::Alteration(variant::Error::BothEmpty)));
    }

    #[test]
    fn it_rejects_wrong_part_counts() {
        assert!(matches!(
            "seq0:100:A"
                .parse::<Variant<dna::Nucleotide>>()
                .unwrap_err(),
            Error::InvalidFormat(_)
        ));
    }

    fn normalized(input: &str) -> String {
        input
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap()
            .normalize()
            .unwrap()
            .to_string()
    }

    #[test]
    fn it_normalizes_mnv_to_snv() {
        assert_eq!(normalized("seq0:+:100:AT:AG"), "seq0:+:101:T:G");
    }

    #[test]
    fn it_normalizes_on_the_negative_strand() {
        // Trimming the shared prefix `A` moves the coordinate in the strand's
        // forward direction, which is downward on the negative strand.
        assert_eq!(normalized("seq0:-:100:AT:AG"), "seq0:-:99:T:G");
    }

    #[test]
    fn it_collapses_to_an_insertion_at_the_correct_boundary() {
        // `A:AT` at base 100 inserts `T` after base 100 -> interbase 100.
        let variant = "seq0:+:100:A:AT"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap();
        let normalized = variant.normalize().unwrap();
        assert_eq!(normalized.kind(), Kind::Insertion);
        assert_eq!(normalized.to_string(), "seq0:+:100:.:T");

        // `A:TA` at base 100 inserts `T` before base 100 -> interbase 99.
        assert_eq!(
            normalized_kind_and_string("seq0:+:100:A:TA"),
            (Kind::Insertion, "seq0:+:99:.:T".to_string())
        );
    }

    #[test]
    fn it_collapses_delins_to_a_deletion() {
        let variant = "seq0:+:100:ATG:AG"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap();
        let normalized = variant.normalize().unwrap();
        assert_eq!(normalized.kind(), Kind::Deletion);
        assert_eq!(normalized.to_string(), "seq0:+:101:T:.");
    }

    #[test]
    fn it_passes_insertions_through_unchanged() {
        assert_eq!(normalized("seq0:+:100:.:AT"), "seq0:+:100:.:AT");
    }

    fn normalized_kind_and_string(input: &str) -> (Kind, String) {
        let normalized = input
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap()
            .normalize()
            .unwrap();
        (normalized.kind(), normalized.to_string())
    }
}
