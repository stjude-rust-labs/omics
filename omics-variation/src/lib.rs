//! Sequence, copy-number, and structural variants.
//!
//! `omics-variation` represents genomic changes that can be written as a
//! reference allele and an alternate allele at one locus. This covers single
//! nucleotide variants (`SNV`s), multi-nucleotide variants (`MNV`s),
//! insertions, deletions, and deletion-insertions (`delins`). Larger events
//! such as inversions, translocations, and breakends need additional structure
//! and are now modeled by the [`structural`] module.
//!
//! The central type is [`Variant`]. A [`Variant`] is generic over the
//! nucleotide type, so the same representation works for DNA and RNA alleles
//! from `omics-molecule`.
//!
//! Beyond small variants, the [`structural`] module models breakends, novel
//! adjacencies, and structural-variant events such as large deletions, large
//! insertions, inversions, tandem duplications, and translocations. A
//! structural variant derives its [`structural::Kind`] from breakend geometry
//! rather than storing it.
//!
//! The [`copy_number`] module models strandless, half-open copy-number
//! observations with typed absolute counts and required reference ploidy as
//! part of variant identity. Top-level aliases such as [`CopyNumberVariant`]
//! keep these types available alongside [`Variant`] and
//! [`StructuralVariant`]. Unlike [`structural`], a copy-number variant records
//! an observed count over an interval, carries the reference ploidy that
//! interprets that count, and does not encode a breakpoint adjacency.
//!
//! ```
//! use omics_molecule::polymer::dna;
//! use omics_variation::Variant;
//! use omics_variation::variant::Kind;
//!
//! let variant = "seq0:+:100(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
//! assert_eq!(variant.kind(), Kind::Snv);
//! assert_eq!(variant.to_string(), "seq0:+:100(b):A:C");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ```
//! use omics_molecule::polymer::dna;
//! use omics_variation::StructuralKind;
//! use omics_variation::StructuralVariant;
//! use omics_variation::structural::adjacency::Adjacency;
//! use omics_variation::structural::breakend::Breakend;
//! use omics_variation::structural::orientation::Orientation;
//!
//! let a = Breakend::try_new("seq0", Orientation::LowerFlank, 100)?;
//! let b = Breakend::try_new("seq0", Orientation::HigherFlank, 200)?;
//! let adjacency = Adjacency::<dna::Nucleotide>::try_new_paired(a, b, ".".parse()?)?;
//! let variant = StructuralVariant::try_new(vec![adjacency])?;
//! assert_eq!(variant.kind(), StructuralKind::Deletion);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ```
//! use omics_variation::CopyNumberChange;
//! use omics_variation::CopyNumberPloidy;
//! use omics_variation::CopyNumberVariant;
//!
//! let reference =
//!     CopyNumberVariant::try_new("seq0", 100, 200, 2, CopyNumberPloidy::DIPLOID)?;
//! assert_eq!(reference.count().get(), reference.ploidy().get());
//! assert_eq!(reference.change(), CopyNumberChange::Reference);
//!
//! let variant =
//!     CopyNumberVariant::try_new("seq0", 100, 200, 3, CopyNumberPloidy::DIPLOID)?;
//! assert_eq!(variant.change(), CopyNumberChange::Gain);
//!
//! let log2 = variant.log2();
//! let log10 = variant.log10();
//!
//! assert_eq!(
//!     CopyNumberVariant::try_from_log2("seq0", 100, 200, log2, CopyNumberPloidy::DIPLOID)?,
//!     variant
//! );
//! assert_eq!(
//!     CopyNumberVariant::try_from_log10("seq0", 100, 200, log10, CopyNumberPloidy::DIPLOID)?,
//!     variant
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ```
//! use omics_variation::CopyNumberVariant;
//!
//! let variant = "seq0:100-200(i):3/2".parse::<CopyNumberVariant>()?;
//! assert_eq!(variant.count().get(), 3);
//! assert_eq!(variant.to_string(), "seq0:100-200(i):3/2");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ```
//! use omics_variation::CopyNumberCount;
//! use omics_variation::CopyNumberPloidy;
//!
//! let count = CopyNumberCount::new(3);
//! let log2 = count.log2(CopyNumberPloidy::DIPLOID);
//!
//! assert_eq!(
//!     CopyNumberCount::try_from_log2(log2, CopyNumberPloidy::DIPLOID)?,
//!     count
//! );
//! # Ok::<(), omics_variation::copy_number::LogarithmicError>(())
//! ```
//!
//! # Variant strings
//!
//! The crate uses a compact, crate-local string format.
//!
//! ```text
//! contig:[strand:]position(system):reference:alternate
//! ```
//!
//! The `strand` field is optional when parsing and defaults to `+`.
//! Serialization always emits the explicit strand. The `position` field must
//! end in `(b)` for base coordinates or `(i)` for interbase coordinates. The
//! `reference` and `alternate` fields are nucleotide sequences; `.` denotes an
//! empty allele. This format is intentionally not [VCF], `SPDI`, or [HGVS].
//! Those formats carry additional metadata and normalization rules that belong
//! in dedicated adapters.
//!
//! | Kind | Example | Meaning |
//! | --- | --- | --- |
//! | [`Kind::Snv`] | `seq0:+:100(b):A:C` | substitute `A` with `C` at base `100` |
//! | [`Kind::Mnv`] | `seq0:+:100(b):AT:GC` | substitute `AT` with `GC` starting at base `100` |
//! | [`Kind::Insertion`] | `seq0:+:100(i):.:AT` | insert `AT` at interbase boundary `100` |
//! | [`Kind::Deletion`] | `seq0:+:100(b):AT:.` | delete `AT` starting at base `100` |
//! | [`Kind::Delins`] | `seq0:+:100(b):AT:G` | replace `AT` with `G` starting at base `100` |
//!
//! Parsing classifies a variant from its allele lengths. Identical
//! `reference` and `alternate` alleles are rejected because they do not
//! describe variation, and `.:.` is rejected because it changes nothing.
//! `SNV`, `MNV`, deletion, and `delins` records require `(b)` because their
//! position identifies an existing base. Insertion records require `(i)`
//! because their position identifies the boundary between bases. The opposite
//! qualifier is rejected for each kind.
//! Use [`Variant::reference_interval`] when comparing a variant to a reference
//! genome or asking which reference bases are consumed. Use
//! [`Variant::alternate_interval`] when asking about the local span of the
//! sequence introduced by the alternate allele. `SNV` and `MNV` variants use a
//! single interval because their reference and alternate alleles occupy the
//! same bases. The alternate interval is not a full query-genome coordinate;
//! downstream variants can shift global query coordinates.
//!
//! ```
//! use omics_molecule::polymer::dna;
//! use omics_variation::Variant;
//! use omics_variation::variant::Kind;
//!
//! let insertion = "seq0:+:100(i):.:AT".parse::<Variant<dna::Nucleotide>>()?;
//! assert_eq!(insertion.kind(), Kind::Insertion);
//!
//! let deletion = "seq0:100(b):AT:.".parse::<Variant<dna::Nucleotide>>()?;
//! assert_eq!(deletion.kind(), Kind::Deletion);
//! assert_eq!(deletion.to_string(), "seq0:+:100(b):AT:.");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Coordinates and spans
//!
//! `SNV`, `MNV`, deletion, and `delins` variants occupy bases, so their typed
//! variants use `Coordinate<Base>` and expose `Interval<Base>` spans.
//! Insertions occur between bases, so [`insertion::Variant`] uses
//! `Coordinate<Interbase>` and exposes a zero-width `Interval<Interbase>` with
//! [`insertion::Variant::interbase_interval`].
//!
//! ```
//! use omics_molecule::polymer::dna;
//! use omics_variation::Variant;
//!
//! let variant = "seq0:+:100(i):.:AT".parse::<Variant<dna::Nucleotide>>()?;
//! let Variant::Insertion(insertion) = variant else {
//!     // SAFETY: an empty reference allele always classifies as an insertion.
//!     unreachable!();
//! };
//!
//! let interval = insertion.interbase_interval();
//! assert_eq!(interval.start().position().get(), 100);
//! assert_eq!(interval.end().position().get(), 100);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Normalization and equality
//!
//! Variants preserve the alleles as written. [`PartialEq`] therefore compares
//! the literal stored form, not biological equivalence. Use
//! [`Variant::normalize`] when you need a parsimonious representation. The
//! normalization implemented here trims shared prefix and suffix bases and
//! moves the coordinate forward in the molecule's strand-aware direction. It
//! does not left-align through repeats because that requires a reference
//! sequence.
//!
//! ```
//! use omics_molecule::polymer::dna;
//! use omics_variation::Variant;
//!
//! let variant = "seq0:+:100(b):AT:AG".parse::<Variant<dna::Nucleotide>>()?;
//! assert_eq!(variant.normalize()?.to_string(), "seq0:+:101(b):T:G");
//!
//! let variant = "seq0:-:100(b):AT:AG".parse::<Variant<dna::Nucleotide>>()?;
//! assert_eq!(variant.normalize()?.to_string(), "seq0:-:99(b):T:G");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Typed variants
//!
//! Use the top-level [`Variant`] enum when parsing or dispatching across all
//! supported small variants. Use the modules under [`variant`] when a caller
//! already knows the expected kind and wants kind-specific accessors.
//!
//! ```
//! use omics_molecule::polymer::dna::Nucleotide;
//! use omics_variation::variant::deletion;
//!
//! let variant = deletion::Variant::<Nucleotide>::try_new("seq0:+:100", "AT")?;
//!
//! assert_eq!(variant.reference().to_string(), "AT");
//! assert_eq!(variant.interval().end().position().get(), 101);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! [`Kind`]: variant::Kind
//! [VCF]: https://samtools.github.io/hts-specs/VCFv4.5.pdf
//! [HGVS]: https://hgvs-nomenclature.org/stable/

use std::str::FromStr;

use omics_coordinate::Coordinate;
use omics_coordinate::Strand;
use omics_coordinate::position::Number;
use omics_coordinate::system::Base;
use omics_coordinate::system::Interbase;
use omics_core::VARIANT_SEPARATOR;
use omics_molecule::compound::Nucleotide;
use omics_molecule::sequence::Sequence;
use thiserror::Error;

pub mod copy_number;
pub mod small;
pub mod structural;

/// The top-level alias for a copy-number change.
pub use copy_number::Change as CopyNumberChange;
/// The top-level alias for a typed copy-number count.
pub use copy_number::Count as CopyNumberCount;
/// The top-level alias for a typed copy-number ploidy.
pub use copy_number::Ploidy as CopyNumberPloidy;
/// The top-level alias for a copy-number variant.
pub use copy_number::Variant as CopyNumberVariant;
/// The small-variant module, retained at its historical path.
pub use small as variant;
use small::Alteration;
use small::Kind;
pub use small::deletion;
pub use small::delins;
pub use small::insertion;
pub use small::mnv;
pub use small::snv;
pub use structural::Join;
pub use structural::Kind as StructuralKind;
pub use structural::Locality;
pub use structural::StructuralVariant;

/// An error related to a top-level [`Variant`].
#[derive(Error, Debug)]
pub enum Error {
    /// The variant string had the wrong number of `:`-separated parts.
    #[error("invalid variant format: `{0}`")]
    InvalidFormat(String),

    /// The coordinate portion failed to parse.
    #[error(transparent)]
    Coordinate(#[from] omics_coordinate::coordinate::Error),

    /// The coordinate system qualifier was missing or invalid.
    #[error(
        "position `{position}` must end with `(b)` for a base coordinate or `(i)` for an \
         interbase coordinate"
    )]
    CoordinateSystemQualifier {
        /// The unqualified or incorrectly qualified position token.
        position: String,
    },

    /// The coordinate system qualifier does not match the variant kind.
    #[error("wrong coordinate system for {kind}; expected {expected}; found {found}")]
    CoordinateSystem {
        /// The variant kind.
        kind: &'static str,

        /// The expected coordinate system qualifier.
        expected: &'static str,

        /// The coordinate system qualifier found in the variant string.
        found: &'static str,
    },

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

/// An interval associated with a variant allele.
///
/// Use the enum arms to preserve whether an allele interval sits on bases or
/// between bases. Callers that inspect reference coverage usually want
/// [`Variant::reference_interval`]. Callers that inspect the replacement
/// sequence usually want [`Variant::alternate_interval`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantInterval {
    /// A base interval.
    Base(omics_coordinate::Interval<Base>),

    /// An interbase interval.
    Interbase(omics_coordinate::Interval<Interbase>),
}

/// A coordinate-system qualifier parsed from a serialized variant position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PositionQualifier {
    /// A base coordinate qualifier.
    Base,

    /// An interbase coordinate qualifier.
    Interbase,
}

impl PositionQualifier {
    /// Gets the required coordinate-system qualifier for a variant kind.
    pub(crate) fn for_kind(kind: Kind) -> Self {
        match kind {
            Kind::Insertion => Self::Interbase,
            Kind::Snv | Kind::Mnv | Kind::Deletion | Kind::Delins => Self::Base,
        }
    }

    /// Gets the serialized suffix for this qualifier.
    pub(crate) fn suffix(self) -> &'static str {
        match self {
            Self::Base => "(b)",
            Self::Interbase => "(i)",
        }
    }

    /// Gets a human-readable description for error messages.
    pub(crate) fn description(self) -> &'static str {
        match self {
            Self::Base => "`(b)`",
            Self::Interbase => "`(i)`",
        }
    }
}

/// Gets a human-readable variant kind name for error messages.
fn kind_name(kind: Kind) -> &'static str {
    match kind {
        Kind::Snv => "snv",
        Kind::Mnv => "mnv",
        Kind::Insertion => "insertion",
        Kind::Deletion => "deletion",
        Kind::Delins => "delins",
    }
}

/// Splits a serialized position token into its numeric value and qualifier.
pub(crate) fn split_qualified_position(position: &str) -> Option<(&str, PositionQualifier)> {
    position
        .strip_suffix(PositionQualifier::Base.suffix())
        .map(|position| (position, PositionQualifier::Base))
        .or_else(|| {
            position
                .strip_suffix(PositionQualifier::Interbase.suffix())
                .map(|position| (position, PositionQualifier::Interbase))
        })
}

/// Formats a numeric position with its coordinate-system qualifier.
pub(crate) fn qualified_position(position: Number, qualifier: PositionQualifier) -> String {
    format!("{position}{}", qualifier.suffix())
}

/// A genomic variant.
///
/// `PartialEq` compares the stored (as-given) form literally: two spellings of
/// the same variant (for example `seq0:+:100(b):AT:AG` and `seq0:+:101(b):T:G`)
/// are **not** equal unless they are first put through
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
    /// let variant = "seq0:+:100(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
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

    /// Gets the interval spanned by the reference allele.
    ///
    /// Insertions return a zero-width interbase interval because their
    /// reference allele is empty. All other variant kinds return a base
    /// interval because they consume one or more reference bases.
    ///
    /// Use this method for reference-facing questions such as overlap against a
    /// reference annotation, determining which reference bases are consumed, or
    /// deciding whether a variant intersects a reference interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::Variant;
    /// use omics_variation::VariantInterval;
    ///
    /// let variant = "seq0:+:100(i):.:AT".parse::<Variant<dna::Nucleotide>>()?;
    /// assert!(matches!(
    ///     variant.reference_interval(),
    ///     VariantInterval::Interbase(_)
    /// ));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_interval(&self) -> VariantInterval {
        match self {
            Variant::Snv(variant) => VariantInterval::Base(variant.interval()),
            Variant::Mnv(variant) => VariantInterval::Base(variant.interval()),
            Variant::Insertion(variant) => VariantInterval::Interbase(variant.reference_interval()),
            Variant::Deletion(variant) => VariantInterval::Base(variant.reference_interval()),
            Variant::Delins(variant) => VariantInterval::Base(variant.reference_interval()),
        }
    }

    /// Gets the interval spanned by the alternate allele.
    ///
    /// Deletions return a zero-width interbase interval because their
    /// alternate allele is empty. Insertions and delins return local alternate
    /// base intervals when that span can be represented.
    ///
    /// Use this method for allele-facing questions such as the local span of an
    /// inserted or replacement sequence. Do not use it as a query-genome
    /// coordinate after applying multiple variants; this crate does not track
    /// cumulative coordinate shifts.
    ///
    /// # Examples
    ///
    /// ```
    /// use omics_molecule::polymer::dna;
    /// use omics_variation::Variant;
    /// use omics_variation::VariantInterval;
    ///
    /// let variant = "seq0:+:100(b):AT:.".parse::<Variant<dna::Nucleotide>>()?;
    /// assert!(matches!(
    ///     variant.alternate_interval(),
    ///     Some(VariantInterval::Interbase(_))
    /// ));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate_interval(&self) -> Option<VariantInterval> {
        match self {
            Variant::Snv(variant) => Some(VariantInterval::Base(variant.interval())),
            Variant::Mnv(variant) => Some(VariantInterval::Base(variant.interval())),
            Variant::Insertion(variant) => variant.alternate_interval().map(VariantInterval::Base),
            Variant::Deletion(variant) => {
                Some(VariantInterval::Interbase(variant.alternate_interval()))
            }
            Variant::Delins(variant) => variant.alternate_interval().map(VariantInterval::Base),
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
    /// let variant = "seq0:+:100(b):AT:AG".parse::<Variant<dna::Nucleotide>>()?;
    /// assert_eq!(variant.normalize()?.to_string(), "seq0:+:101(b):T:G");
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
                Variant::Insertion(insertion::Variant::try_from((boundary, trimmed))?)
            }
            Kind::Snv => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Snv(snv::Variant::try_from((c, trimmed))?)
            }
            Kind::Mnv => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Mnv(mnv::Variant::try_from((c, trimmed))?)
            }
            Kind::Deletion => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Deletion(deletion::Variant::try_from((c, trimmed))?)
            }
            Kind::Delins => {
                let c = start
                    .into_move_forward(prefix)
                    .ok_or(Error::NormalizeOverflow)?;
                Variant::Delins(delins::Variant::try_from((c, trimmed))?)
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

        let (position, found_qualifier) =
            split_qualified_position(position).ok_or_else(|| Error::CoordinateSystemQualifier {
                position: position.to_owned(),
            })?;

        let reference = reference
            .parse::<Sequence<N>>()
            .map_err(Error::ReferenceSequence)?;
        let alternate = alternate
            .parse::<Sequence<N>>()
            .map_err(Error::AlternateSequence)?;
        let alteration = Alteration::try_new(reference, alternate)?;
        let kind = alteration.kind();
        let expected_qualifier = PositionQualifier::for_kind(kind);

        if found_qualifier != expected_qualifier {
            return Err(Error::CoordinateSystem {
                kind: kind_name(kind),
                expected: expected_qualifier.description(),
                found: found_qualifier.description(),
            });
        }

        let coord = [contig, strand.as_str(), position].join(VARIANT_SEPARATOR);

        Ok(match kind {
            Kind::Insertion => {
                let coordinate = coord.parse::<Coordinate<Interbase>>()?;
                Variant::Insertion(insertion::Variant::try_from((coordinate, alteration))?)
            }
            Kind::Snv => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Snv(snv::Variant::try_from((coordinate, alteration))?)
            }
            Kind::Mnv => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Mnv(mnv::Variant::try_from((coordinate, alteration))?)
            }
            Kind::Deletion => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Deletion(deletion::Variant::try_from((coordinate, alteration))?)
            }
            Kind::Delins => {
                let coordinate = coord.parse::<Coordinate<Base>>()?;
                Variant::Delins(delins::Variant::try_from((coordinate, alteration))?)
            }
        })
    }
}

impl<N: Nucleotide> std::fmt::Display for Variant<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // The coordinate types differ (insertion is interbase), so read the
        // display primitives per arm rather than unifying the coordinate type.
        let (contig, strand, position, qualifier, alteration) = match self {
            Variant::Snv(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                PositionQualifier::Base,
                v.alteration(),
            ),
            Variant::Mnv(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                PositionQualifier::Base,
                v.alteration(),
            ),
            Variant::Deletion(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                PositionQualifier::Base,
                v.alteration(),
            ),
            Variant::Delins(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                PositionQualifier::Base,
                v.alteration(),
            ),
            Variant::Insertion(v) => (
                v.coordinate().contig(),
                v.coordinate().strand(),
                v.coordinate().position().get(),
                PositionQualifier::Interbase,
                v.alteration(),
            ),
        };
        let position = qualified_position(position, qualifier);

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
            "seq0:+:100(b):A:C"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Snv
        );
        assert_eq!(
            "seq0:+:100(b):AT:GC"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Mnv
        );
        assert_eq!(
            "seq0:+:100(i):.:AT"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Insertion
        );
        assert_eq!(
            "seq0:+:100(b):AT:."
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Deletion
        );
        assert_eq!(
            "seq0:+:100(b):AT:G"
                .parse::<Variant<dna::Nucleotide>>()?
                .kind(),
            Kind::Delins
        );
        Ok(())
    }

    #[test]
    fn it_defaults_missing_strand_to_positive() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:100(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
        match variant {
            Variant::Snv(snv) => assert_eq!(snv.coordinate().strand(), Strand::Positive),
            _ => panic!("expected SNV"),
        }
        Ok(())
    }

    #[test]
    fn display_round_trips_with_explicit_strand() -> Result<(), Box<dyn std::error::Error>> {
        for input in [
            "seq0:+:100(b):A:C",
            "seq0:+:100(b):AT:GC",
            "seq0:+:100(i):.:AT",
            "seq0:+:100(b):AT:.",
            "seq0:+:100(b):AT:G",
            "seq0:-:100(b):A:C",
            "seq0:-:100(b):AT:GC",
            "seq0:-:100(i):.:AT",
            "seq0:-:100(b):AT:.",
            "seq0:-:100(b):AT:G",
        ] {
            let variant = input.parse::<Variant<dna::Nucleotide>>()?;
            assert_eq!(variant.to_string(), input);
        }
        // Missing strand canonicalizes to explicit `+`.
        let variant = "seq0:100(b):A:C".parse::<Variant<dna::Nucleotide>>()?;
        assert_eq!(variant.to_string(), "seq0:+:100(b):A:C");

        let variant = "seq0:100(i):.:AT".parse::<Variant<dna::Nucleotide>>()?;
        assert_eq!(variant.to_string(), "seq0:+:100(i):.:AT");
        Ok(())
    }

    #[test]
    fn it_rejects_missing_coordinate_system_qualifiers() {
        let err = "seq0:+:100:A:C"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(
            err,
            Error::CoordinateSystemQualifier { ref position } if position == "100"
        ));
        assert_eq!(
            err.to_string(),
            "position `100` must end with `(b)` for a base coordinate or `(i)` for an interbase \
             coordinate"
        );
    }

    #[test]
    fn it_rejects_coordinate_system_mismatches() {
        for input in [
            "seq0:+:100(i):A:C",
            "seq0:+:100(i):AT:GC",
            "seq0:+:100(b):.:AT",
            "seq0:+:100(i):AT:.",
            "seq0:+:100(i):AT:G",
        ] {
            assert!(matches!(
                input.parse::<Variant<dna::Nucleotide>>().unwrap_err(),
                Error::CoordinateSystem { .. }
            ));
        }

        let err = "seq0:+:100(i):A:C"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "wrong coordinate system for snv; expected `(b)`; found `(i)`"
        );
    }

    #[test]
    fn it_reports_reference_and_alternate_intervals() {
        // SAFETY: this canonical insertion parses successfully.
        let insertion = "seq0:+:100(i):.:AT"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap();
        assert!(matches!(
            insertion.reference_interval(),
            VariantInterval::Interbase(_)
        ));
        // SAFETY: insertions always expose an alternate interval.
        match insertion.alternate_interval().unwrap() {
            VariantInterval::Base(interval) => {
                assert_eq!(interval.start().position().get(), 101);
                assert_eq!(interval.end().position().get(), 102);
            }
            VariantInterval::Interbase(_) => panic!("expected base alternate interval"),
        }

        // SAFETY: this canonical deletion parses successfully.
        let deletion = "seq0:+:100(b):AT:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap();
        assert!(matches!(
            deletion.reference_interval(),
            VariantInterval::Base(_)
        ));
        // SAFETY: deletions always expose an alternate interval.
        match deletion.alternate_interval().unwrap() {
            VariantInterval::Interbase(interval) => {
                assert_eq!(interval.start().position().get(), 99);
                assert_eq!(interval.end().position().get(), 99);
            }
            VariantInterval::Base(_) => panic!("expected interbase alternate interval"),
        }
    }

    #[test]
    fn it_rejects_both_empty_alleles() {
        let err = "seq0:+:100(b):.:."
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, Error::Alteration(variant::Error::BothEmpty)));
    }

    #[test]
    fn it_rejects_wrong_part_counts() {
        assert!(matches!(
            "seq0:100(b):A"
                .parse::<Variant<dna::Nucleotide>>()
                .unwrap_err(),
            Error::InvalidFormat(_)
        ));
    }

    fn normalized(input: &str) -> String {
        // SAFETY: this helper only receives canonical small-variant strings.
        let variant = input.parse::<Variant<dna::Nucleotide>>().unwrap();
        // SAFETY: these test inputs normalize without overflowing coordinates.
        variant.normalize().unwrap().to_string()
    }

    #[test]
    fn it_normalizes_mnv_to_snv() {
        assert_eq!(normalized("seq0:+:100(b):AT:AG"), "seq0:+:101(b):T:G");
    }

    #[test]
    fn it_normalizes_on_the_negative_strand() {
        // Trimming the shared prefix `A` moves the coordinate in the strand's
        // forward direction, which is downward on the negative strand.
        assert_eq!(normalized("seq0:-:100(b):AT:AG"), "seq0:-:99(b):T:G");
    }

    #[test]
    fn it_collapses_to_an_insertion_at_the_correct_boundary() {
        // `A:AT` at base 100 inserts `T` after base 100 -> interbase 100.
        // SAFETY: this canonical delins parses successfully.
        let variant = "seq0:+:100(b):A:AT"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap();
        // SAFETY: this test case normalizes without overflowing coordinates.
        let normalized = variant.normalize().unwrap();
        assert_eq!(normalized.kind(), Kind::Insertion);
        assert_eq!(normalized.to_string(), "seq0:+:100(i):.:T");

        // `A:TA` at base 100 inserts `T` before base 100 -> interbase 99.
        assert_eq!(
            normalized_kind_and_string("seq0:+:100(b):A:TA"),
            (Kind::Insertion, "seq0:+:99(i):.:T".to_string())
        );
    }

    #[test]
    fn it_collapses_delins_to_a_deletion() {
        // SAFETY: this canonical delins parses successfully.
        let variant = "seq0:+:100(b):ATG:AG"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap();
        // SAFETY: this test case normalizes without overflowing coordinates.
        let normalized = variant.normalize().unwrap();
        assert_eq!(normalized.kind(), Kind::Deletion);
        assert_eq!(normalized.to_string(), "seq0:+:101(b):T:.");
    }

    #[test]
    fn it_passes_insertions_through_unchanged() {
        assert_eq!(normalized("seq0:+:100(i):.:AT"), "seq0:+:100(i):.:AT");
    }

    fn normalized_kind_and_string(input: &str) -> (Kind, String) {
        // SAFETY: this helper only receives canonical small-variant strings.
        let variant = input.parse::<Variant<dna::Nucleotide>>().unwrap();
        // SAFETY: these test inputs normalize without overflowing coordinates.
        let normalized = variant.normalize().unwrap();
        (normalized.kind(), normalized.to_string())
    }

    #[test]
    fn it_collapses_to_an_insertion_on_the_negative_strand() {
        // Trimming moves in the strand's forward direction (downward on `-`).
        assert_eq!(
            normalized_kind_and_string("seq0:-:100(b):A:AT"),
            (Kind::Insertion, "seq0:-:99(i):.:T".to_string())
        );
        assert_eq!(
            normalized_kind_and_string("seq0:-:100(b):A:TA"),
            (Kind::Insertion, "seq0:-:100(i):.:T".to_string())
        );
    }

    #[test]
    fn it_parses_a_direct_interbase_insertion_at_zero() -> Result<(), Box<dyn std::error::Error>> {
        // Interbase 0 (start of contig) is a valid insertion anchor.
        let variant = "seq0:+:0(i):.:AT".parse::<Variant<dna::Nucleotide>>()?;
        assert_eq!(variant.kind(), Kind::Insertion);
        assert_eq!(variant.to_string(), "seq0:+:0(i):.:AT");
        Ok(())
    }

    #[test]
    fn it_rejects_a_negative_strand_span_underflow() {
        // A two-base MNV at base 1 on the negative strand would need base 0.
        let err = "seq0:-:1(b):AT:GC"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();
        assert!(matches!(err, Error::Kind(variant::KindError::SpanOverflow)));
    }

    #[test]
    fn it_round_trips_every_kind_on_the_negative_strand() -> Result<(), Box<dyn std::error::Error>>
    {
        for input in [
            "seq0:-:100(b):A:C",
            "seq0:-:100(b):AT:GC",
            "seq0:-:100(i):.:AT",
            "seq0:-:100(b):AT:.",
            "seq0:-:100(b):AT:G",
        ] {
            let variant = input.parse::<Variant<dna::Nucleotide>>()?;
            assert_eq!(variant.to_string(), input);
        }
        Ok(())
    }
}
