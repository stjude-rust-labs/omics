//! Acceptance tests for the six required structural-variant events.
//!
//! Each event is constructed through the real, opaque `Adjacency` API and its
//! derived [`Kind`] is asserted. The swap-invariance tests additionally build
//! every single-adjacency event with its breakends supplied in both orders,
//! confirming that canonicalization makes the two adjacencies equal and that
//! their enclosing structural variants classify identically.

use omics_coordinate::position::Number;
use omics_molecule::polymer::dna;
use omics_variation::StructuralVariant;
use omics_variation::structural::Kind;
use omics_variation::structural::Locality;
use omics_variation::structural::adjacency::Adjacency;
use omics_variation::structural::breakend::Breakend;
use omics_variation::structural::orientation::Orientation;

/// A DNA structural variant used by these acceptance tests.
type Sv = StructuralVariant<dna::Nucleotide>;

/// Builds a breakend on the named contig.
fn bnd(contig: &str, orientation: Orientation, position: Number) -> Breakend {
    Breakend::try_new(contig, orientation, position).unwrap()
}

/// Builds a paired adjacency from two breakends and an insertion token.
fn paired(x: Breakend, y: Breakend, insertion: &str) -> Adjacency<dna::Nucleotide> {
    Adjacency::try_new_paired(x, y, insertion.parse().unwrap()).unwrap()
}

/// Asserts that a single-adjacency event built in both breakend orders is
/// order-invariant.
///
/// The junction is supplied once as `(x, y)` with `insertion_forward` and once
/// as `(y, x)` with `insertion_reversed`. Because canonicalization reverse
/// complements the insertion of whichever call it reorders, the two adjacencies
/// must be equal, and both enclosing structural variants must classify as
/// `expected`.
fn assert_swap_invariant(
    x: Breakend,
    y: Breakend,
    insertion_forward: &str,
    insertion_reversed: &str,
    expected: Kind,
) {
    let forward = Adjacency::<dna::Nucleotide>::try_new_paired(
        x.clone(),
        y.clone(),
        insertion_forward.parse().unwrap(),
    )
    .unwrap();
    let backward =
        Adjacency::<dna::Nucleotide>::try_new_paired(y, x, insertion_reversed.parse().unwrap())
            .unwrap();

    assert_eq!(forward, backward, "swapped-order adjacencies must be equal");

    let forward_variant = Sv::try_new(vec![forward]).unwrap();
    let backward_variant = Sv::try_new(vec![backward]).unwrap();

    assert_eq!(forward_variant, backward_variant);
    assert_eq!(forward_variant.kind(), expected);
    assert_eq!(backward_variant.kind(), expected);
}

#[test]
fn large_deletion() {
    let adjacency = paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 5000),
        ".",
    );
    assert_eq!(Sv::try_new(vec![adjacency]).unwrap().kind(), Kind::Deletion);
}

#[test]
fn large_insertion() {
    let adjacency = paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 100),
        "ACGTACGTACGT",
    );
    assert_eq!(
        Sv::try_new(vec![adjacency]).unwrap().kind(),
        Kind::Insertion
    );
}

#[test]
fn interchromosomal_translocation() {
    let adjacency = paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq1", Orientation::HigherFlank, 900),
        ".",
    );
    assert_eq!(
        Sv::try_new(vec![adjacency]).unwrap().kind(),
        Kind::Translocation(Locality::Interchromosomal)
    );
}

#[test]
fn intrachromosomal_translocation() {
    // A balanced forward relocation: three co-linear junctions over the three
    // boundaries {100, 200, 400}, each appearing once as a lower flank and once
    // as a higher flank across the junctions.
    let origin = paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 200),
        ".",
    );
    let target_left = paired(
        bnd("seq0", Orientation::LowerFlank, 400),
        bnd("seq0", Orientation::HigherFlank, 100),
        ".",
    );
    let target_right = paired(
        bnd("seq0", Orientation::LowerFlank, 200),
        bnd("seq0", Orientation::HigherFlank, 400),
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
fn inversion() {
    let lower = paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::LowerFlank, 300),
        ".",
    );
    let higher = paired(
        bnd("seq0", Orientation::HigherFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 300),
        ".",
    );
    assert_eq!(
        Sv::try_new(vec![lower, higher]).unwrap().kind(),
        Kind::Inversion
    );
}

#[test]
fn internal_tandem_duplication() {
    let adjacency = paired(
        bnd("seq0", Orientation::HigherFlank, 100),
        bnd("seq0", Orientation::LowerFlank, 160),
        "GAT",
    );
    assert_eq!(
        Sv::try_new(vec![adjacency]).unwrap().kind(),
        Kind::TandemDuplication
    );
}

#[test]
fn deletion_is_swap_invariant() {
    assert_swap_invariant(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 5000),
        ".",
        ".",
        Kind::Deletion,
    );
}

#[test]
fn insertion_is_swap_invariant() {
    // `AAAC` forward reverse-complements to `GTTT`, so the swapped call must be
    // supplied `GTTT` to land on the same canonical insertion.
    assert_swap_invariant(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 100),
        "AAAC",
        "GTTT",
        Kind::Insertion,
    );
}

#[test]
fn tandem_duplication_is_swap_invariant() {
    assert_swap_invariant(
        bnd("seq0", Orientation::HigherFlank, 100),
        bnd("seq0", Orientation::LowerFlank, 160),
        "AAAC",
        "GTTT",
        Kind::TandemDuplication,
    );
}

#[test]
fn interchromosomal_translocation_is_swap_invariant() {
    assert_swap_invariant(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq1", Orientation::HigherFlank, 900),
        ".",
        ".",
        Kind::Translocation(Locality::Interchromosomal),
    );
}

#[test]
fn reverse_complement_lands_in_the_canonical_frame() {
    // The same physical junction supplied in both orders must classify the same,
    // and the insertion must be reverse-complemented into the canonical frame.
    // `AAAC` reverse-complements to `GTTT`, so the two orders differ observably
    // in their stored insertion when the reversed call is supplied `AAAC`.
    let lo = bnd("seq0", Orientation::LowerFlank, 100);
    let hi = bnd("seq0", Orientation::HigherFlank, 200);

    let forward = Adjacency::<dna::Nucleotide>::try_new_paired(
        lo.clone(),
        hi.clone(),
        "AAAC".parse().unwrap(),
    )
    .unwrap();
    let backward =
        Adjacency::<dna::Nucleotide>::try_new_paired(hi, lo, "AAAC".parse().unwrap()).unwrap();

    let (_, _, forward_insertion) = forward.paired().expect("a paired adjacency");
    let (_, _, backward_insertion) = backward.paired().expect("a paired adjacency");
    assert_eq!(forward_insertion.to_string(), "AAAC");
    assert_eq!(backward_insertion.to_string(), "GTTT");
}
