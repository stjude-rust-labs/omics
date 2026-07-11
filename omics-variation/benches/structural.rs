//! Benchmarks for structural variants.
#![expect(
    missing_docs,
    reason = "criterion_group generates undocumented registration functions"
)]

use std::hint::black_box;

use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use omics_molecule::polymer::dna;
use omics_variation::StructuralVariant;
use omics_variation::structural::adjacency::Adjacency;
use omics_variation::structural::breakend::Breakend;
use omics_variation::structural::orientation::Orientation;

/// A DNA structural variant used by these benchmarks.
type Sv = StructuralVariant<dna::Nucleotide>;

/// Builds a breakend on the named contig.
fn bnd(contig: &str, orientation: Orientation, position: u32) -> Breakend {
    Breakend::try_new(contig, orientation, position).unwrap()
}

/// Builds a paired adjacency from two breakends and an insertion token.
fn paired(x: Breakend, y: Breakend, insertion: &str) -> Adjacency<dna::Nucleotide> {
    Adjacency::try_new_paired(x, y, insertion.parse().unwrap()).unwrap()
}

/// Builds a large-deletion fixture.
fn deletion() -> Sv {
    Sv::try_new(vec![paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 5000),
        ".",
    )])
    .unwrap()
}

/// Builds a large-insertion fixture.
fn insertion() -> Sv {
    Sv::try_new(vec![paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq0", Orientation::HigherFlank, 100),
        "ACGTACGTACGT",
    )])
    .unwrap()
}

/// Builds an interchromosomal-translocation fixture.
fn interchromosomal_translocation() -> Sv {
    Sv::try_new(vec![paired(
        bnd("seq0", Orientation::LowerFlank, 100),
        bnd("seq1", Orientation::HigherFlank, 900),
        ".",
    )])
    .unwrap()
}

/// Builds an intrachromosomal-translocation fixture.
fn intrachromosomal_translocation() -> Sv {
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
    Sv::try_new(vec![origin, target_left, target_right]).unwrap()
}

/// Builds an inversion fixture.
fn inversion() -> Sv {
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
    Sv::try_new(vec![lower, higher]).unwrap()
}

/// Builds an internal-tandem-duplication fixture.
fn tandem_duplication() -> Sv {
    Sv::try_new(vec![paired(
        bnd("seq0", Orientation::HigherFlank, 100),
        bnd("seq0", Orientation::LowerFlank, 160),
        "GAT",
    )])
    .unwrap()
}

/// The six required events, paired with a stable benchmark label.
fn fixtures() -> Vec<(&'static str, Sv)> {
    vec![
        ("deletion", deletion()),
        ("insertion", insertion()),
        (
            "interchromosomal_translocation",
            interchromosomal_translocation(),
        ),
        (
            "intrachromosomal_translocation",
            intrachromosomal_translocation(),
        ),
        ("inversion", inversion()),
        ("tandem_duplication", tandem_duplication()),
    ]
}

/// Parses a serialized fixture back into a structural variant.
fn parse(input: &str) -> Sv {
    match input.parse() {
        Ok(variant) => variant,
        Err(err) => panic!("benchmark fixture failed to parse; {err}"),
    }
}

/// Registers structural-variant parsing benchmarks over all six events.
fn parse_benches(c: &mut Criterion) {
    for (label, variant) in fixtures() {
        let rendered = variant.to_string();
        c.bench_function(&format!("structural::parse::{label}"), |b| {
            b.iter(|| parse(black_box(&rendered)))
        });
    }
}

/// Registers structural-variant classification benchmarks over all six events.
fn classify_benches(c: &mut Criterion) {
    for (label, variant) in fixtures() {
        c.bench_function(&format!("structural::classify::{label}"), |b| {
            b.iter(|| black_box(&variant).kind())
        });
    }
}

criterion_group!(benches, parse_benches, classify_benches);
criterion_main!(benches);
