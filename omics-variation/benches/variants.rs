//! Benchmarks for small variants.
#![expect(
    missing_docs,
    reason = "criterion_group generates undocumented registration functions"
)]

use std::hint::black_box;

use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use omics_molecule::polymer::dna;
use omics_variation::Variant;
use omics_variation::VariantInterval;

/// A DNA variant used by these benchmarks.
type DnaVariant = Variant<dna::Nucleotide>;

/// A base-qualified single nucleotide variant fixture.
const SNV: &str = "seq0:+:100(b):A:C";

/// A base-qualified multi-nucleotide variant fixture.
const MNV: &str = "seq0:+:100(b):AT:GC";

/// An interbase-qualified insertion fixture.
const INSERTION: &str = "seq0:+:100(i):.:AT";

/// A base-qualified deletion fixture.
const DELETION: &str = "seq0:+:100(b):AT:.";

/// A base-qualified deletion-insertion fixture.
const DELINS: &str = "seq0:+:100(b):AT:G";

/// A fixture that normalizes to a single nucleotide variant.
const NORMALIZE_TO_SNV: &str = "seq0:+:100(b):AT:AG";

/// A fixture that normalizes to an insertion.
const NORMALIZE_TO_INSERTION: &str = "seq0:+:100(b):G:GAT";

/// A fixture missing the required coordinate-system qualifier.
const MISSING_QUALIFIER: &str = "seq0:+:100:A:C";

/// A fixture with a coordinate-system qualifier that does not match its kind.
const MISMATCHED_QUALIFIER: &str = "seq0:+:100(i):A:C";

/// Parses a benchmark fixture into a DNA variant.
fn parse_variant(input: &str) -> DnaVariant {
    match input.parse() {
        Ok(variant) => variant,
        Err(err) => panic!("benchmark fixture failed to parse; {err}"),
    }
}

/// Parses a benchmark fixture and reports whether it was rejected.
fn rejects_variant(input: &str) -> bool {
    black_box(input).parse::<DnaVariant>().is_err()
}

/// Parses the SNV fixture.
fn parse_snv() -> DnaVariant {
    parse_variant(black_box(SNV))
}

/// Parses the MNV fixture.
fn parse_mnv() -> DnaVariant {
    parse_variant(black_box(MNV))
}

/// Parses the insertion fixture.
fn parse_insertion() -> DnaVariant {
    parse_variant(black_box(INSERTION))
}

/// Parses the deletion fixture.
fn parse_deletion() -> DnaVariant {
    parse_variant(black_box(DELETION))
}

/// Parses the delins fixture.
fn parse_delins() -> DnaVariant {
    parse_variant(black_box(DELINS))
}

/// Displays a variant.
fn display_variant(variant: &DnaVariant) -> String {
    black_box(variant).to_string()
}

/// Normalizes a variant.
fn normalize_variant(variant: &DnaVariant) -> Result<DnaVariant, omics_variation::Error> {
    black_box(variant).normalize()
}

/// Gets a variant reference interval.
fn reference_interval(variant: &DnaVariant) -> VariantInterval {
    black_box(variant).reference_interval()
}

/// Gets a variant alternate interval.
fn alternate_interval(variant: &DnaVariant) -> Option<VariantInterval> {
    black_box(variant).alternate_interval()
}

/// Registers variant parsing benchmarks.
fn parse_benches(c: &mut Criterion) {
    c.bench_function("variants::parse::snv_base_qualified", |b| b.iter(parse_snv));
    c.bench_function("variants::parse::mnv_base_qualified", |b| b.iter(parse_mnv));
    c.bench_function("variants::parse::insertion_interbase_qualified", |b| {
        b.iter(parse_insertion)
    });
    c.bench_function("variants::parse::deletion_base_qualified", |b| {
        b.iter(parse_deletion)
    });
    c.bench_function("variants::parse::delins_base_qualified", |b| {
        b.iter(parse_delins)
    });
    c.bench_function("variants::parse::reject_missing_qualifier", |b| {
        b.iter(|| rejects_variant(MISSING_QUALIFIER))
    });
    c.bench_function("variants::parse::reject_mismatched_qualifier", |b| {
        b.iter(|| rejects_variant(MISMATCHED_QUALIFIER))
    });
}

/// Registers variant display benchmarks.
fn display_benches(c: &mut Criterion) {
    let snv = parse_variant(SNV);
    let insertion = parse_variant(INSERTION);

    c.bench_function("variants::display::snv_base_qualified", |b| {
        b.iter(|| display_variant(&snv))
    });
    c.bench_function("variants::display::insertion_interbase_qualified", |b| {
        b.iter(|| display_variant(&insertion))
    });
}

/// Registers variant normalization benchmarks.
fn normalize_benches(c: &mut Criterion) {
    let to_snv = parse_variant(NORMALIZE_TO_SNV);
    let to_insertion = parse_variant(NORMALIZE_TO_INSERTION);

    c.bench_function("variants::normalize::trim_to_snv", |b| {
        b.iter(|| normalize_variant(&to_snv))
    });
    c.bench_function("variants::normalize::collapse_to_insertion", |b| {
        b.iter(|| normalize_variant(&to_insertion))
    });
}

/// Registers variant interval benchmarks.
fn interval_benches(c: &mut Criterion) {
    let deletion = parse_variant(DELETION);
    let delins = parse_variant(DELINS);
    let insertion = parse_variant(INSERTION);

    c.bench_function("variants::intervals::reference_deletion", |b| {
        b.iter(|| reference_interval(&deletion))
    });
    c.bench_function("variants::intervals::reference_insertion", |b| {
        b.iter(|| reference_interval(&insertion))
    });
    c.bench_function("variants::intervals::alternate_deletion", |b| {
        b.iter(|| alternate_interval(&deletion))
    });
    c.bench_function("variants::intervals::alternate_delins", |b| {
        b.iter(|| alternate_interval(&delins))
    });
    c.bench_function("variants::intervals::alternate_insertion", |b| {
        b.iter(|| alternate_interval(&insertion))
    });
}

criterion_group!(
    benches,
    parse_benches,
    display_benches,
    normalize_benches,
    interval_benches,
);
criterion_main!(benches);
