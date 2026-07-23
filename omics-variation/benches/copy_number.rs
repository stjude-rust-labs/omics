//! Benchmarks for copy-number variants.
#![expect(
    missing_docs,
    reason = "criterion_group generates undocumented registration functions"
)]

use std::hint::black_box;

use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use omics_variation::CopyNumberCount;
use omics_variation::CopyNumberPloidy;
use omics_variation::CopyNumberVariant;

/// A canonical copy-number variant fixture.
const VARIANT: &str = "seq0:100-200(i):3";

/// The observed copy count used by the logarithmic benchmarks.
const COUNT: CopyNumberCount = CopyNumberCount::new(3);

/// The baseline count used by the change benchmark.
const BASELINE: CopyNumberCount = CopyNumberCount::new(2);

/// Constructs the canonical copy-number fixture.
fn construct_variant() -> CopyNumberVariant {
    match CopyNumberVariant::try_new(
        black_box("seq0"),
        black_box(100),
        black_box(200),
        black_box(3),
    ) {
        Ok(variant) => variant,
        Err(err) => panic!("benchmark fixture failed to construct; {err}"),
    }
}

/// Parses the canonical copy-number fixture.
fn parse_variant(input: &str) -> CopyNumberVariant {
    match input.parse() {
        Ok(variant) => variant,
        Err(err) => panic!("benchmark fixture failed to parse; {err}"),
    }
}

/// Displays a copy-number variant.
fn display_variant(variant: &CopyNumberVariant) -> String {
    black_box(variant).to_string()
}

/// Classifies a copy-number variant against the baseline count.
fn classify_change(variant: &CopyNumberVariant) -> omics_variation::CopyNumberChange {
    black_box(variant).change(black_box(BASELINE))
}

/// Converts an absolute count into a base-2 logarithmic ratio.
fn to_log2(count: CopyNumberCount) -> f64 {
    black_box(count).log2(CopyNumberPloidy::DIPLOID)
}

/// Converts an absolute count into a base-10 logarithmic ratio.
fn to_log10(count: CopyNumberCount) -> f64 {
    black_box(count).log10(CopyNumberPloidy::DIPLOID)
}

/// Converts a base-2 logarithmic ratio into an absolute count.
fn from_log2(
    value: f64,
) -> Result<CopyNumberCount, omics_variation::copy_number::LogarithmicError> {
    CopyNumberCount::try_from_log2(black_box(value), CopyNumberPloidy::DIPLOID)
}

/// Converts a base-10 logarithmic ratio into an absolute count.
fn from_log10(
    value: f64,
) -> Result<CopyNumberCount, omics_variation::copy_number::LogarithmicError> {
    CopyNumberCount::try_from_log10(black_box(value), CopyNumberPloidy::DIPLOID)
}

/// Registers copy-number construction benchmarks.
fn construct_benches(c: &mut Criterion) {
    c.bench_function("copy_number::construct::variant", |b| {
        b.iter(construct_variant)
    });
}

/// Registers copy-number parsing benchmarks.
fn parse_benches(c: &mut Criterion) {
    c.bench_function("copy_number::parse::canonical", |b| {
        b.iter(|| parse_variant(black_box(VARIANT)))
    });
}

/// Registers copy-number display benchmarks.
fn display_benches(c: &mut Criterion) {
    let variant = parse_variant(VARIANT);

    c.bench_function("copy_number::display::canonical", |b| {
        b.iter(|| display_variant(&variant))
    });
}

/// Registers copy-number change benchmarks.
fn change_benches(c: &mut Criterion) {
    let variant = parse_variant(VARIANT);

    c.bench_function("copy_number::change::classify", |b| {
        b.iter(|| classify_change(&variant))
    });
}

/// Registers copy-number logarithmic conversion benchmarks.
fn logarithmic_benches(c: &mut Criterion) {
    let log2 = to_log2(COUNT);
    let log10 = to_log10(COUNT);

    c.bench_function("copy_number::count::to_log2", |b| b.iter(|| to_log2(COUNT)));
    c.bench_function("copy_number::count::to_log10", |b| {
        b.iter(|| to_log10(COUNT))
    });
    c.bench_function("copy_number::count::from_log2", |b| {
        b.iter(|| from_log2(black_box(log2)))
    });
    c.bench_function("copy_number::count::from_log10", |b| {
        b.iter(|| from_log10(black_box(log10)))
    });
}

criterion_group!(
    benches,
    construct_benches,
    parse_benches,
    display_benches,
    change_benches,
    logarithmic_benches,
);
criterion_main!(benches);
