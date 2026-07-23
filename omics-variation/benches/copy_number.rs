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
const VARIANT: &str = "seq0:100-200(i):3/2";

/// The observed copy count used by the logarithmic benchmarks.
const COUNT: CopyNumberCount = CopyNumberCount::new(3);

/// The reference ploidy used by the variant benchmarks.
const PLOIDY: CopyNumberPloidy = CopyNumberPloidy::DIPLOID;

/// Constructs the canonical copy-number fixture.
fn construct_variant() -> CopyNumberVariant {
    match CopyNumberVariant::try_new(
        black_box("seq0"),
        black_box(100),
        black_box(200),
        black_box(3),
        black_box(PLOIDY),
    ) {
        Ok(variant) => variant,
        Err(err) => panic!("benchmark fixture failed to construct; {err}"),
    }
}

/// Constructs the canonical copy-number fixture from a base-2 logarithmic
/// ratio.
fn construct_variant_from_log2(
    value: f64,
) -> Result<CopyNumberVariant, omics_variation::copy_number::LogarithmicVariantError> {
    CopyNumberVariant::try_from_log2(
        black_box("seq0"),
        black_box(100),
        black_box(200),
        black_box(value),
        black_box(PLOIDY),
    )
}

/// Constructs the canonical copy-number fixture from a base-10 logarithmic
/// ratio.
fn construct_variant_from_log10(
    value: f64,
) -> Result<CopyNumberVariant, omics_variation::copy_number::LogarithmicVariantError> {
    CopyNumberVariant::try_from_log10(
        black_box("seq0"),
        black_box(100),
        black_box(200),
        black_box(value),
        black_box(PLOIDY),
    )
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

/// Classifies a copy-number variant against its stored ploidy.
fn classify_change(variant: &CopyNumberVariant) -> omics_variation::CopyNumberChange {
    black_box(variant).change()
}

/// Gets the base-2 logarithmic ratio from a copy-number variant.
fn variant_log2(variant: &CopyNumberVariant) -> f64 {
    black_box(variant).log2()
}

/// Gets the base-10 logarithmic ratio from a copy-number variant.
fn variant_log10(variant: &CopyNumberVariant) -> f64 {
    black_box(variant).log10()
}

/// Converts an absolute count into a base-2 logarithmic ratio.
fn to_log2(count: CopyNumberCount) -> f64 {
    black_box(count).log2(PLOIDY)
}

/// Converts an absolute count into a base-10 logarithmic ratio.
fn to_log10(count: CopyNumberCount) -> f64 {
    black_box(count).log10(PLOIDY)
}

/// Converts a base-2 logarithmic ratio into an absolute count.
fn from_log2(
    value: f64,
) -> Result<CopyNumberCount, omics_variation::copy_number::LogarithmicError> {
    CopyNumberCount::try_from_log2(black_box(value), PLOIDY)
}

/// Converts a base-10 logarithmic ratio into an absolute count.
fn from_log10(
    value: f64,
) -> Result<CopyNumberCount, omics_variation::copy_number::LogarithmicError> {
    CopyNumberCount::try_from_log10(black_box(value), PLOIDY)
}

/// Registers copy-number construction benchmarks.
fn construct_benches(c: &mut Criterion) {
    let log2 = to_log2(COUNT);
    let log10 = to_log10(COUNT);

    c.bench_function("copy_number::construct::variant", |b| {
        b.iter(construct_variant)
    });
    c.bench_function("copy_number::construct::variant_from_log2", |b| {
        b.iter(|| construct_variant_from_log2(black_box(log2)))
    });
    c.bench_function("copy_number::construct::variant_from_log10", |b| {
        b.iter(|| construct_variant_from_log10(black_box(log10)))
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
    c.bench_function("copy_number::variant::log2", |b| {
        b.iter(|| variant_log2(&variant))
    });
    c.bench_function("copy_number::variant::log10", |b| {
        b.iter(|| variant_log10(&variant))
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
