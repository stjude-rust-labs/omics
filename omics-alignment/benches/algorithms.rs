//! Benchmarks for global and local alignment algorithms.
#![expect(
    missing_docs,
    reason = "criterion_group generates undocumented registration functions"
)]

use std::hint::black_box;

use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::Throughput;
use criterion::criterion_group;
use criterion::criterion_main;
use omics_alignment::algorithm::Scoring;
use omics_alignment::algorithm::global;
use omics_alignment::algorithm::local;
use omics_alignment::algorithm::simd;

/// Equal input lengths used to measure quadratic scaling.
const SIZES: [usize; 4] = [16, 64, 256, 1024];

/// Input shapes that exercise common alignment workloads.
#[derive(Clone, Copy)]
enum Fixture {
    /// Equal-length inputs with regular substitutions.
    Balanced,
    /// Inputs with a shorter query sequence.
    Asymmetric,
    /// Inputs with repeated query-only regions.
    GapHeavy,
}

/// Input shapes registered for every benchmark size and algorithm.
const FIXTURES: [Fixture; 3] = [Fixture::Balanced, Fixture::Asymmetric, Fixture::GapHeavy];

impl Fixture {
    /// Returns the stable Criterion label for this input shape.
    const fn label(self) -> &'static str {
        match self {
            Self::Balanced => "balanced",
            Self::Asymmetric => "asymmetric",
            Self::GapHeavy => "gap-heavy",
        }
    }

    /// Builds benchmark inputs with the requested base reference length.
    fn inputs(self, size: usize) -> (Vec<u8>, Vec<u8>) {
        match self {
            Self::Balanced => balanced_fixture(size),
            Self::Asymmetric => asymmetric_fixture(size),
            Self::GapHeavy => gap_heavy_fixture(size),
        }
    }
}

/// Builds a repeating DNA-like sequence with the requested length.
fn patterned_bytes(size: usize) -> Vec<u8> {
    (0..size)
        .map(|index| b"ACGT"[index % 4])
        .collect::<Vec<_>>()
}

/// Builds equal-length, mostly matching DNA-like benchmark inputs.
fn balanced_fixture(size: usize) -> (Vec<u8>, Vec<u8>) {
    let reference = patterned_bytes(size);
    let mut query = reference.clone();

    for chunk in query.chunks_mut(16) {
        if chunk.len() > 3 {
            chunk[3] = b'A';
        }
        if chunk.len() > 8 {
            let rotation_end = chunk.len().min(12);
            chunk[8..rotation_end].rotate_left(1);
        }
    }

    (reference, query)
}

/// Builds inputs whose query is half the reference length.
fn asymmetric_fixture(size: usize) -> (Vec<u8>, Vec<u8>) {
    let reference = patterned_bytes(size);
    let mut query = patterned_bytes(size / 2);

    for index in (5..query.len()).step_by(19) {
        query[index] = b'N';
    }

    (reference, query)
}

/// Builds inputs containing a four-byte query-only region per reference block.
fn gap_heavy_fixture(size: usize) -> (Vec<u8>, Vec<u8>) {
    let reference = patterned_bytes(size);
    let mut query = Vec::with_capacity(size + 4 * reference.len().div_ceil(24));

    for block in reference.chunks(24) {
        let midpoint = block.len() / 2;
        query.extend_from_slice(&block[..midpoint]);
        query.extend_from_slice(b"NNNN");
        query.extend_from_slice(&block[midpoint..]);
    }

    (reference, query)
}

/// Constructs the scoring configuration shared by every benchmark.
fn scoring() -> Scoring {
    // SAFETY: the benchmark scores satisfy all Scoring sign constraints.
    Scoring::try_new(2, -3, -2, -1).unwrap()
}

/// Returns the number of dynamic-programming cells for two input lengths.
fn cells(reference_length: usize, query_length: usize) -> u64 {
    ((reference_length + 1) * (query_length + 1)) as u64
}

/// Registers end-to-end global-alignment scaling benchmarks.
fn global_benches(c: &mut Criterion) {
    let scoring = scoring();

    for fixture in FIXTURES {
        let mut group = c.benchmark_group(format!("algorithms/global/{}", fixture.label()));

        for size in SIZES {
            let (reference, query) = fixture.inputs(size);
            group.throughput(Throughput::Elements(cells(reference.len(), query.len())));
            group.bench_with_input(BenchmarkId::new("scalar", size), &size, |b, _| {
                b.iter(|| {
                    match global(
                        black_box(reference.as_slice()),
                        black_box(query.as_slice()),
                        black_box(scoring),
                    ) {
                        Ok(outcome) => black_box(outcome),
                        Err(error) => panic!("scalar global benchmark failed; {error}"),
                    }
                })
            });
            group.bench_with_input(BenchmarkId::new("simd", size), &size, |b, _| {
                b.iter(|| {
                    match simd::global(
                        black_box(reference.as_slice()),
                        black_box(query.as_slice()),
                        black_box(scoring),
                    ) {
                        Ok(outcome) => black_box(outcome),
                        Err(error) => panic!("SIMD global benchmark failed; {error}"),
                    }
                })
            });
        }

        group.finish();
    }
}

/// Registers end-to-end local-alignment scaling benchmarks.
fn local_benches(c: &mut Criterion) {
    let scoring = scoring();

    for fixture in FIXTURES {
        let mut group = c.benchmark_group(format!("algorithms/local/{}", fixture.label()));

        for size in SIZES {
            let (reference, query) = fixture.inputs(size);
            group.throughput(Throughput::Elements(cells(reference.len(), query.len())));
            group.bench_with_input(BenchmarkId::new("scalar", size), &size, |b, _| {
                b.iter(|| {
                    match local(
                        black_box(reference.as_slice()),
                        black_box(query.as_slice()),
                        black_box(scoring),
                    ) {
                        Ok(Some(outcome)) => black_box(outcome),
                        Ok(None) => panic!("scalar local benchmark produced no positive alignment"),
                        Err(error) => panic!("scalar local benchmark failed; {error}"),
                    }
                })
            });
            group.bench_with_input(BenchmarkId::new("simd", size), &size, |b, _| {
                b.iter(|| {
                    match simd::local(
                        black_box(reference.as_slice()),
                        black_box(query.as_slice()),
                        black_box(scoring),
                    ) {
                        Ok(Some(outcome)) => black_box(outcome),
                        Ok(None) => panic!("SIMD local benchmark produced no positive alignment"),
                        Err(error) => panic!("SIMD local benchmark failed; {error}"),
                    }
                })
            });
        }

        group.finish();
    }
}

criterion_group!(benches, global_benches, local_benches);
criterion_main!(benches);
