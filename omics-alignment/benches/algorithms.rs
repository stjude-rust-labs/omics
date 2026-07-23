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

/// Equal input lengths used to measure quadratic scaling.
const SIZES: [usize; 3] = [16, 64, 256];

/// Builds equal-length, mostly matching DNA-like benchmark inputs.
fn fixture(size: usize) -> (Vec<u8>, Vec<u8>) {
    let reference = (0..size)
        .map(|index| b"ACGT"[index % 4])
        .collect::<Vec<_>>();
    let mut query = reference.clone();

    for start in (0..size).step_by(16) {
        query[start + 3] = b'A';
        query[start + 8..start + 12].rotate_left(1);
    }

    (reference, query)
}

/// Constructs the scoring configuration shared by every benchmark.
fn scoring() -> Scoring {
    // SAFETY: the benchmark scores satisfy all Scoring sign constraints.
    Scoring::try_new(2, -3, -2, -1).unwrap()
}

/// Returns the number of dynamic-programming cells for equal-length inputs.
fn cells(size: usize) -> u64 {
    ((size + 1) * (size + 1)) as u64
}

/// Registers end-to-end global-alignment scaling benchmarks.
fn global_benches(c: &mut Criterion) {
    let scoring = scoring();
    let mut group = c.benchmark_group("algorithms/global");

    for size in SIZES {
        let (reference, query) = fixture(size);
        group.throughput(Throughput::Elements(cells(size)));
        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
            b.iter(|| {
                match global(
                    black_box(reference.as_slice()),
                    black_box(query.as_slice()),
                    scoring,
                ) {
                    Ok(outcome) => black_box(outcome),
                    Err(error) => panic!("global benchmark failed; {error}"),
                }
            })
        });
    }

    group.finish();
}

/// Registers end-to-end local-alignment scaling benchmarks.
fn local_benches(c: &mut Criterion) {
    let scoring = scoring();
    let mut group = c.benchmark_group("algorithms/local");

    for size in SIZES {
        let (reference, query) = fixture(size);
        group.throughput(Throughput::Elements(cells(size)));
        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
            b.iter(|| {
                match local(
                    black_box(reference.as_slice()),
                    black_box(query.as_slice()),
                    scoring,
                ) {
                    Ok(Some(outcome)) => black_box(outcome),
                    Ok(None) => panic!("local benchmark produced no positive alignment"),
                    Err(error) => panic!("local benchmark failed; {error}"),
                }
            })
        });
    }

    group.finish();
}

criterion_group!(benches, global_benches, local_benches);
criterion_main!(benches);
