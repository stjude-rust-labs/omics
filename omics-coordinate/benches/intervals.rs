//! Benchmarks for intervals.
#![allow(missing_docs)]

use criterion::criterion_group;
use criterion::criterion_main;

////////////////////////////////////////////////////////////////////////////////////////
// Interbase intervals
////////////////////////////////////////////////////////////////////////////////////////

pub mod interbase {

    use std::hint::black_box;

    use criterion::BatchSize;
    use criterion::Criterion;
    use omics_coordinate::Coordinate;
    use omics_coordinate::Interval;
    use omics_coordinate::system::Interbase;

    /// Benchmarks interval construction from raw coordinate parts.
    fn from_raw_coordinates() {
        // SAFETY: the benchmark inputs form a valid interbase coordinate.
        let start =
            Coordinate::<Interbase>::try_new(black_box("seq0"), black_box("+"), black_box(10))
                .unwrap();
        // SAFETY: the benchmark inputs form a valid interbase coordinate.
        let end =
            Coordinate::<Interbase>::try_new(black_box("seq0"), black_box("+"), black_box(20))
                .unwrap();
        // SAFETY: the endpoints share metadata and form a positively sized interval.
        black_box(Interval::try_from((start, end)).unwrap());
    }

    pub fn benches(c: &mut Criterion) {
        c.bench_function("intervals::interbase::from_raw_coordinates", |b| {
            b.iter(from_raw_coordinates)
        });

        // SAFETY: the benchmark inputs form a valid interbase coordinate.
        let start = Coordinate::<Interbase>::try_new("seq0", "+", 10).unwrap();
        // SAFETY: the benchmark inputs form a valid interbase coordinate.
        let end = Coordinate::<Interbase>::try_new("seq0", "+", 20).unwrap();

        c.bench_function("intervals::interbase::try_new_prepared", |b| {
            b.iter_batched(
                || (start.clone(), end.clone()),
                |(start, end)| {
                    // SAFETY: the prepared endpoints form a valid interval.
                    black_box(Interval::try_from((start, end)).unwrap())
                },
                BatchSize::SmallInput,
            )
        });
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Registration
////////////////////////////////////////////////////////////////////////////////////////

criterion_group!(benches, interbase::benches);
criterion_main!(benches);
