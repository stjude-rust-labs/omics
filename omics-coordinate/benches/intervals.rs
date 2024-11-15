//! Benchmarks for intervals.
#![allow(missing_docs)]

use criterion::criterion_group;
use criterion::criterion_main;

////////////////////////////////////////////////////////////////////////////////////////
// Interbase intervals
////////////////////////////////////////////////////////////////////////////////////////

pub mod interbase {

    use criterion::Criterion;
    use omics_coordinate::Coordinate;
    use omics_coordinate::Interval;
    use omics_coordinate::system::Interbase;

    /// Benchmarks for the [`Interval::try_new()`] method.
    ///
    /// This benchmark intentionally includes the creation of the constituent
    /// coordinates. This is because creation a new interval _requires_ taking
    /// ownership of the underlying coordinates, so each new interval requires
    /// reallocation of all parts.
    ///
    /// In addition, though we might could preallocate the coordinates required,
    /// we're not really interested in benchmarking how long it takes to
    /// validate the interval—we want to know how quickly intervals can be
    /// constructed, say, for the creation of an interval tree.
    fn try_new() {
        let start = Coordinate::<Interbase>::try_new("seq0", "+", 10).unwrap();
        let end = Coordinate::<Interbase>::try_new("seq0", "+", 20).unwrap();
        Interval::try_new(start, end).unwrap();
    }

    pub fn benches(c: &mut Criterion) {
        c.bench_function("intervals::interbase::try_new", |b| b.iter(try_new));
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Registration
////////////////////////////////////////////////////////////////////////////////////////

criterion_group!(benches, interbase::benches);
criterion_main!(benches);
