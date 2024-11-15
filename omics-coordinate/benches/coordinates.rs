//! Benchmarks for coordinates.
#![allow(missing_docs)]

use criterion::criterion_group;
use criterion::criterion_main;

////////////////////////////////////////////////////////////////////////////////////////
// Interbase coordinates
////////////////////////////////////////////////////////////////////////////////////////

pub mod interbase {
    use criterion::Criterion;
    use omics_coordinate::Coordinate;
    use omics_coordinate::system::Interbase;

    /// Benchmarks for the [`Coordinate::try_new()`] method.
    fn try_new() {
        Coordinate::<Interbase>::try_new("seq0", "+", 10).unwrap();
    }

    pub fn benches(c: &mut Criterion) {
        c.bench_function("coordinates::interbase::try_new", |b| b.iter(try_new));
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Registration
////////////////////////////////////////////////////////////////////////////////////////

criterion_group!(benches, interbase::benches);
criterion_main!(benches);
