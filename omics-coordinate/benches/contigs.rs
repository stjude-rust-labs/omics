//! Benchmarks for contigs.
#![allow(missing_docs)]

use criterion::criterion_group;
use criterion::criterion_main;

////////////////////////////////////////////////////////////////////////////////////////
// Contigs
////////////////////////////////////////////////////////////////////////////////////////

pub mod interbase {
    use std::hint::black_box;

    use criterion::Criterion;
    use omics_coordinate::Contig;

    /// Benchmarks for the [`Contig::new()`] method.
    fn new() {
        let _ = black_box(Contig::new("seq0"));
    }

    pub fn benches(c: &mut Criterion) {
        c.bench_function("contig::new", |b| b.iter(new));
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Registration
////////////////////////////////////////////////////////////////////////////////////////

criterion_group!(benches, interbase::benches);
criterion_main!(benches);
