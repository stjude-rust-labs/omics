//! Benchmarks for positions.
#![allow(missing_docs)]

use criterion::criterion_group;
use criterion::criterion_main;

////////////////////////////////////////////////////////////////////////////////////////
// Interbase positions
////////////////////////////////////////////////////////////////////////////////////////

pub mod interbase {
    use std::hint::black_box;

    use criterion::Criterion;
    use omics_coordinate::Position;
    use omics_coordinate::position::Number;
    use omics_coordinate::system::Interbase;

    /// A constantly allocated position to remove creation time from the
    /// benchmarks that use it.
    const POSITION: Position<Interbase> = Position::<Interbase>::new(20);

    /// Benchmarks for the [`Position::new()`] method.
    fn new() {
        let _ = Position::<Interbase>::new(black_box(20));
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Basic arithmetic
    ////////////////////////////////////////////////////////////////////////////////////

    /// Benchmarks for the [`Position::checked_add()`] method.
    fn checked_add() -> Number {
        // NOTE: `black_box()` to ensure that the result isn't just inlined.
        POSITION.checked_add(black_box(10)).unwrap().get()
    }

    /// Benchmarks for the [`Position::checked_sub()`] method.
    fn checked_sub() -> Number {
        // NOTE: `black_box()` to ensure that the result isn't just inlined.
        POSITION.checked_sub(black_box(10)).unwrap().get()
    }

    pub fn benches(c: &mut Criterion) {
        c.bench_function("positions::interbase::new", |b| b.iter(new));
        c.bench_function("positions::interbase::checked_add", |b| b.iter(checked_add));
        c.bench_function("positions::interbase::checked_sub", |b| b.iter(checked_sub));
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Base positions
////////////////////////////////////////////////////////////////////////////////////////

pub mod base {
    use std::cell::LazyCell;
    use std::hint::black_box;

    use criterion::Criterion;
    use omics_coordinate::Position;
    use omics_coordinate::position::Number;
    use omics_coordinate::system::Base;

    /// Benchmarks for the [`Position::try_new()`] method.
    fn try_new() {
        let _ = Position::<Base>::try_new(black_box(20)).unwrap();
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Basic arithmetic
    ////////////////////////////////////////////////////////////////////////////////////

    /// Benchmarks for the [`Position::checked_add()`] method.
    fn checked_add() -> Number {
        let position = LazyCell::new(|| Position::<Base>::try_new(20).unwrap());
        position.checked_add(black_box(10)).unwrap().get()
    }

    /// Benchmarks for the [`Position::checked_sub()`] method.
    fn checked_sub() -> Number {
        let position = LazyCell::new(|| Position::<Base>::try_new(20).unwrap());
        position.checked_sub(black_box(10)).unwrap().get()
    }

    pub fn benches(c: &mut Criterion) {
        c.bench_function("positions::base::try_new", |b| b.iter(try_new));
        c.bench_function("positions::base::checked_add", |b| b.iter(checked_add));
        c.bench_function("positions::base::checked_sub", |b| b.iter(checked_sub));
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Registration
////////////////////////////////////////////////////////////////////////////////////////

criterion_group!(benches, interbase::benches, base::benches);
criterion_main!(benches);
