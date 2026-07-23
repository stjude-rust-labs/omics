//! Benchmarks for spans and intervals.
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
    use omics_coordinate::Position;
    use omics_coordinate::Span;
    use omics_coordinate::system::Base;
    use omics_coordinate::system::Interbase;

    fn ascending_span() -> Span<Interbase> {
        // SAFETY: the benchmark endpoints form a valid interbase span.
        Span::<Interbase>::try_new(10, 20).unwrap()
    }

    fn descending_span() -> Span<Interbase> {
        // SAFETY: the benchmark endpoints form a valid interbase span.
        Span::<Interbase>::try_new(20, 10).unwrap()
    }

    fn descending_operand_span() -> Span<Interbase> {
        // SAFETY: the benchmark endpoints form a valid interbase span.
        Span::<Interbase>::try_new(18, 8).unwrap()
    }

    fn ascending_operand_span() -> Span<Interbase> {
        // SAFETY: the benchmark endpoints form a valid interbase span.
        Span::<Interbase>::try_new(15, 25).unwrap()
    }

    fn empty_span() -> Span<Interbase> {
        // SAFETY: equal interbase endpoints are valid and represent an empty span.
        Span::<Interbase>::try_new(10, 10).unwrap()
    }

    fn interval() -> Interval<Interbase> {
        // SAFETY: the positive strand accepts ascending spans.
        Interval::try_new("seq0", "+", ascending_span()).unwrap()
    }

    fn operand_interval() -> Interval<Interbase> {
        // SAFETY: the positive strand accepts ascending spans.
        Interval::try_new("seq0", "+", ascending_operand_span()).unwrap()
    }

    /// Benchmarks span construction from raw position numbers.
    fn span_try_new() {
        // SAFETY: the benchmark endpoints form a valid interbase span.
        black_box(Span::<Interbase>::try_new(black_box(10), black_box(20)).unwrap());
    }

    /// Benchmarks span parsing.
    fn span_parse() {
        // SAFETY: the benchmark input is a valid interbase span string.
        black_box(black_box("20-10").parse::<Span<Interbase>>().unwrap());
    }

    /// Benchmarks span display formatting.
    fn span_display(span: &Span<Interbase>) {
        black_box(black_box(span).to_string());
    }

    /// Benchmarks span direction queries.
    fn span_direction(span: &Span<Interbase>) {
        black_box(black_box(span).direction());
    }

    /// Benchmarks span emptiness queries.
    fn span_is_empty(span: &Span<Interbase>) {
        black_box(black_box(span).is_empty());
    }

    /// Benchmarks span position containment.
    fn span_contains_position(span: &Span<Interbase>, position: &Position<Interbase>) {
        black_box(black_box(span).contains_position(black_box(position)));
    }

    /// Benchmarks span offset lookup from the start.
    fn span_position_offset(span: &Span<Interbase>, position: &Position<Interbase>) {
        black_box(black_box(span).position_offset(black_box(position)));
    }

    /// Benchmarks span offset lookup into the span.
    fn span_position_at_offset(span: &Span<Interbase>) {
        black_box(black_box(span).position_at_offset(black_box(5)));
    }

    /// Benchmarks interval construction from raw coordinate parts.
    fn coordinate_pair_try_from_raw() {
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

    /// Benchmarks interval construction from contextual data and a span.
    fn interval_try_new() {
        // SAFETY: the benchmark endpoints form a valid interbase span.
        let span = Span::<Interbase>::try_new(black_box(10), black_box(20)).unwrap();
        // SAFETY: the positive strand accepts ascending spans.
        black_box(Interval::try_new(black_box("seq0"), black_box("+"), span).unwrap());
    }

    /// Benchmarks interval coordinate containment.
    fn interval_contains_coordinate(
        interval: &Interval<Interbase>,
        coordinate: &Coordinate<Interbase>,
    ) {
        black_box(black_box(interval).contains_coordinate(black_box(coordinate)));
    }

    /// Benchmarks interval entity containment.
    fn interval_contains_entity(interval: &Interval<Interbase>, coordinate: &Coordinate<Base>) {
        black_box(black_box(interval).contains_entity(black_box(coordinate)));
    }

    /// Benchmarks interval entity counting.
    fn interval_count_entities(interval: &Interval<Interbase>) {
        black_box(black_box(interval).count_entities());
    }

    /// Benchmarks interval coordinate offsets.
    fn interval_coordinate_offset(
        interval: &Interval<Interbase>,
        coordinate: &Coordinate<Interbase>,
    ) {
        black_box(black_box(interval).coordinate_offset(black_box(coordinate)));
    }

    /// Benchmarks interval coordinate lookup by offset.
    fn interval_coordinate_at_offset(interval: &Interval<Interbase>) {
        black_box(black_box(interval).coordinate_at_offset(black_box(5)));
    }

    pub fn benches(c: &mut Criterion) {
        let ascending = ascending_span();
        let descending = descending_span();
        let descending_operand = descending_operand_span();
        let empty = empty_span();
        let position = Position::<Interbase>::new(15);
        let interval = interval();
        let operand = operand_interval();
        // SAFETY: the benchmark coordinate points to a valid entity on the interval
        // contig.
        let entity = Coordinate::<Base>::try_new("seq0", "+", 15).unwrap();
        // SAFETY: the benchmark coordinate shares the interval context and lies within
        // the span.
        let coordinate = Coordinate::<Interbase>::try_new("seq0", "+", 15).unwrap();

        c.bench_function("intervals::interbase::span_try_new", |b| {
            b.iter(span_try_new)
        });
        c.bench_function("intervals::interbase::span_parse", |b| b.iter(span_parse));
        c.bench_function("intervals::interbase::span_display", |b| {
            b.iter(|| span_display(&descending))
        });
        c.bench_function("intervals::interbase::span_direction", |b| {
            b.iter(|| span_direction(&descending))
        });
        c.bench_function("intervals::interbase::span_is_empty", |b| {
            b.iter(|| span_is_empty(&empty))
        });
        c.bench_function("intervals::interbase::span_contains_position", |b| {
            b.iter(|| span_contains_position(&ascending, &position))
        });
        c.bench_function("intervals::interbase::span_position_offset", |b| {
            b.iter(|| span_position_offset(&descending, &position))
        });
        c.bench_function("intervals::interbase::span_position_at_offset", |b| {
            b.iter(|| span_position_at_offset(&descending))
        });
        c.bench_function("intervals::interbase::span_reversed", |b| {
            b.iter_batched(
                descending_span,
                |span| black_box(span.reversed()),
                BatchSize::SmallInput,
            )
        });
        c.bench_function("intervals::interbase::span_clamp", |b| {
            b.iter_batched(
                || (descending.clone(), descending_operand.clone()),
                |(span, operand)| {
                    // SAFETY: the prepared spans overlap and move in compatible directions.
                    black_box(span.clamp(operand).unwrap())
                },
                BatchSize::SmallInput,
            )
        });
        c.bench_function("intervals::interbase::coordinate_pair_try_from_raw", |b| {
            b.iter(coordinate_pair_try_from_raw)
        });
        c.bench_function(
            "intervals::interbase::coordinate_pair_try_from_prepared",
            |b| {
                b.iter_batched(
                    || {
                        // SAFETY: the benchmark inputs form a valid interbase coordinate.
                        let start = Coordinate::<Interbase>::try_new("seq0", "+", 10).unwrap();
                        // SAFETY: the benchmark inputs form a valid interbase coordinate.
                        let end = Coordinate::<Interbase>::try_new("seq0", "+", 20).unwrap();
                        (start, end)
                    },
                    |(start, end)| {
                        // SAFETY: the prepared endpoints form a valid interval.
                        black_box(Interval::try_from((start, end)).unwrap())
                    },
                    BatchSize::SmallInput,
                )
            },
        );
        c.bench_function("intervals::interbase::interval_try_new", |b| {
            b.iter(interval_try_new)
        });
        c.bench_function("intervals::interbase::try_new_prepared", |b| {
            b.iter_batched(
                ascending_span,
                |span| {
                    // SAFETY: the positive strand accepts the prepared ascending span.
                    black_box(Interval::try_new("seq0", "+", span).unwrap())
                },
                BatchSize::SmallInput,
            )
        });
        c.bench_function("intervals::interbase::contains_coordinate", |b| {
            b.iter(|| interval_contains_coordinate(&interval, &coordinate))
        });
        c.bench_function("intervals::interbase::contains_entity", |b| {
            b.iter(|| interval_contains_entity(&interval, &entity))
        });
        c.bench_function("intervals::interbase::count_entities", |b| {
            b.iter(|| interval_count_entities(&interval))
        });
        c.bench_function("intervals::interbase::coordinate_offset", |b| {
            b.iter(|| interval_coordinate_offset(&interval, &coordinate))
        });
        c.bench_function("intervals::interbase::coordinate_at_offset", |b| {
            b.iter(|| interval_coordinate_at_offset(&interval))
        });
        c.bench_function("intervals::interbase::interval_clamp", |b| {
            b.iter_batched(
                || (interval.clone(), operand.clone()),
                |(interval, operand)| {
                    // SAFETY: the prepared intervals share context and overlap.
                    black_box(interval.clamp(operand).unwrap())
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
