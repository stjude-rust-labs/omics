#![expect(missing_docs, reason = "integration tests do not need crate docs")]

use omics_coordinate::Coordinate;
use omics_coordinate::Interval;
use omics_coordinate::Position;
use omics_coordinate::Span;
use omics_coordinate::Strand;
use omics_coordinate::interval;
use omics_coordinate::span::ClampError;
use omics_coordinate::span::Direction;
use omics_coordinate::system::Base;
use omics_coordinate::system::Interbase;

#[test]
fn preserves_all_endpoint_directions() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let stationary = Span::<Interbase>::try_new(10, 10)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert_eq!(ascending.direction(), Direction::Ascending);
    assert_eq!(stationary.direction(), Direction::Stationary);
    assert_eq!(descending.direction(), Direction::Descending);
    assert_eq!(descending.reversed(), ascending);

    Ok(())
}

#[test]
fn keeps_coordinate_system_semantics() -> Result<(), Box<dyn std::error::Error>> {
    let base = Span::<Base>::try_new(10, 10)?;
    let interbase = Span::<Interbase>::try_new(10, 10)?;

    assert!(!base.is_empty());
    assert_eq!(base.count_entities(), 1);
    assert!(interbase.is_empty());
    assert_eq!(interbase.count_entities(), 0);
    assert!(Span::<Base>::try_new(0, 1).is_err());

    Ok(())
}

#[test]
fn parses_and_displays_directed_spans() -> Result<(), Box<dyn std::error::Error>> {
    let span = "20-10".parse::<Span<Interbase>>()?;

    assert_eq!(span.direction(), Direction::Descending);
    assert_eq!(span.to_string(), "20-10");

    Ok(())
}

#[test]
fn converts_from_typed_positions() {
    let start = Position::<Interbase>::new(20);
    let end = Position::<Interbase>::new(10);
    let span = Span::from((start, end));

    assert_eq!(span.into_positions(), (start, end));
}

#[test]
fn contains_positions_in_both_directions() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert!(ascending.contains_position(&Position::<Interbase>::new(10)));
    assert!(ascending.contains_position(&Position::<Interbase>::new(15)));
    assert!(ascending.contains_position(&Position::<Interbase>::new(20)));
    assert!(!ascending.contains_position(&Position::<Interbase>::new(9)));
    assert!(!ascending.contains_position(&Position::<Interbase>::new(21)));

    assert!(descending.contains_position(&Position::<Interbase>::new(20)));
    assert!(descending.contains_position(&Position::<Interbase>::new(15)));
    assert!(descending.contains_position(&Position::<Interbase>::new(10)));
    assert!(!descending.contains_position(&Position::<Interbase>::new(9)));
    assert!(!descending.contains_position(&Position::<Interbase>::new(21)));

    Ok(())
}

#[test]
fn contains_entities_with_base_semantics() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Base>::try_new(10, 20)?;
    let descending = Span::<Base>::try_new(20, 10)?;

    assert!(!ascending.contains_entity(&Position::<Base>::try_new(9)?));
    assert!(ascending.contains_entity(&Position::<Base>::try_new(10)?));
    assert!(ascending.contains_entity(&Position::<Base>::try_new(15)?));
    assert!(ascending.contains_entity(&Position::<Base>::try_new(20)?));
    assert!(!ascending.contains_entity(&Position::<Base>::try_new(21)?));

    assert!(!descending.contains_entity(&Position::<Base>::try_new(9)?));
    assert!(descending.contains_entity(&Position::<Base>::try_new(10)?));
    assert!(descending.contains_entity(&Position::<Base>::try_new(15)?));
    assert!(descending.contains_entity(&Position::<Base>::try_new(20)?));
    assert!(!descending.contains_entity(&Position::<Base>::try_new(21)?));

    Ok(())
}

#[test]
fn contains_entities_with_interbase_traversal_semantics() -> Result<(), Box<dyn std::error::Error>>
{
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let stationary = Span::<Interbase>::try_new(10, 10)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert!(!ascending.contains_entity(&Position::<Base>::try_new(10)?));
    assert!(ascending.contains_entity(&Position::<Base>::try_new(11)?));
    assert!(ascending.contains_entity(&Position::<Base>::try_new(15)?));
    assert!(ascending.contains_entity(&Position::<Base>::try_new(20)?));
    assert!(!ascending.contains_entity(&Position::<Base>::try_new(21)?));

    assert!(!stationary.contains_entity(&Position::<Base>::try_new(10)?));

    assert!(!descending.contains_entity(&Position::<Base>::try_new(10)?));
    assert!(descending.contains_entity(&Position::<Base>::try_new(11)?));
    assert!(descending.contains_entity(&Position::<Base>::try_new(15)?));
    assert!(descending.contains_entity(&Position::<Base>::try_new(20)?));
    assert!(!descending.contains_entity(&Position::<Base>::try_new(21)?));

    Ok(())
}

#[test]
fn exposes_system_aliases() -> Result<(), Box<dyn std::error::Error>> {
    let base = omics_coordinate::span::base::Span::try_new(10, 12)?;
    let interbase = omics_coordinate::span::interbase::Span::try_new(10, 12)?;

    assert_eq!(base.start().get(), 10);
    assert_eq!(interbase.end().get(), 12);

    Ok(())
}

#[test]
fn resolves_offsets_in_traversal_direction() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert_eq!(ascending.position_offset(&Position::new(15)), Some(5));
    assert_eq!(descending.position_offset(&Position::new(15)), Some(5));
    assert_eq!(ascending.position_offset(&Position::new(9)), None);
    assert_eq!(descending.position_offset(&Position::new(21)), None);
    assert_eq!(ascending.position_at_offset(5), Some(Position::new(15)));
    assert_eq!(descending.position_at_offset(5), Some(Position::new(15)));
    assert_eq!(descending.position_at_offset(11), None);

    Ok(())
}

#[test]
fn clamps_same_direction_overlaps() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let ascending_operand = Span::<Interbase>::try_new(15, 25)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;
    let descending_operand = Span::<Interbase>::try_new(18, 8)?;

    assert_eq!(
        ascending.clamp(ascending_operand)?,
        Span::<Interbase>::try_new(15, 20)?
    );
    assert_eq!(
        descending.clamp(descending_operand)?,
        Span::<Interbase>::try_new(18, 10)?
    );

    Ok(())
}

#[test]
fn clamps_stationary_spans_against_either_direction() -> Result<(), Box<dyn std::error::Error>> {
    let point = Span::<Interbase>::try_new(15, 15)?;
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert_eq!(point.clone().clamp(ascending)?, point);
    assert_eq!(point.clone().clamp(descending)?, point);

    Ok(())
}

#[test]
fn clamps_when_operand_is_stationary() -> Result<(), Box<dyn std::error::Error>> {
    let point = Span::<Interbase>::try_new(15, 15)?;
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert_eq!(ascending.clamp(point.clone())?, point);
    assert_eq!(descending.clamp(point.clone())?, point);

    Ok(())
}

#[test]
fn rejects_disjoint_spans() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let operand = Span::<Interbase>::try_new(21, 30)?;

    assert_eq!(
        ascending.clamp(operand).unwrap_err(),
        ClampError::Disjoint {
            original_start: 10,
            original_end: 20,
            operand_start: 21,
            operand_end: 30,
        }
    );

    Ok(())
}

#[test]
fn rejects_incompatible_directions() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let descending = Span::<Interbase>::try_new(20, 10)?;

    assert_eq!(
        ascending.clone().clamp(descending.clone()).unwrap_err(),
        ClampError::DirectionMismatch {
            original_start: 10,
            original_end: 20,
            original_direction: Direction::Ascending,
            operand_start: 20,
            operand_end: 10,
            operand_direction: Direction::Descending,
        }
    );
    assert_eq!(
        descending.clamp(ascending).unwrap_err(),
        ClampError::DirectionMismatch {
            original_start: 20,
            original_end: 10,
            original_direction: Direction::Descending,
            operand_start: 10,
            operand_end: 20,
            operand_direction: Direction::Ascending,
        }
    );

    Ok(())
}

#[test]
fn localizes_spans_with_compatible_strands() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let positive = Interval::try_new("seq0", "+", ascending.clone())?;

    assert_eq!(positive.contig().as_str(), "seq0");
    assert_eq!(positive.strand(), Strand::Positive);
    assert_eq!(positive.span(), &ascending);

    let descending = ascending.reversed();
    let negative = Interval::try_new("seq0", "-", descending.clone())?;

    assert_eq!(negative.strand(), Strand::Negative);
    assert_eq!(negative.span(), &descending);

    Ok(())
}

#[test]
fn rejects_direction_and_strand_mismatches() -> Result<(), Box<dyn std::error::Error>> {
    let ascending = Span::<Interbase>::try_new(10, 20)?;
    let descending = ascending.clone().reversed();

    assert_eq!(
        Interval::try_new("seq0", "-", ascending).unwrap_err(),
        interval::Error::DirectionMismatch {
            strand: Strand::Negative,
            direction: Direction::Ascending,
        }
    );
    assert_eq!(
        Interval::try_new("seq0", "+", descending).unwrap_err(),
        interval::Error::DirectionMismatch {
            strand: Strand::Positive,
            direction: Direction::Descending,
        }
    );

    Ok(())
}

#[test]
fn accepts_stationary_spans_on_both_strands() -> Result<(), Box<dyn std::error::Error>> {
    let span = Span::<Interbase>::try_new(10, 10)?;
    let positive = Interval::try_new("seq0", "+", span.clone())?;
    let negative = Interval::try_new("seq0", "-", span.clone())?;

    assert_eq!(positive.span(), &span);
    assert_eq!(negative.span(), &span);

    Ok(())
}

#[test]
fn converts_coordinate_pairs() -> Result<(), Box<dyn std::error::Error>> {
    let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    let end = Coordinate::<Interbase>::try_new("seq0", "+", 20)?;
    let interval = Interval::try_from((start, end))?;

    assert_eq!(interval.span().to_string(), "10-20");

    Ok(())
}

#[test]
fn rejects_coordinate_pairs_with_mismatched_context() -> Result<(), Box<dyn std::error::Error>> {
    let start = Coordinate::<Interbase>::try_new("seq0", "+", 10)?;
    let other_contig = Coordinate::<Interbase>::try_new("seq1", "+", 20)?;
    let other_strand = Coordinate::<Interbase>::try_new("seq0", "-", 10)?;

    assert!(matches!(
        Interval::try_from((start.clone(), other_contig)),
        Err(interval::Error::Nonsensical(
            interval::NonsensicalError::MismatchedContigs { .. }
        ))
    ));
    assert!(matches!(
        Interval::try_from((start, other_strand)),
        Err(interval::Error::Nonsensical(
            interval::NonsensicalError::MismatchedStrands { .. }
        ))
    ));

    Ok(())
}

#[test]
fn consumes_interval_into_span() -> Result<(), Box<dyn std::error::Error>> {
    let span = Span::<Interbase>::try_new(10, 20)?;
    let interval = Interval::try_new("seq0", "+", span.clone())?;

    assert_eq!(interval.into_span(), span);

    Ok(())
}

#[test]
fn interval_strings_round_trip_with_colon_contigs() -> Result<(), Box<dyn std::error::Error>> {
    for value in ["seq0:+:10-20", "assembly:chr:1:-:20-10"] {
        let interval = value.parse::<Interval<Interbase>>()?;

        assert_eq!(interval.to_string(), value);
        assert_eq!(
            interval.to_string().parse::<Interval<Interbase>>()?,
            interval
        );
    }

    Ok(())
}

#[test]
fn delegates_geometry_to_spans() -> Result<(), Box<dyn std::error::Error>> {
    let span = Span::<Interbase>::try_new(10, 20)?;
    let interval = Interval::try_new("seq0", "+", span.clone())?;
    let coordinate = Coordinate::<Interbase>::try_new("seq0", "+", 15)?;
    let entity = Position::<Base>::try_new(15)?;

    assert_eq!(
        interval.contains_coordinate(&coordinate),
        span.contains_position(coordinate.position())
    );
    assert_eq!(
        interval.contains_entity(&Coordinate::try_new("seq0", "+", 15)?),
        span.contains_entity(&entity)
    );
    assert_eq!(interval.count_entities(), span.count_entities());
    assert_eq!(
        interval.coordinate_offset(&coordinate),
        span.position_offset(coordinate.position())
    );
    assert_eq!(
        interval
            .coordinate_at_offset(5)
            .map(|coordinate| coordinate.into_position()),
        span.position_at_offset(5)
    );

    let operand_span = Span::<Interbase>::try_new(15, 25)?;
    let operand = Interval::try_new("seq0", "+", operand_span.clone())?;
    let clamped = interval.clone().clamp(operand)?;

    assert_eq!(clamped.span(), &span.clone().clamp(operand_span)?);

    let reversed = interval.reverse_complement();
    assert_eq!(reversed.strand(), Strand::Negative);
    assert_eq!(reversed.span(), &span.reversed());

    Ok(())
}
