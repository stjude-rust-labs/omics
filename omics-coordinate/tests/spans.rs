#![expect(missing_docs, reason = "integration tests do not need crate docs")]

use omics_coordinate::Position;
use omics_coordinate::Span;
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
