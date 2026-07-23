#![expect(missing_docs, reason = "integration tests do not need crate docs")]

use omics_variation::copy_number::{Change, Error, ParseError, Variant};

#[test]
fn constructs_and_classifies_a_copy_number_variant() {
    // SAFETY: the constructor is expected to accept this valid contig and span.
    let variant = Variant::try_new("seq0", 100, 200, 3).unwrap();
    assert_eq!(variant.contig().as_str(), "seq0");
    assert_eq!(variant.start().get(), 100);
    assert_eq!(variant.end().get(), 200);
    assert_eq!(variant.copies(), 3);
    assert_eq!(variant.change(4), Change::Loss);
    assert_eq!(variant.change(3), Change::Baseline);
    assert_eq!(variant.change(2), Change::Gain);
}

#[test]
fn accepts_complete_copy_loss() {
    // SAFETY: the constructor is expected to accept a complete copy loss.
    let variant = Variant::try_new("seq0", 100, 200, 0).unwrap();
    assert_eq!(variant.copies(), 0);
}

#[test]
fn rejects_empty_and_reversed_regions() {
    assert!(matches!(
        Variant::try_new("seq0", 100, 100, 2),
        Err(Error::EmptyRegion)
    ));
    assert!(matches!(
        Variant::try_new("seq0", 200, 100, 2),
        Err(Error::ReversedRegion)
    ));
}

#[test]
fn parses_and_displays_the_canonical_form() {
    // SAFETY: the parser is expected to accept the canonical copy-number form.
    let variant = "seq0:100-200(i):3".parse::<Variant>().unwrap();
    assert_eq!(variant.to_string(), "seq0:100-200(i):3");
}

#[test]
fn rejects_a_missing_interbase_qualifier() {
    let err = "seq0:100-200:3".parse::<Variant>().unwrap_err();
    assert!(matches!(err, ParseError::Qualifier { .. }));
}

#[test]
fn rejects_a_malformed_position_range() {
    let err = "seq0:100200(i):3".parse::<Variant>().unwrap_err();
    assert!(matches!(err, ParseError::Range { .. }));
}

#[test]
fn rejects_a_non_integral_count() {
    let err = "seq0:100-200(i):3.5".parse::<Variant>().unwrap_err();
    assert!(matches!(err, ParseError::Copies { .. }));
}
