#![expect(missing_docs, reason = "integration tests do not need crate docs")]

use omics_variation::copy_number::Change;
use omics_variation::copy_number::Count;
use omics_variation::copy_number::Error;
use omics_variation::copy_number::LogarithmicError;
use omics_variation::copy_number::ParseError;
use omics_variation::copy_number::Ploidy;
use omics_variation::copy_number::PloidyError;
use omics_variation::copy_number::Variant;

#[test]
fn constructs_count_and_ploidy_values() {
    assert_eq!(Count::new(3).get(), 3);
    assert_eq!(Ploidy::HAPLOID.get(), 1);
    assert_eq!(Ploidy::DIPLOID.get(), 2);
    // SAFETY: `4` is a valid nonzero ploidy.
    assert_eq!(Ploidy::try_new(4).unwrap().get(), 4);
    // SAFETY: `4` is a valid nonzero ploidy.
    assert_eq!(Ploidy::try_from(4).unwrap().get(), 4);
    assert_eq!(Ploidy::try_new(0), Err(PloidyError::Zero));
    assert_eq!(Ploidy::try_from(0), Err(PloidyError::Zero));
}

#[test]
fn converts_exact_counts_to_and_from_logarithmic_ratios() {
    let count = Count::new(3);
    let log2 = (1.5_f64).log2();
    let log10 = (1.5_f64).log10();

    // SAFETY: `log2(3 / 2)` maps exactly back to three copies at diploid ploidy.
    let from_log2 = Count::try_from_log2(log2, Ploidy::DIPLOID).unwrap();
    // SAFETY: `log10(3 / 2)` maps exactly back to three copies at diploid ploidy.
    let from_log10 = Count::try_from_log10(log10, Ploidy::DIPLOID).unwrap();

    assert_eq!(from_log2, count);
    assert_eq!(from_log10, count);
    assert!((count.log2(Ploidy::DIPLOID) - log2).abs() < f64::EPSILON);
    assert!((count.log10(Ploidy::DIPLOID) - log10).abs() < f64::EPSILON);
}

#[test]
fn constructs_and_classifies_a_copy_number_variant() {
    // SAFETY: the constructor is expected to accept this valid contig and span.
    let variant = Variant::try_new("seq0", 100, 200, 3).unwrap();
    assert_eq!(variant.contig().as_str(), "seq0");
    assert_eq!(variant.start().get(), 100);
    assert_eq!(variant.end().get(), 200);
    assert_eq!(variant.count(), Count::new(3));
    assert_eq!(variant.change(Count::new(4)), Change::Loss);
    assert_eq!(variant.change(Count::new(3)), Change::Baseline);
    assert_eq!(variant.change(Count::new(2)), Change::Gain);
}

#[test]
fn accepts_complete_copy_loss() {
    // SAFETY: the constructor is expected to accept a complete copy loss.
    let variant = Variant::try_new("seq0", 100, 200, 0).unwrap();
    assert_eq!(variant.count(), Count::new(0));
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

#[test]
fn rejects_invalid_logarithmic_values() {
    assert_eq!(
        Count::try_from_log2(f64::NAN, Ploidy::DIPLOID),
        Err(LogarithmicError::NotANumber)
    );
    assert_eq!(
        Count::try_from_log10(f64::INFINITY, Ploidy::DIPLOID),
        Err(LogarithmicError::PositiveInfinity)
    );
    assert_eq!(
        Count::try_from_log2(-1076.0, Ploidy::DIPLOID),
        Err(LogarithmicError::Underflow)
    );
    assert_eq!(
        Count::try_from_log2(((f64::from(u32::MAX) + 1.0) / 2.0).log2(), Ploidy::DIPLOID),
        Err(LogarithmicError::Overflow)
    );
    assert_eq!(
        Count::try_from_log2((1.25_f64).log2(), Ploidy::DIPLOID),
        Err(LogarithmicError::NonIntegral)
    );
}

#[test]
fn accepts_negative_infinity_and_tolerant_rounding() {
    // SAFETY: negative infinity is defined to round-trip to zero copies.
    let zero_from_log2 = Count::try_from_log2(f64::NEG_INFINITY, Ploidy::DIPLOID).unwrap();
    // SAFETY: negative infinity is defined to round-trip to zero copies.
    let zero_from_log10 = Count::try_from_log10(f64::NEG_INFINITY, Ploidy::DIPLOID).unwrap();
    let within_tolerance = ((3.0 + 12.0 * f64::EPSILON) / 2.0).log2();

    // SAFETY: this logarithmic ratio stays within the approved integral tolerance.
    let tolerated = Count::try_from_log2(within_tolerance, Ploidy::DIPLOID).unwrap();

    assert_eq!(zero_from_log2, Count::new(0));
    assert_eq!(zero_from_log10, Count::new(0));
    assert_eq!(tolerated, Count::new(3));
    assert!(Count::new(0).log2(Ploidy::DIPLOID).is_infinite());
    assert!(Count::new(0).log2(Ploidy::DIPLOID).is_sign_negative());
    assert!(Count::new(0).log10(Ploidy::DIPLOID).is_infinite());
    assert!(Count::new(0).log10(Ploidy::DIPLOID).is_sign_negative());
}
