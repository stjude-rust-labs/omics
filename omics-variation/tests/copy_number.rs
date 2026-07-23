#![expect(missing_docs, reason = "integration tests do not need crate docs")]

use omics_variation::CopyNumberChange;
use omics_variation::CopyNumberCount;
use omics_variation::CopyNumberPloidy;
use omics_variation::CopyNumberVariant;
use omics_variation::copy_number::Change;
use omics_variation::copy_number::Count;
use omics_variation::copy_number::Error;
use omics_variation::copy_number::LogarithmicError;
use omics_variation::copy_number::LogarithmicVariantError;
use omics_variation::copy_number::ParseError;
use omics_variation::copy_number::Ploidy;
use omics_variation::copy_number::PloidyError;
use omics_variation::copy_number::Variant;

#[test]
fn imports_top_level_copy_number_aliases() {
    let variant = CopyNumberVariant::try_new("seq0", 100, 200, 3, CopyNumberPloidy::DIPLOID);
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // ploidy.
    let variant = variant.unwrap();
    let count = variant.count();
    let log2 = variant.log2();
    let log10 = variant.log10();
    let parsed = "seq0:100-200(i):3/2".parse::<CopyNumberVariant>();
    // SAFETY: the canonical top-level alias form contains a valid region, count,
    // and reference ploidy.
    let parsed = parsed.unwrap();

    assert_eq!(count, CopyNumberCount::new(3));
    assert_eq!(variant.ploidy(), CopyNumberPloidy::DIPLOID);
    let count_from_log2 = CopyNumberCount::try_from_log2(log2, CopyNumberPloidy::DIPLOID);
    // SAFETY: `log2(3 / 2)` round-trips to the original typed copy-number count.
    let count_from_log2 = count_from_log2.unwrap();
    let count_from_log10 = CopyNumberCount::try_from_log10(log10, CopyNumberPloidy::DIPLOID);
    // SAFETY: `log10(3 / 2)` round-trips to the original typed copy-number count.
    let count_from_log10 = count_from_log10.unwrap();
    let variant_from_log2 =
        CopyNumberVariant::try_from_log2("seq0", 100, 200, log2, CopyNumberPloidy::DIPLOID);
    // SAFETY: `log2(3 / 2)` round-trips to the original top-level copy-number
    // variant.
    let variant_from_log2 = variant_from_log2.unwrap();
    let variant_from_log10 =
        CopyNumberVariant::try_from_log10("seq0", 100, 200, log10, CopyNumberPloidy::DIPLOID);
    // SAFETY: `log10(3 / 2)` round-trips to the original top-level copy-number
    // variant.
    let variant_from_log10 = variant_from_log10.unwrap();

    assert_eq!(count_from_log2, count);
    assert_eq!(count_from_log10, count);
    assert_eq!(variant.change(), CopyNumberChange::Gain);
    assert_eq!(parsed, variant);
    assert_eq!(variant_from_log2, variant);
    assert_eq!(variant_from_log10, variant);
    assert_eq!(variant.to_string(), "seq0:100-200(i):3/2");
}

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
fn stores_reference_ploidy_and_derives_change() {
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // diploid ploidy.
    let variant = Variant::try_new("seq0", 100, 200, 3, Ploidy::DIPLOID).unwrap();

    assert_eq!(variant.contig().as_str(), "seq0");
    assert_eq!(variant.start().get(), 100);
    assert_eq!(variant.end().get(), 200);
    assert_eq!(variant.count(), Count::new(3));
    assert_eq!(variant.ploidy(), Ploidy::DIPLOID);
    assert_eq!(variant.change(), Change::Gain);
    assert_eq!(variant.log2(), Count::new(3).log2(Ploidy::DIPLOID));
    assert_eq!(variant.log10(), Count::new(3).log10(Ploidy::DIPLOID));
}

#[test]
fn constructs_variants_from_exact_logarithmic_ratios() {
    let diploid_log2 = Variant::try_from_log2("seq0", 100, 200, (1.5_f64).log2(), Ploidy::DIPLOID);
    // SAFETY: `log2(3 / 2)` must reconstruct three copies at diploid ploidy.
    let diploid_log2 = diploid_log2.unwrap();
    let haploid_log10 = Variant::try_from_log10("seq0", 300, 400, 0.0, Ploidy::HAPLOID);
    // SAFETY: `log10(1 / 1)` must reconstruct one copy at haploid ploidy.
    let haploid_log10 = haploid_log10.unwrap();
    // SAFETY: `4` is a valid nonzero ploidy.
    let tetraploid = Ploidy::try_new(4).unwrap();
    let tetraploid_log2 = Variant::try_from_log2("seq0", 500, 600, (1.5_f64).log2(), tetraploid);
    // SAFETY: `log2(6 / 4)` must reconstruct six copies while preserving tetraploid
    // identity.
    let tetraploid_log2 = tetraploid_log2.unwrap();

    assert_eq!(diploid_log2.count(), Count::new(3));
    assert_eq!(diploid_log2.ploidy(), Ploidy::DIPLOID);
    assert_eq!(diploid_log2.change(), Change::Gain);
    assert_eq!(haploid_log10.count(), Count::new(1));
    assert_eq!(haploid_log10.ploidy(), Ploidy::HAPLOID);
    assert_eq!(haploid_log10.change(), Change::Reference);
    assert_eq!(tetraploid_log2.count(), Count::new(6));
    assert_eq!(tetraploid_log2.ploidy(), tetraploid);
    assert_eq!(tetraploid_log2.change(), Change::Gain);
}

#[test]
fn constructs_logarithmic_variants_for_complete_copy_loss() {
    let zero_from_log2 =
        Variant::try_from_log2("seq0", 100, 200, f64::NEG_INFINITY, Ploidy::DIPLOID);
    // SAFETY: negative infinity is defined to map to zero copies.
    let zero_from_log2 = zero_from_log2.unwrap();
    let zero_from_log10 =
        Variant::try_from_log10("seq0", 300, 400, f64::NEG_INFINITY, Ploidy::DIPLOID);
    // SAFETY: negative infinity is defined to map to zero copies.
    let zero_from_log10 = zero_from_log10.unwrap();

    assert_eq!(zero_from_log2.count(), Count::new(0));
    assert_eq!(zero_from_log2.ploidy(), Ploidy::DIPLOID);
    assert_eq!(zero_from_log2.change(), Change::Loss);
    assert_eq!(zero_from_log10.count(), Count::new(0));
    assert_eq!(zero_from_log10.ploidy(), Ploidy::DIPLOID);
    assert_eq!(zero_from_log10.change(), Change::Loss);
}

#[test]
fn rejects_non_integral_logarithmic_variant_counts() {
    let from_log2 = Variant::try_from_log2("seq0", 100, 200, (1.25_f64).log2(), Ploidy::DIPLOID);
    let from_log10 = Variant::try_from_log10("seq0", 100, 200, (1.25_f64).log10(), Ploidy::DIPLOID);

    assert_eq!(
        from_log2,
        Err(LogarithmicVariantError::Logarithmic(
            LogarithmicError::NonIntegral,
        ))
    );
    assert_eq!(
        from_log10,
        Err(LogarithmicVariantError::Logarithmic(
            LogarithmicError::NonIntegral,
        ))
    );
}

#[test]
fn rejects_overflow_and_underflow_in_logarithmic_variant_construction() {
    let overflow = Variant::try_from_log2(
        "seq0",
        100,
        200,
        ((f64::from(u32::MAX) + 1.0) / 2.0).log2(),
        Ploidy::DIPLOID,
    );
    let underflow =
        Variant::try_from_log10("seq0", 100, 200, (5e-301_f64).log10(), Ploidy::DIPLOID);

    assert_eq!(
        overflow,
        Err(LogarithmicVariantError::Logarithmic(
            LogarithmicError::Overflow,
        ))
    );
    assert_eq!(
        underflow,
        Err(LogarithmicVariantError::Logarithmic(
            LogarithmicError::Underflow,
        ))
    );
}

#[test]
fn rejects_invalid_contigs_in_logarithmic_variant_construction() {
    let err = Variant::try_from_log2("", 100, 200, (1.5_f64).log2(), Ploidy::DIPLOID);

    assert!(matches!(
        err,
        Err(LogarithmicVariantError::Variant(Error::Contig(_)))
    ));
}

#[test]
fn reference_ploidy_participates_in_identity() {
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // diploid ploidy.
    let diploid = Variant::try_new("seq0", 100, 200, 3, Ploidy::DIPLOID).unwrap();
    // SAFETY: `3` is a valid nonzero ploidy.
    let triploid_ploidy = Ploidy::try_new(3).unwrap();
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // triploid ploidy.
    let triploid = Variant::try_new("seq0", 100, 200, 3, triploid_ploidy).unwrap();

    assert_ne!(diploid, triploid);
}

#[test]
fn round_trips_a_non_diploid_variant_through_the_canonical_form() {
    // SAFETY: `4` is a valid nonzero ploidy.
    let tetraploid = Ploidy::try_new(4).unwrap();
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // tetraploid ploidy.
    let variant = Variant::try_new("seq0", 100, 200, 3, tetraploid).unwrap();
    let rendered = variant.to_string();

    assert_eq!(rendered, "seq0:100-200(i):3/4");
    // SAFETY: parsing a displayed valid variant should reconstruct the original
    // value.
    let reparsed = rendered.parse::<Variant>().unwrap();

    assert_eq!(reparsed, variant);
}

#[test]
fn classifies_loss_reference_gain_and_complete_loss() {
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // diploid ploidy.
    let loss = Variant::try_new("seq0", 100, 200, 1, Ploidy::DIPLOID).unwrap();
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // diploid ploidy.
    let reference = Variant::try_new("seq0", 100, 200, 2, Ploidy::DIPLOID).unwrap();
    // SAFETY: the constructor is expected to accept this valid contig, span, and
    // diploid ploidy.
    let gain = Variant::try_new("seq0", 100, 200, 3, Ploidy::DIPLOID).unwrap();
    // SAFETY: the constructor is expected to accept a complete copy loss with a
    // valid diploid ploidy.
    let complete_loss = Variant::try_new("seq0", 100, 200, 0, Ploidy::DIPLOID).unwrap();

    assert_eq!(loss.change(), Change::Loss);
    assert_eq!(reference.change(), Change::Reference);
    assert_eq!(gain.change(), Change::Gain);
    assert_eq!(complete_loss.count(), Count::new(0));
    assert_eq!(complete_loss.change(), Change::Loss);
}

#[test]
fn rejects_empty_and_reversed_regions() {
    assert!(matches!(
        Variant::try_new("seq0", 100, 100, 2, Ploidy::DIPLOID),
        Err(Error::EmptyRegion)
    ));
    assert!(matches!(
        Variant::try_new("seq0", 200, 100, 2, Ploidy::DIPLOID),
        Err(Error::ReversedRegion)
    ));
}

#[test]
fn parses_and_displays_the_canonical_form() {
    // SAFETY: the parser is expected to accept the canonical copy-number form.
    let variant = "seq0:100-200(i):3/2".parse::<Variant>().unwrap();
    assert_eq!(variant.ploidy(), Ploidy::DIPLOID);
    assert_eq!(variant.to_string(), "seq0:100-200(i):3/2");
}

#[test]
fn round_trips_a_variant_whose_contig_contains_separators() {
    // SAFETY: the constructor is expected to accept this non-empty contig, span,
    // and diploid ploidy.
    let variant = Variant::try_new("assembly:chr:1", 100, 200, 3, Ploidy::DIPLOID).unwrap();
    let rendered = variant.to_string();

    assert_eq!(rendered, "assembly:chr:1:100-200(i):3/2");
    // SAFETY: parsing a displayed valid variant should reconstruct the original
    // value.
    let reparsed = rendered.parse::<Variant>().unwrap();

    assert_eq!(reparsed, variant);
}

#[test]
fn rejects_a_missing_interbase_qualifier() {
    let err = "seq0:100-200:3/2".parse::<Variant>().unwrap_err();
    assert!(matches!(err, ParseError::Qualifier { .. }));
}

#[test]
fn rejects_a_malformed_position_range() {
    let err = "seq0:100200(i):3/2".parse::<Variant>().unwrap_err();
    assert!(matches!(err, ParseError::Range { .. }));
}

#[test]
fn rejects_a_non_integral_count() {
    let err = "seq0:100-200(i):3.5/2".parse::<Variant>().unwrap_err();
    assert!(matches!(err, ParseError::Copies { .. }));
}

#[test]
fn rejects_a_count_only_copy_number_string() {
    let err = "seq0:100-200(i):3".parse::<Variant>().unwrap_err();

    assert_eq!(err, ParseError::Format("seq0:100-200(i):3".to_string()));
}

#[test]
fn rejects_an_invalid_ploidy() {
    let err = "seq0:100-200(i):3/not-a-number"
        .parse::<Variant>()
        .unwrap_err();

    assert!(matches!(err, ParseError::Ploidy { .. }));
}

#[test]
fn rejects_a_zero_ploidy() {
    let err = "seq0:100-200(i):3/0".parse::<Variant>().unwrap_err();

    assert_eq!(err, ParseError::ZeroPloidy);
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
fn final_review_accepts_a_maximum_count_log2_round_trip_within_tolerance() {
    let max_count = Count::new(u32::MAX);
    let log2 = max_count.log2(Ploidy::DIPLOID);

    // SAFETY: a logarithmic round trip from a valid maximum count must reconstruct
    // the same count.
    let round_trip = Count::try_from_log2(log2, Ploidy::DIPLOID).unwrap();

    assert_eq!(round_trip, max_count);
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

#[test]
fn rejects_finite_inputs_that_round_to_zero() {
    let tiny_log2 = (5e-301_f64).log2();

    assert_eq!(
        Count::try_from_log2(tiny_log2, Ploidy::DIPLOID),
        Err(LogarithmicError::Underflow)
    );
}

#[test]
fn final_review_maps_an_invalid_contig_to_the_dedicated_parse_error() {
    let err = ":100-200(i):3/2".parse::<Variant>().unwrap_err();

    assert!(matches!(err, ParseError::Contig(_)));
}

#[test]
fn final_review_maps_an_empty_region_to_the_dedicated_parse_error() {
    let err = "seq0:100-100(i):3/2".parse::<Variant>().unwrap_err();

    assert_eq!(err, ParseError::EmptyRegion);
}

#[test]
fn final_review_maps_a_reversed_region_to_the_dedicated_parse_error() {
    let err = "seq0:200-100(i):3/2".parse::<Variant>().unwrap_err();

    assert_eq!(err, ParseError::ReversedRegion);
}
