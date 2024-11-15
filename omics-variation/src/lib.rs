//! Genomic variation.

use std::str::FromStr;

use omics_molecule::compound::Nucleotide;

pub mod snv;

/// An error related to a [`Variant`].
#[derive(Debug)]
pub enum Error {
    /// Unsuccessfully attempted to parse a [`Variant`] from a string.
    ParseError(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::ParseError(v) => write!(f, "unable to parse a variant from string: {v}"),
        }
    }
}

impl std::error::Error for Error {}

/// A variant.
#[derive(Debug)]
pub enum Variant<N: Nucleotide> {
    /// A single nucleotide substitution.
    SingleNucleotideVariation(snv::Variant<N>),
}

impl<N: Nucleotide> std::str::FromStr for Variant<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(snv) = s.parse::<snv::Variant<N>>() {
            return Ok(Variant::SingleNucleotideVariation(snv));
        }

        Err(Error::ParseError(s.to_string()))
    }
}

impl<N: Nucleotide> std::fmt::Display for Variant<N>
where
    <N as FromStr>::Err: std::fmt::Debug + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Variant::SingleNucleotideVariation(variant) => write!(f, "{}", variant),
        }
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::Strand;
    use omics_molecule::polymer::dna;

    use super::*;

    #[test]
    fn it_parses_dna_snvs_correctly() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:+:1:A:C".parse::<Variant<dna::Nucleotide>>()?;
        assert!(matches!(variant, Variant::SingleNucleotideVariation(_)));

        match variant {
            Variant::SingleNucleotideVariation(snv) => {
                assert_eq!(snv.coordinate().contig().as_str(), "seq0");
                assert_eq!(snv.coordinate().strand(), Strand::Positive);
                assert_eq!(snv.coordinate().position().get(), 1);
                assert_eq!(snv.reference(), &dna::Nucleotide::A);
                assert_eq!(snv.alternate(), &dna::Nucleotide::C);
            }
        }

        Ok(())
    }

    #[test]
    fn it_errors_when_attempting_to_parse_invalid_variants()
    -> Result<(), Box<dyn std::error::Error>> {
        let err = "seq0:1:A".parse::<Variant<dna::Nucleotide>>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "unable to parse a variant from string: seq0:1:A"
        );

        let err = "seq0:1:A:".parse::<Variant<dna::Nucleotide>>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "unable to parse a variant from string: seq0:1:A:"
        );

        let err = "seq0:A:C:1"
            .parse::<Variant<dna::Nucleotide>>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "unable to parse a variant from string: seq0:A:C:1"
        );

        Ok(())
    }

    #[test]
    fn it_serializes_correctly() -> Result<(), Box<dyn std::error::Error>> {
        let variant = "seq0:+:1:A:C".parse::<Variant<dna::Nucleotide>>()?;
        assert_eq!(variant.to_string(), "seq0:+:1:A:C");

        let variant = "seq0:+:1:A:C".parse::<Variant<dna::Nucleotide>>()?;
        assert_eq!(variant.to_string(), "seq0:+:1:A:C");

        Ok(())
    }
}
