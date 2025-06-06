//! Nucleotides in RNA.

use thiserror::Error;

use crate::compound::Kind;
use crate::compound::nucleotide::Analogous;
use crate::compound::nucleotide::ReverseTranscribe;
use crate::polymer::dna;

// NOTE: there are two invalid nucleotide errors—one under [`ParseError`] and
// one under [`Error`]. I felt this was appropriate to represent the following
// two cases respectively:
//
// * A case where a nucleotide is failed to be parsed from a [`str`].
// * A case where a nucleotide is failed to be converted from a [`char`]
//   directly.
//
// Just wanted to point out to any future folks that this was intentional though
// the signatures look identical (they will be converted to string differently).

/// An error when parsing a nucleotide.
#[derive(Error, Debug)]
pub enum ParseError {
    /// An invalid format was attempted to be parsed.
    #[error("invalid nucleotide format `{0}`")]
    InvalidFormat(String),

    /// An invalid nucleotide was attempted to be parsed.
    #[error("invalid nucleotide `{0}`")]
    InvalidNucleotide(char),
}

/// An error related to a [`Nucleotide`].
#[derive(Error, Debug)]
pub enum Error {
    /// An invalid nucleotide was attempted to be created from a [`char`].
    #[error("invalid nucleotide `{0}`")]
    InvalidNucleotide(char),

    /// A parse error.
    #[error(transparent)]
    ParseError(#[from] ParseError),
}

/// A nucleotide in an RNA context.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Nucleotide {
    /// Adenine.
    A,

    /// Cytosine.
    C,

    /// Guanine.
    G,

    /// Uracil.
    U,
}

impl crate::compound::Nucleotide for Nucleotide {
    fn kind(&self) -> Kind {
        match self {
            Nucleotide::A => Kind::Purine,
            Nucleotide::C => Kind::Pyrimidine,
            Nucleotide::G => Kind::Purine,
            Nucleotide::U => Kind::Pyrimidine,
        }
    }
}

impl Analogous<dna::Nucleotide> for Nucleotide {
    fn analogous(&self) -> dna::Nucleotide {
        match self {
            Nucleotide::A => dna::Nucleotide::A,
            Nucleotide::C => dna::Nucleotide::C,
            Nucleotide::G => dna::Nucleotide::G,
            Nucleotide::U => dna::Nucleotide::T,
        }
    }
}

impl ReverseTranscribe<dna::Nucleotide> for Nucleotide {
    fn reverse_transcribe(&self) -> dna::Nucleotide {
        match self {
            Nucleotide::A => dna::Nucleotide::T,
            Nucleotide::C => dna::Nucleotide::G,
            Nucleotide::G => dna::Nucleotide::C,
            Nucleotide::U => dna::Nucleotide::A,
        }
    }
}

impl std::fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Nucleotide::A => write!(f, "A"),
            Nucleotide::C => write!(f, "C"),
            Nucleotide::G => write!(f, "G"),
            Nucleotide::U => write!(f, "U"),
        }
    }
}

impl TryFrom<char> for Nucleotide {
    type Error = Error;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' | 'a' => Ok(Nucleotide::A),
            'C' | 'c' => Ok(Nucleotide::C),
            'G' | 'g' => Ok(Nucleotide::G),
            'U' | 'u' => Ok(Nucleotide::U),
            _ => Err(Error::InvalidNucleotide(c)),
        }
    }
}

impl std::str::FromStr for Nucleotide {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != 1 {
            return Err(Error::ParseError(ParseError::InvalidFormat(s.to_string())));
        }

        // SAFETY: we just ensured that the length is one, so this must unwrap.
        let c = s.chars().next().unwrap();

        match c {
            'A' | 'a' => Ok(Nucleotide::A),
            'C' | 'c' => Ok(Nucleotide::C),
            'G' | 'g' => Ok(Nucleotide::G),
            'U' | 'u' => Ok(Nucleotide::U),
            _ => Err(Error::ParseError(ParseError::InvalidNucleotide(c))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_correctly_defines_analogous_nucleotides() {
        assert_eq!(Nucleotide::A.analogous(), dna::Nucleotide::A);
        assert_eq!(Nucleotide::C.analogous(), dna::Nucleotide::C);
        assert_eq!(Nucleotide::G.analogous(), dna::Nucleotide::G);
        assert_eq!(Nucleotide::U.analogous(), dna::Nucleotide::T);
    }

    #[test]
    fn it_correctly_defines_reverse_transcribed_nucleotides() {
        assert_eq!(Nucleotide::A.reverse_transcribe(), dna::Nucleotide::T);
        assert_eq!(Nucleotide::C.reverse_transcribe(), dna::Nucleotide::G);
        assert_eq!(Nucleotide::G.reverse_transcribe(), dna::Nucleotide::C);
        assert_eq!(Nucleotide::U.reverse_transcribe(), dna::Nucleotide::A);
    }

    #[test]
    fn it_correctly_serializes_nucleotides() {
        assert_eq!(Nucleotide::A.to_string(), "A");
        assert_eq!(Nucleotide::C.to_string(), "C");
        assert_eq!(Nucleotide::G.to_string(), "G");
        assert_eq!(Nucleotide::U.to_string(), "U");
    }

    #[test]
    fn it_correctly_creates_a_nucleotide_from_a_char() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(Nucleotide::try_from('a')?, Nucleotide::A);
        assert_eq!(Nucleotide::try_from('A')?, Nucleotide::A);
        assert_eq!(Nucleotide::try_from('c')?, Nucleotide::C);
        assert_eq!(Nucleotide::try_from('C')?, Nucleotide::C);
        assert_eq!(Nucleotide::try_from('g')?, Nucleotide::G);
        assert_eq!(Nucleotide::try_from('G')?, Nucleotide::G);
        assert_eq!(Nucleotide::try_from('u')?, Nucleotide::U);
        assert_eq!(Nucleotide::try_from('U')?, Nucleotide::U);

        let err = Nucleotide::try_from('t').unwrap_err();
        assert_eq!(err.to_string(), "invalid nucleotide `t`");

        let err = Nucleotide::try_from('T').unwrap_err();
        assert_eq!(err.to_string(), "invalid nucleotide `T`");

        let err = Nucleotide::try_from('q').unwrap_err();
        assert_eq!(err.to_string(), "invalid nucleotide `q`");

        Ok(())
    }

    #[test]
    fn it_correctly_deserializes_nucleotides() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!("a".parse::<Nucleotide>()?, Nucleotide::A);
        assert_eq!("A".parse::<Nucleotide>()?, Nucleotide::A);
        assert_eq!("c".parse::<Nucleotide>()?, Nucleotide::C);
        assert_eq!("C".parse::<Nucleotide>()?, Nucleotide::C);
        assert_eq!("g".parse::<Nucleotide>()?, Nucleotide::G);
        assert_eq!("G".parse::<Nucleotide>()?, Nucleotide::G);
        assert_eq!("u".parse::<Nucleotide>()?, Nucleotide::U);
        assert_eq!("U".parse::<Nucleotide>()?, Nucleotide::U);

        let err = "word".parse::<Nucleotide>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidFormat(_))
        ));
        assert_eq!(err.to_string(), "invalid nucleotide format `word`");

        let err = "q".parse::<Nucleotide>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidNucleotide(_))
        ));
        assert_eq!(err.to_string(), "invalid nucleotide `q`");

        let err = "t".parse::<Nucleotide>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidNucleotide(_))
        ));
        assert_eq!(err.to_string(), "invalid nucleotide `t`");

        let err = "T".parse::<Nucleotide>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidNucleotide(_))
        ));
        assert_eq!(err.to_string(), "invalid nucleotide `T`");

        Ok(())
    }
}
