//! Nucleotides in DNA.

use thiserror::Error;

use crate::compound::Kind;
use crate::compound::nucleotide::Analogous;
use crate::compound::nucleotide::Transcribe;
use crate::polymer::rna;

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

/// A nucleotide in an DNA context.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Nucleotide {
    /// Adenine.
    A,

    /// Cytosine.
    C,

    /// Guanine.
    G,

    /// Thymine.
    T,
}

impl crate::compound::Nucleotide for Nucleotide {
    fn kind(&self) -> Kind {
        match self {
            Nucleotide::A => Kind::Purine,
            Nucleotide::C => Kind::Pyrimidine,
            Nucleotide::G => Kind::Purine,
            Nucleotide::T => Kind::Pyrimidine,
        }
    }
}

impl Analogous<rna::Nucleotide> for Nucleotide {
    fn analogous(&self) -> rna::Nucleotide {
        match self {
            Nucleotide::A => rna::Nucleotide::A,
            Nucleotide::C => rna::Nucleotide::C,
            Nucleotide::G => rna::Nucleotide::G,
            Nucleotide::T => rna::Nucleotide::U,
        }
    }
}

impl Transcribe<rna::Nucleotide> for Nucleotide {
    fn transcribe(&self) -> rna::Nucleotide {
        match self {
            Nucleotide::A => rna::Nucleotide::U,
            Nucleotide::C => rna::Nucleotide::G,
            Nucleotide::G => rna::Nucleotide::C,
            Nucleotide::T => rna::Nucleotide::A,
        }
    }
}

impl std::fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Nucleotide::A => write!(f, "A"),
            Nucleotide::C => write!(f, "C"),
            Nucleotide::G => write!(f, "G"),
            Nucleotide::T => write!(f, "T"),
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
            'T' | 't' => Ok(Nucleotide::T),
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
            'T' | 't' => Ok(Nucleotide::T),
            _ => Err(Error::ParseError(ParseError::InvalidNucleotide(c))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_correctly_defines_analogous_nucleotides() {
        assert_eq!(Nucleotide::A.analogous(), rna::Nucleotide::A);
        assert_eq!(Nucleotide::C.analogous(), rna::Nucleotide::C);
        assert_eq!(Nucleotide::G.analogous(), rna::Nucleotide::G);
        assert_eq!(Nucleotide::T.analogous(), rna::Nucleotide::U);
    }

    #[test]
    fn it_correctly_defines_transcribed_nucleotides() {
        assert_eq!(Nucleotide::A.transcribe(), rna::Nucleotide::U);
        assert_eq!(Nucleotide::C.transcribe(), rna::Nucleotide::G);
        assert_eq!(Nucleotide::G.transcribe(), rna::Nucleotide::C);
        assert_eq!(Nucleotide::T.transcribe(), rna::Nucleotide::A);
    }

    #[test]
    fn it_correctly_serializes_nucleotides() {
        assert_eq!(Nucleotide::A.to_string(), "A");
        assert_eq!(Nucleotide::C.to_string(), "C");
        assert_eq!(Nucleotide::G.to_string(), "G");
        assert_eq!(Nucleotide::T.to_string(), "T");
    }

    #[test]
    fn it_correctly_creates_a_nucleotide_from_a_char() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(Nucleotide::try_from('a')?, Nucleotide::A);
        assert_eq!(Nucleotide::try_from('A')?, Nucleotide::A);
        assert_eq!(Nucleotide::try_from('c')?, Nucleotide::C);
        assert_eq!(Nucleotide::try_from('C')?, Nucleotide::C);
        assert_eq!(Nucleotide::try_from('g')?, Nucleotide::G);
        assert_eq!(Nucleotide::try_from('G')?, Nucleotide::G);
        assert_eq!(Nucleotide::try_from('t')?, Nucleotide::T);
        assert_eq!(Nucleotide::try_from('T')?, Nucleotide::T);

        let err = Nucleotide::try_from('u').unwrap_err();
        assert_eq!(err.to_string(), "invalid nucleotide `u`");

        let err = Nucleotide::try_from('U').unwrap_err();
        assert_eq!(err.to_string(), "invalid nucleotide `U`");

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
        assert_eq!("t".parse::<Nucleotide>()?, Nucleotide::T);
        assert_eq!("T".parse::<Nucleotide>()?, Nucleotide::T);

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

        let err = "u".parse::<Nucleotide>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidNucleotide(_))
        ));
        assert_eq!(err.to_string(), "invalid nucleotide `u`");

        let err = "U".parse::<Nucleotide>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidNucleotide(_))
        ));
        assert_eq!(err.to_string(), "invalid nucleotide `U`");

        Ok(())
    }
}
