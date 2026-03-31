//! Amino acids in proteins.

use thiserror::Error;

use crate::compound::amino_acid::Classification;

/// An error when parsing an amino acid.
#[derive(Error, Debug)]
pub enum ParseError {
    /// An invalid format was attempted to be parsed.
    #[error("invalid amino acid format `{0}`")]
    InvalidFormat(String),

    /// An invalid amino acid was attempted to be parsed.
    #[error("invalid amino acid `{0}`")]
    InvalidAminoAcid(char),
}

/// An error related to an [`AminoAcid`].
#[derive(Error, Debug)]
pub enum Error {
    /// An invalid amino acid was attempted to be created from a [`char`].
    #[error("invalid amino acid `{0}`")]
    InvalidAminoAcid(char),

    /// A parse error.
    #[error(transparent)]
    ParseError(#[from] ParseError),
}

/// An amino acid in a protein context.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum AminoAcid {
    /// Alanine.
    A,

    /// Arginine.
    R,

    /// Asparagine.
    N,

    /// Aspartic acid.
    D,

    /// Cysteine.
    C,

    /// Glutamic acid.
    E,

    /// Glutamine.
    Q,

    /// Glycine.
    G,

    /// Histidine.
    H,

    /// Isoleucine.
    I,

    /// Leucine.
    L,

    /// Lysine.
    K,

    /// Methionine.
    M,

    /// Phenylalanine.
    F,

    /// Proline.
    P,

    /// Serine.
    S,

    /// Threonine.
    T,

    /// Tryptophan.
    W,

    /// Tyrosine.
    Y,

    /// Valine.
    V,
}

impl crate::compound::AminoAcid for AminoAcid {
    fn classification(&self) -> Classification {
        match self {
            AminoAcid::G
            | AminoAcid::A
            | AminoAcid::V
            | AminoAcid::L
            | AminoAcid::I
            | AminoAcid::P
            | AminoAcid::M => Classification::NonpolarAliphatic,
            AminoAcid::F | AminoAcid::W | AminoAcid::Y => Classification::Aromatic,
            AminoAcid::S | AminoAcid::T | AminoAcid::C | AminoAcid::N | AminoAcid::Q => {
                Classification::PolarUncharged
            }
            AminoAcid::K | AminoAcid::R | AminoAcid::H => Classification::PositivelyCharged,
            AminoAcid::D | AminoAcid::E => Classification::NegativelyCharged,
        }
    }
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AminoAcid::A => write!(f, "A"),
            AminoAcid::R => write!(f, "R"),
            AminoAcid::N => write!(f, "N"),
            AminoAcid::D => write!(f, "D"),
            AminoAcid::C => write!(f, "C"),
            AminoAcid::E => write!(f, "E"),
            AminoAcid::Q => write!(f, "Q"),
            AminoAcid::G => write!(f, "G"),
            AminoAcid::H => write!(f, "H"),
            AminoAcid::I => write!(f, "I"),
            AminoAcid::L => write!(f, "L"),
            AminoAcid::K => write!(f, "K"),
            AminoAcid::M => write!(f, "M"),
            AminoAcid::F => write!(f, "F"),
            AminoAcid::P => write!(f, "P"),
            AminoAcid::S => write!(f, "S"),
            AminoAcid::T => write!(f, "T"),
            AminoAcid::W => write!(f, "W"),
            AminoAcid::Y => write!(f, "Y"),
            AminoAcid::V => write!(f, "V"),
        }
    }
}

impl TryFrom<char> for AminoAcid {
    type Error = Error;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' | 'a' => Ok(AminoAcid::A),
            'R' | 'r' => Ok(AminoAcid::R),
            'N' | 'n' => Ok(AminoAcid::N),
            'D' | 'd' => Ok(AminoAcid::D),
            'C' | 'c' => Ok(AminoAcid::C),
            'E' | 'e' => Ok(AminoAcid::E),
            'Q' | 'q' => Ok(AminoAcid::Q),
            'G' | 'g' => Ok(AminoAcid::G),
            'H' | 'h' => Ok(AminoAcid::H),
            'I' | 'i' => Ok(AminoAcid::I),
            'L' | 'l' => Ok(AminoAcid::L),
            'K' | 'k' => Ok(AminoAcid::K),
            'M' | 'm' => Ok(AminoAcid::M),
            'F' | 'f' => Ok(AminoAcid::F),
            'P' | 'p' => Ok(AminoAcid::P),
            'S' | 's' => Ok(AminoAcid::S),
            'T' | 't' => Ok(AminoAcid::T),
            'W' | 'w' => Ok(AminoAcid::W),
            'Y' | 'y' => Ok(AminoAcid::Y),
            'V' | 'v' => Ok(AminoAcid::V),
            _ => Err(Error::InvalidAminoAcid(c)),
        }
    }
}

impl std::str::FromStr for AminoAcid {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != 1 {
            return Err(Error::ParseError(ParseError::InvalidFormat(s.to_string())));
        }

        // SAFETY: we just ensured that the length is one, so this must unwrap.
        let c = s.chars().next().unwrap();

        match c {
            'A' | 'a' => Ok(AminoAcid::A),
            'R' | 'r' => Ok(AminoAcid::R),
            'N' | 'n' => Ok(AminoAcid::N),
            'D' | 'd' => Ok(AminoAcid::D),
            'C' | 'c' => Ok(AminoAcid::C),
            'E' | 'e' => Ok(AminoAcid::E),
            'Q' | 'q' => Ok(AminoAcid::Q),
            'G' | 'g' => Ok(AminoAcid::G),
            'H' | 'h' => Ok(AminoAcid::H),
            'I' | 'i' => Ok(AminoAcid::I),
            'L' | 'l' => Ok(AminoAcid::L),
            'K' | 'k' => Ok(AminoAcid::K),
            'M' | 'm' => Ok(AminoAcid::M),
            'F' | 'f' => Ok(AminoAcid::F),
            'P' | 'p' => Ok(AminoAcid::P),
            'S' | 's' => Ok(AminoAcid::S),
            'T' | 't' => Ok(AminoAcid::T),
            'W' | 'w' => Ok(AminoAcid::W),
            'Y' | 'y' => Ok(AminoAcid::Y),
            'V' | 'v' => Ok(AminoAcid::V),
            _ => Err(Error::ParseError(ParseError::InvalidAminoAcid(c))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_correctly_classifies_nonpolar_aliphatic_amino_acids() {
        use crate::compound::AminoAcid as _;

        assert_eq!(
            AminoAcid::G.classification(),
            Classification::NonpolarAliphatic
        );
        assert_eq!(
            AminoAcid::A.classification(),
            Classification::NonpolarAliphatic
        );
        assert_eq!(
            AminoAcid::V.classification(),
            Classification::NonpolarAliphatic
        );
        assert_eq!(
            AminoAcid::L.classification(),
            Classification::NonpolarAliphatic
        );
        assert_eq!(
            AminoAcid::I.classification(),
            Classification::NonpolarAliphatic
        );
        assert_eq!(
            AminoAcid::P.classification(),
            Classification::NonpolarAliphatic
        );
        assert_eq!(
            AminoAcid::M.classification(),
            Classification::NonpolarAliphatic
        );
    }

    #[test]
    fn it_correctly_classifies_aromatic_amino_acids() {
        use crate::compound::AminoAcid as _;

        assert_eq!(AminoAcid::F.classification(), Classification::Aromatic);
        assert_eq!(AminoAcid::W.classification(), Classification::Aromatic);
        assert_eq!(AminoAcid::Y.classification(), Classification::Aromatic);
    }

    #[test]
    fn it_correctly_classifies_polar_uncharged_amino_acids() {
        use crate::compound::AminoAcid as _;

        assert_eq!(
            AminoAcid::S.classification(),
            Classification::PolarUncharged
        );
        assert_eq!(
            AminoAcid::T.classification(),
            Classification::PolarUncharged
        );
        assert_eq!(
            AminoAcid::C.classification(),
            Classification::PolarUncharged
        );
        assert_eq!(
            AminoAcid::N.classification(),
            Classification::PolarUncharged
        );
        assert_eq!(
            AminoAcid::Q.classification(),
            Classification::PolarUncharged
        );
    }

    #[test]
    fn it_correctly_classifies_positively_charged_amino_acids() {
        use crate::compound::AminoAcid as _;

        assert_eq!(
            AminoAcid::K.classification(),
            Classification::PositivelyCharged
        );
        assert_eq!(
            AminoAcid::R.classification(),
            Classification::PositivelyCharged
        );
        assert_eq!(
            AminoAcid::H.classification(),
            Classification::PositivelyCharged
        );
    }

    #[test]
    fn it_correctly_classifies_negatively_charged_amino_acids() {
        use crate::compound::AminoAcid as _;

        assert_eq!(
            AminoAcid::D.classification(),
            Classification::NegativelyCharged
        );
        assert_eq!(
            AminoAcid::E.classification(),
            Classification::NegativelyCharged
        );
    }

    #[test]
    fn it_correctly_serializes_amino_acids() {
        assert_eq!(AminoAcid::A.to_string(), "A");
        assert_eq!(AminoAcid::R.to_string(), "R");
        assert_eq!(AminoAcid::N.to_string(), "N");
        assert_eq!(AminoAcid::D.to_string(), "D");
        assert_eq!(AminoAcid::C.to_string(), "C");
        assert_eq!(AminoAcid::E.to_string(), "E");
        assert_eq!(AminoAcid::Q.to_string(), "Q");
        assert_eq!(AminoAcid::G.to_string(), "G");
        assert_eq!(AminoAcid::H.to_string(), "H");
        assert_eq!(AminoAcid::I.to_string(), "I");
        assert_eq!(AminoAcid::L.to_string(), "L");
        assert_eq!(AminoAcid::K.to_string(), "K");
        assert_eq!(AminoAcid::M.to_string(), "M");
        assert_eq!(AminoAcid::F.to_string(), "F");
        assert_eq!(AminoAcid::P.to_string(), "P");
        assert_eq!(AminoAcid::S.to_string(), "S");
        assert_eq!(AminoAcid::T.to_string(), "T");
        assert_eq!(AminoAcid::W.to_string(), "W");
        assert_eq!(AminoAcid::Y.to_string(), "Y");
        assert_eq!(AminoAcid::V.to_string(), "V");
    }

    #[test]
    fn it_correctly_creates_an_amino_acid_from_a_char() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(AminoAcid::try_from('A')?, AminoAcid::A);
        assert_eq!(AminoAcid::try_from('a')?, AminoAcid::A);
        assert_eq!(AminoAcid::try_from('R')?, AminoAcid::R);
        assert_eq!(AminoAcid::try_from('r')?, AminoAcid::R);
        assert_eq!(AminoAcid::try_from('N')?, AminoAcid::N);
        assert_eq!(AminoAcid::try_from('n')?, AminoAcid::N);
        assert_eq!(AminoAcid::try_from('D')?, AminoAcid::D);
        assert_eq!(AminoAcid::try_from('d')?, AminoAcid::D);
        assert_eq!(AminoAcid::try_from('C')?, AminoAcid::C);
        assert_eq!(AminoAcid::try_from('c')?, AminoAcid::C);
        assert_eq!(AminoAcid::try_from('E')?, AminoAcid::E);
        assert_eq!(AminoAcid::try_from('e')?, AminoAcid::E);
        assert_eq!(AminoAcid::try_from('Q')?, AminoAcid::Q);
        assert_eq!(AminoAcid::try_from('q')?, AminoAcid::Q);
        assert_eq!(AminoAcid::try_from('G')?, AminoAcid::G);
        assert_eq!(AminoAcid::try_from('g')?, AminoAcid::G);
        assert_eq!(AminoAcid::try_from('H')?, AminoAcid::H);
        assert_eq!(AminoAcid::try_from('h')?, AminoAcid::H);
        assert_eq!(AminoAcid::try_from('I')?, AminoAcid::I);
        assert_eq!(AminoAcid::try_from('i')?, AminoAcid::I);
        assert_eq!(AminoAcid::try_from('L')?, AminoAcid::L);
        assert_eq!(AminoAcid::try_from('l')?, AminoAcid::L);
        assert_eq!(AminoAcid::try_from('K')?, AminoAcid::K);
        assert_eq!(AminoAcid::try_from('k')?, AminoAcid::K);
        assert_eq!(AminoAcid::try_from('M')?, AminoAcid::M);
        assert_eq!(AminoAcid::try_from('m')?, AminoAcid::M);
        assert_eq!(AminoAcid::try_from('F')?, AminoAcid::F);
        assert_eq!(AminoAcid::try_from('f')?, AminoAcid::F);
        assert_eq!(AminoAcid::try_from('P')?, AminoAcid::P);
        assert_eq!(AminoAcid::try_from('p')?, AminoAcid::P);
        assert_eq!(AminoAcid::try_from('S')?, AminoAcid::S);
        assert_eq!(AminoAcid::try_from('s')?, AminoAcid::S);
        assert_eq!(AminoAcid::try_from('T')?, AminoAcid::T);
        assert_eq!(AminoAcid::try_from('t')?, AminoAcid::T);
        assert_eq!(AminoAcid::try_from('W')?, AminoAcid::W);
        assert_eq!(AminoAcid::try_from('w')?, AminoAcid::W);
        assert_eq!(AminoAcid::try_from('Y')?, AminoAcid::Y);
        assert_eq!(AminoAcid::try_from('y')?, AminoAcid::Y);
        assert_eq!(AminoAcid::try_from('V')?, AminoAcid::V);
        assert_eq!(AminoAcid::try_from('v')?, AminoAcid::V);

        let err = AminoAcid::try_from('Z').unwrap_err();
        assert_eq!(err.to_string(), "invalid amino acid `Z`");

        let err = AminoAcid::try_from('1').unwrap_err();
        assert_eq!(err.to_string(), "invalid amino acid `1`");

        Ok(())
    }

    #[test]
    fn it_correctly_deserializes_amino_acids() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!("A".parse::<AminoAcid>()?, AminoAcid::A);
        assert_eq!("a".parse::<AminoAcid>()?, AminoAcid::A);
        assert_eq!("R".parse::<AminoAcid>()?, AminoAcid::R);
        assert_eq!("r".parse::<AminoAcid>()?, AminoAcid::R);
        assert_eq!("W".parse::<AminoAcid>()?, AminoAcid::W);
        assert_eq!("w".parse::<AminoAcid>()?, AminoAcid::W);
        assert_eq!("V".parse::<AminoAcid>()?, AminoAcid::V);
        assert_eq!("v".parse::<AminoAcid>()?, AminoAcid::V);

        let err = "word".parse::<AminoAcid>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidFormat(_))
        ));
        assert_eq!(err.to_string(), "invalid amino acid format `word`");

        let err = "Z".parse::<AminoAcid>().unwrap_err();
        assert!(matches!(
            err,
            Error::ParseError(ParseError::InvalidAminoAcid(_))
        ));
        assert_eq!(err.to_string(), "invalid amino acid `Z`");

        Ok(())
    }
}
