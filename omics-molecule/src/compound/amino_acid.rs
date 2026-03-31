//! Amino acids.

/// A classification of an amino acid based on the biochemical properties of its
/// side chain.
///
/// This classification follows the standard Lehninger grouping of the twenty
/// canonical amino acids by side-chain chemistry.
#[derive(Debug, Eq, PartialEq)]
pub enum Classification {
    /// A nonpolar amino acid with an aliphatic side chain (Gly, Ala, Val, Leu,
    /// Ile, Pro, Met).
    NonpolarAliphatic,

    /// A nonpolar amino acid with an aromatic side chain (Phe, Trp, Tyr).
    Aromatic,

    /// A polar amino acid with an uncharged side chain at physiological pH
    /// (Ser, Thr, Cys, Asn, Gln).
    PolarUncharged,

    /// An amino acid with a positively charged side chain at physiological pH
    /// (Lys, Arg, His).
    PositivelyCharged,

    /// An amino acid with a negatively charged side chain at physiological pH
    /// (Asp, Glu).
    NegativelyCharged,
}

/// A marker trait that denotes a type of amino acid.
pub trait AminoAcid:
    std::fmt::Debug + std::fmt::Display + Copy + Eq + PartialEq + std::str::FromStr
{
    /// Gets the [`Classification`] for a given [`AminoAcid`].
    fn classification(&self) -> Classification;
}
