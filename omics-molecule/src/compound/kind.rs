//! Kinds of compounds for nucleotides.

/// A kind of compound for a nucleotide.
#[derive(Debug, Eq, PartialEq)]
pub enum Kind {
    /// Purine nucleotide.
    Purine,

    /// Pyrimidine nucleotide.
    Pyrimidine,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compound::Nucleotide;
    use crate::polymer::dna;
    use crate::polymer::rna;

    #[test]
    fn it_correctly_assigns_a_type_to_dna_nucleotides() {
        assert_eq!(dna::Nucleotide::A.kind(), Kind::Purine);
        assert_eq!(dna::Nucleotide::G.kind(), Kind::Purine);
        assert_eq!(dna::Nucleotide::C.kind(), Kind::Pyrimidine);
        assert_eq!(dna::Nucleotide::T.kind(), Kind::Pyrimidine);
    }

    #[test]
    fn it_correctly_assigns_a_type_to_rna_nucleotides() {
        assert_eq!(rna::Nucleotide::A.kind(), Kind::Purine);
        assert_eq!(rna::Nucleotide::G.kind(), Kind::Purine);
        assert_eq!(rna::Nucleotide::C.kind(), Kind::Pyrimidine);
        assert_eq!(rna::Nucleotide::U.kind(), Kind::Pyrimidine);
    }
}
