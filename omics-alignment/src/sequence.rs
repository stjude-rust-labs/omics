//! Encoded sequence representations for alignment.
//!
//! [`EncodedSequence`] provides a compact integer-encoded representation of
//! biological sequences suitable for use with scoring matrices. Conversions
//! from `omics-molecule` DNA, RNA, and protein types are provided.

use omics_molecule::polymer::dna;
use omics_molecule::polymer::protein;
use omics_molecule::polymer::rna;

/// A sequence of residues encoded as compact integer indices.
///
/// Each residue is represented as a `u8` index into a scoring matrix. The
/// encoding is determined by the source molecule type:
///
/// - **DNA**: A=0, C=1, G=2, T=3
/// - **RNA**: A=0, C=1, G=2, U=3
/// - **Protein**: BLOSUM order — A=0, R=1, N=2, D=3, C=4, Q=5, E=6, G=7, H=8,
///   I=9, L=10, K=11, M=12, F=13, P=14, S=15, T=16, W=17, Y=18, V=19
///
/// # Examples
///
/// ```
/// use omics_alignment::sequence::EncodedSequence;
/// use omics_molecule::polymer::dna;
///
/// let mol = "ACGT".parse::<dna::Molecule>()?;
/// let seq = EncodedSequence::from(&mol);
/// assert_eq!(seq.len(), 4);
/// assert_eq!(seq.get(0), 0); // A
/// assert_eq!(seq.get(2), 2); // G
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct EncodedSequence {
    /// The encoded residue data.
    data: Vec<u8>,
}

impl EncodedSequence {
    /// Returns the length of the sequence.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns `true` if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns the encoded residue at the given index.
    pub fn get(&self, index: usize) -> u8 {
        self.data[index]
    }
}

impl From<&dna::Molecule> for EncodedSequence {
    fn from(mol: &dna::Molecule) -> Self {
        let data = mol
            .inner()
            .iter()
            .map(|n| match n {
                dna::Nucleotide::A => 0,
                dna::Nucleotide::C => 1,
                dna::Nucleotide::G => 2,
                dna::Nucleotide::T => 3,
            })
            .collect();
        Self { data }
    }
}

impl From<&rna::Molecule> for EncodedSequence {
    fn from(mol: &rna::Molecule) -> Self {
        let data = mol
            .inner()
            .iter()
            .map(|n| match n {
                rna::Nucleotide::A => 0,
                rna::Nucleotide::C => 1,
                rna::Nucleotide::G => 2,
                rna::Nucleotide::U => 3,
            })
            .collect();
        Self { data }
    }
}

impl From<&protein::Molecule> for EncodedSequence {
    fn from(mol: &protein::Molecule) -> Self {
        let data = mol
            .inner()
            .iter()
            .map(|aa| match aa {
                protein::AminoAcid::A => 0,
                protein::AminoAcid::R => 1,
                protein::AminoAcid::N => 2,
                protein::AminoAcid::D => 3,
                protein::AminoAcid::C => 4,
                protein::AminoAcid::Q => 5,
                protein::AminoAcid::E => 6,
                protein::AminoAcid::G => 7,
                protein::AminoAcid::H => 8,
                protein::AminoAcid::I => 9,
                protein::AminoAcid::L => 10,
                protein::AminoAcid::K => 11,
                protein::AminoAcid::M => 12,
                protein::AminoAcid::F => 13,
                protein::AminoAcid::P => 14,
                protein::AminoAcid::S => 15,
                protein::AminoAcid::T => 16,
                protein::AminoAcid::W => 17,
                protein::AminoAcid::Y => 18,
                protein::AminoAcid::V => 19,
            })
            .collect();
        Self { data }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_encodes_dna_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let mol = "ACGT".parse::<dna::Molecule>()?;
        let seq = EncodedSequence::from(&mol);
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), 0);
        assert_eq!(seq.get(1), 1);
        assert_eq!(seq.get(2), 2);
        assert_eq!(seq.get(3), 3);
        Ok(())
    }

    #[test]
    fn it_encodes_rna_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let mol = "ACGU".parse::<rna::Molecule>()?;
        let seq = EncodedSequence::from(&mol);
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(3), 3); // U
        Ok(())
    }

    #[test]
    fn it_encodes_protein_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let mol = "ARNDCQE".parse::<protein::Molecule>()?;
        let seq = EncodedSequence::from(&mol);
        assert_eq!(seq.len(), 7);
        assert_eq!(seq.get(0), 0); // A
        assert_eq!(seq.get(1), 1); // R
        assert_eq!(seq.get(5), 5); // Q (BLOSUM order)
        assert_eq!(seq.get(6), 6); // E (BLOSUM order)
        Ok(())
    }

    #[test]
    fn it_reports_empty_sequences() {
        let mol = dna::Molecule::from(vec![]);
        let seq = EncodedSequence::from(&mol);
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
    }
}
