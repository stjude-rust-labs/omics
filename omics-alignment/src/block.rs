//! Neutral aligned blocks pairing equal-length reference and query intervals.
//!
//! An [`AlignedBlock`] records a one-to-one alignment between a reference
//! interbase interval and a query interbase interval of the same entity count.
//! It carries no coordinate translation, liftover, or format-specific logic.

use omics_coordinate::Interval;
use omics_coordinate::position::Number;
use omics_coordinate::system::Interbase;
use thiserror::Error;

/// An interbase interval used in aligned-block positions.
type InterbaseInterval = Interval<Interbase>;

/// An error when constructing an [`AlignedBlock`].
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum Error {
    /// The reference and query intervals do not have the same entity count.
    #[error("reference entity count ({reference}) does not match query entity count ({query})")]
    EntityCountsDontMatch {
        /// The entity count of the reference interval.
        reference: Number,
        /// The entity count of the query interval.
        query: Number,
    },
}

/// A neutral aligned block pairing a reference interval with an equal-length
/// query interval.
///
/// Both intervals must span the same number of entities. No coordinate
/// translation or liftover is performed; this type is a plain data carrier.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignedBlock {
    /// The reference interbase interval.
    reference: InterbaseInterval,
    /// The query interbase interval.
    query: InterbaseInterval,
}

impl AlignedBlock {
    /// Constructs a new [`AlignedBlock`] from a reference and query interval.
    ///
    /// Returns [`Error::EntityCountsDontMatch`] when the two intervals span
    /// different numbers of entities.
    pub fn try_new(reference: InterbaseInterval, query: InterbaseInterval) -> Result<Self, Error> {
        let ref_count = reference.count_entities();
        let qry_count = query.count_entities();
        if ref_count != qry_count {
            return Err(Error::EntityCountsDontMatch {
                reference: ref_count,
                query: qry_count,
            });
        }
        Ok(Self { reference, query })
    }

    /// Returns a reference to the reference interval.
    pub fn reference(&self) -> &InterbaseInterval {
        &self.reference
    }

    /// Returns a reference to the query interval.
    pub fn query(&self) -> &InterbaseInterval {
        &self.query
    }

    /// Returns the number of entities spanned by this block.
    ///
    /// Because both intervals are required to be equal in length, this is
    /// equivalent to calling `count_entities()` on either interval.
    pub fn length(&self) -> Number {
        self.reference.count_entities()
    }

    /// Consumes the block and returns the reference interval.
    pub fn into_reference(self) -> InterbaseInterval {
        self.reference
    }

    /// Consumes the block and returns the query interval.
    pub fn into_query(self) -> InterbaseInterval {
        self.query
    }

    /// Consumes the block and returns the reference and query intervals.
    pub fn into_parts(self) -> (InterbaseInterval, InterbaseInterval) {
        (self.reference, self.query)
    }
}

#[cfg(test)]
mod tests {
    use omics_coordinate::interval::interbase::Interval;

    use super::*;

    #[test]
    fn equal_length_mixed_strands() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "ref:+:0-3".parse::<Interval>()?;
        let query = "query:-:5-2".parse::<Interval>()?;
        let block = AlignedBlock::try_new(reference.clone(), query.clone())?;
        assert_eq!(block.length(), 3);
        assert_eq!(block.reference(), &reference);
        assert_eq!(block.query(), &query);
        Ok(())
    }

    #[test]
    fn unequal_lengths_rejected() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "ref:+:0-3".parse::<Interval>()?;
        let query = "query:+:0-4".parse::<Interval>()?;
        assert!(matches!(
            AlignedBlock::try_new(reference, query),
            Err(Error::EntityCountsDontMatch {
                reference: 3,
                query: 4
            })
        ));
        Ok(())
    }

    #[test]
    fn equal_length_positive_strands() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "chr1:+:100-200".parse::<Interval>()?;
        let query = "read1:+:0-100".parse::<Interval>()?;
        let block = AlignedBlock::try_new(reference.clone(), query.clone())?;
        assert_eq!(block.length(), 100);
        assert_eq!(block.reference(), &reference);
        assert_eq!(block.query(), &query);
        Ok(())
    }

    #[test]
    fn into_reference_returns_original() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "ref:+:0-5".parse::<Interval>()?;
        let query = "query:+:10-15".parse::<Interval>()?;
        let block = AlignedBlock::try_new(reference.clone(), query)?;
        assert_eq!(block.into_reference(), reference);
        Ok(())
    }

    #[test]
    fn into_query_returns_original() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "ref:+:0-5".parse::<Interval>()?;
        let query = "query:+:10-15".parse::<Interval>()?;
        let block = AlignedBlock::try_new(reference, query.clone())?;
        assert_eq!(block.into_query(), query);
        Ok(())
    }

    #[test]
    fn into_parts_returns_both() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "ref:+:0-7".parse::<Interval>()?;
        let query = "query:-:20-13".parse::<Interval>()?;
        let block = AlignedBlock::try_new(reference.clone(), query.clone())?;
        assert_eq!(block.into_parts(), (reference, query));
        Ok(())
    }

    #[test]
    fn both_negative_strands() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "chr1:-:300-200".parse::<Interval>()?;
        let query = "read1:-:150-50".parse::<Interval>()?;
        let block = AlignedBlock::try_new(reference.clone(), query.clone())?;
        assert_eq!(block.length(), 100);
        Ok(())
    }
}
