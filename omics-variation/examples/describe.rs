//! Attempts to describe a 1-based variant.
//!
//! # Examples
//!
//! * A single-nucleotide variation: `cargo run --example variation_describe
//!   chr1:+:1:A:T`

use std::env;

use omics_coordinate::system::One;
use omics_molecule::polymer::dna;
use omics_variation::Variant;

fn main() -> anyhow::Result<()> {
    let variant = env::args().nth(1).expect("missing variant");

    let variant = variant.parse::<Variant<dna::Nucleotide, One>>()?;
    println!("{:#?}", variant);

    Ok(())
}
