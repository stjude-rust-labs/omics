//! Attempts to move a 1-based coordinate forward by the specified magnitude.
//!
//! # Examples
//!
//! * A positive-stranded coordinate: `cargo run --example
//!   coordinate_move_forward chr1:+:1 10`
//! * A negative-stranded coordinate: `cargo run --example
//!   coordinate_move_forward chr1:-:10 10`

use std::env::args;

use anyhow::Context;
use anyhow::anyhow;
use omics_coordinate::Coordinate;
use omics_coordinate::system::One;

fn main() -> anyhow::Result<()> {
    let coordinate = args()
        .nth(1)
        .expect("missing coordinate")
        .parse::<Coordinate<One>>()
        .with_context(|| "parsing 1-based coordinate")?;
    let magnitude = args()
        .nth(2)
        .expect("missing magnitude")
        .parse::<usize>()
        .with_context(|| "parsing magnitude")?;

    let result = coordinate
        .move_forward(magnitude)?
        .ok_or(anyhow!("out of bounds"))?;

    println!("{:#}", result);

    Ok(())
}
