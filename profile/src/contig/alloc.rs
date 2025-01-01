//! Profiling of contig allocations.

use std::hint::black_box;

use clap::Parser;
use omics_coordinate::Contig;

/// Profiling of contig allocations.
#[derive(Parser, Debug)]
pub struct Args {
    /// The number of allocations to use.
    #[arg(default_value_t = 5_000_000_000)]
    n: usize,
}

/// The main method.
pub fn main(args: Args) {
    for _ in 0..args.n {
        let _ = black_box(Contig::new("seq0"));
    }
}
