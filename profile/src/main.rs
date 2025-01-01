//! A profiling tool for developing `omics`.

use clap::Parser;
use clap::Subcommand;

pub mod contig;

/// The different entities that can be profiled.
#[derive(Subcommand)]
pub enum Command {
    /// Contig profiling.
    Contig(contig::Args),
}

/// Profiling of `omics` entities.
#[derive(Parser)]
pub struct Args {
    /// The subcommand.
    #[command(subcommand)]
    command: Command,
}

/// The main method.
pub fn main() {
    let args = Args::parse();

    match args.command {
        Command::Contig(args) => contig::main(args),
    }
}
