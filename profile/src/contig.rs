//! Profiling of contigs.

pub mod alloc;

use clap::Parser;
use clap::Subcommand;

/// The different aspects of contigs that can be profiled.
#[derive(Subcommand)]
pub enum Command {
    /// Allocations of contigs.
    Alloc(alloc::Args),
}

/// Profiling of contigs.
#[derive(Parser)]
pub struct Args {
    /// The subcommand.
    #[command(subcommand)]
    command: Command,
}

/// The main method.
pub fn main(args: Args) {
    match args.command {
        Command::Alloc(args) => alloc::main(args),
    }
}
