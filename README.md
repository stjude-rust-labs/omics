<p align="center">
  <h1 align="center">
   <code>omics</code>
  </h1>

  <p align="center">
    <a href="https://github.com/stjude-rust-labs/omics/actions/workflows/CI.yml" target="_blank">
      <img alt="CI: Status" src="https://github.com/stjude-rust-labs/omics/actions/workflows/CI.yml/badge.svg" />
    </a>
    <a href="https://crates.io/crates/omics" target="_blank">
      <img alt="crates.io version" src="https://img.shields.io/crates/v/omics">
    </a>
    <img alt="crates.io downloads" src="https://img.shields.io/crates/d/omics">
    <a href="https://github.com/stjude-rust-labs/omics/blob/main/LICENSE-APACHE" target="_blank">
      <img alt="License: Apache 2.0" src="https://img.shields.io/badge/license-Apache 2.0-blue.svg" />
    </a>
    <a href="https://github.com/stjude-rust-labs/omics/blob/main/LICENSE-MIT" target="_blank">
      <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg" />
    </a>
  </p>

  <p align="center">
    Foundations for the Rust omics ecosystem.
    <br />
    <a href="https://docs.rs/omics"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/stjude-rust-labs/omics/issues/new?assignees=&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    ·
    <a href="https://github.com/stjude-rust-labs/omics/issues/new?assignees=&title=Descriptive%20Title&labels=bug">Report Bug</a>
    ·
    ⭐ Consider starring the repo! ⭐
    <br />
  </p>
</p>

**NOTE: this crate is highly experimental and is not ready for production or even really 
experimenting with. There are many things I already don't like about or are arguably
up incorrect within the existing implementation. As I find time, I'll continue to update
the crate with these changes.**

## 🎨 Features

The `omics` crate provides foundational data structures for working with omics data. The
following features are included and are top-level goals of the project.

- **Coordinate systems (`omics::coordinate`)**. Provides multiple genomic coordinate
  systems including both a 0-based, half-open coordinate system (also known as the
  _interbase_ coordinate system) and 1-based, fully-closed coordinate system (also known
  as the _in-base_ or just _base_ coordinate system).
- **Biologically relevant molecules (`omics::molecule`).** Representations of molecules
  relevant in omics including smaller compounds (e.g., nucleotides) and larger
  polymers (e.g., DNA, RNA).
- **Variation (`omics::variation`).** Facilities for expressing common types of
  variation including single nucleotide variations (SNVs), insertions/deletions
  (INDELs), structural variations (SVs), and copy number variations (CNVs).

## 🖥️ Development

To bootstrap a development environment, please use the following commands.

```bash
# Clone the repository
git clone git@github.com:stjude-rust-labs/omics.git
cd omics

# Build the crate in release mode
cargo build --release

# List out the examples
cargo run --release --example
```

## 🚧️ Tests

Before submitting any pull requests, please make sure the code passes the following checks (from the
root directory).

```bash
# Run the project's tests.
cargo test --all-features

# Run the tests for the examples.
cargo test --examples --all-features

# Ensure the project doesn't have any linting warnings.
cargo clippy --all-features

# Ensure the project passes `cargo fmt`.
cargo fmt --check

# Ensure the docs build.
cargo doc
```

## 🤝 Contributing

Contributions, issues and feature requests are welcome! Feel free to check [issues
page](https://github.com/stjude-rust-labs/omics/issues).

## 📝 License

This project is licensed as either [Apache 2.0][license-apache] or
[MIT][license-mit] at your discretion. Additionally, please see [the
disclaimer](https://github.com/stjude-rust-labs#disclaimer) that applies to all
crates and command line tools made available by St. Jude Rust Labs.

Copyright © 2024-Present [St. Jude Children's Research Hospital](https://github.com/stjude).

[license-apache]: https://github.com/stjude-rust-labs/omics/blob/main/LICENSE-APACHE
[license-mit]: https://github.com/stjude-rust-labs/omics/blob/main/LICENSE-MIT
