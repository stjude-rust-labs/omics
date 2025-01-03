# Release

The release process for the `omics` family of crates is intentionally disjoint
across the repository. When releasing any new packages:

- You should first release each of the component crates (all crates that are not
  `omics`) in a sequential fashion using the [component crates](#component-crates)
  section.
- Next (if desired), you should release the convenience crate (the `omics` crate)
  by following the [convenience crate](#convenience-crate) section.

Notably, updates to the files listed below should be grouped in a single commit
**per crate**, but updates to these files across crates should not be contained
within a single commit.

## Component Crates

**Note:** in this example, we will be using the `omics-core` crate. Please
substitute the name of the crate that you are working on.

**Note:** crates should be released sequentially based on the dependency tree.
For example, `omics-core` should be released before `omics-coordinate`, as
`omics-coordinate` depends on `omics-core`.

For every component crate that has changes:

- [ ] Update version in `Cargo.toml`.
- [ ] Update versions of any `omics-*` crate dependencies in `Cargo.toml`.

  - If any of the crates do _not_ have updated versions, be sure they also
    don't have changes to release. If they do, you may be releasing the crates
    out of order!

- [ ] Update versions of any `omics-*` crates that depend on this crate. Be sure
      that you also add an entry to each crates `CHANGELOG.md` in the unreleased
      section (you'll later release those changes when the downstream crates are
      released).

      ```text
      ### Crate Updates

      - `omics-core`: bumped to v0.2.0
        ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-core-v0.2.0))
      ```

- [ ] Update `CHANGELOG.md` with version and publication date.
  - To get the changes to the crate since the last release, you can use a
    command like the following:
    ```bash
    git log omics-core-v0.1.0..HEAD --oneline -- omics-core
    ```
- [ ] Run tests: `cargo test --all-features`.
- [ ] Run tests for examples: `cargo test --examples --all-features`.
- [ ] Run linting: `cargo clippy --all-features`.
- [ ] Run fmt: `cargo fmt --check`.
- [ ] Run doc: `cargo doc`.
- [ ] Stage changes: `git add Cargo.toml CHANGELOG.md`.
- [ ] Create git commit:
  ```
  git commit -m "release: bumps `omics-core` version to v0.1.0"
  ```
- [ ] Create git tag:
  ```
  git tag omics-core-v0.1.0
  ```
- [ ] Push release: `git push && git push --tags`.
- [ ] Publish the component crate: `cargo publish --all-features`.
- [ ] Go to the Releases page in Github, create a Release for this tag, and
      copy the notes from the `CHANGELOG.md` file.

## Convenience Crate

From the root directory:

- [ ] Update the version of the top-level crate in the root `Cargo.toml`.
  - **Note:** changes to the version number will be automatically reflected in
    `omics/Cargo.toml`, as the version there is specified as `version.workspace =
true`.
- [ ] Run tests: `cargo test --all-features`.
- [ ] Run tests for examples: `cargo test --examples --all-features`.
- [ ] Run linting: `cargo clippy --all-features`.
- [ ] Run fmt: `cargo fmt --check`.
- [ ] Run doc: `cargo doc`.
- [ ] Stage changes: `git add Cargo.toml`.
- [ ] Create git commit:

  ```
  git commit
  ```

  The commit message should have a body conforming to this style:

  ```
  release: bumps `omics` version to v0.1.0

  Component Crate Updates
  -----------------------

  * `omics-core`: introduced at v0.1.0 ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-core-v0.1.0))
  * `omics-fictitous`: bumped from v0.1.0 to v0.2.0 ([release](https://github.com/stjude-rust-labs/omics/releases/tag/omics-fictitous-v0.2.0))
  ```

- [ ] Create git tag: `git tag omics-v0.1.0`.
- [ ] Push release: `git push && git push --tags`.
- [ ] Publish the new crate: `cargo publish --all-features -p omics`.
- [ ] Go to the Releases page in Github, create a Release for this tag, and
      copy the body from the commit message that describes the package version
      updates.
  - Ensure that you change the heading style from
    ```
    Component Crate Updates
    -----------------------
    ```
    to
    ```
    ## Component Crate Updates
    ```
