# Release

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
