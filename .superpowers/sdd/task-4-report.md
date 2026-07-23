# Task 4 Report: Apple Silicon NEON

Task 4 adds the macOS AArch64 NEON implementation for the private
anti-diagonal wavefront engine. `algorithm::simd::{global, local}` now select
the narrowest proven NEON lane width on Apple Silicon and retain the scalar
engine when no `i16` or `i32` proof exists.

## Interface reconciliation

The reviewed engine exposes the private generic `Kernel` contract rather than
a production `Backend` interface. The brief's `Backend::Neon` test entry point
therefore becomes a test-only selector that invokes the actual
`NeonI16` and `NeonI32` `Kernel` implementations. No production adapter or
new public API was added.

## RED evidence

I added the macOS AArch64-gated `neon_matches_scalar` differential entry point
before creating the architecture module or either native kernel.

```text
$ cargo test -p omics-alignment neon_matches_scalar -- --nocapture
error[E0433]: cannot find `aarch64` in `super`
  --> omics-alignment/src/algorithm/simd/wavefront.rs:1513:53
...
error: could not compile `omics-alignment` (lib test) due to 6 previous errors
```

After adding the direct dispatch test, the broader native selector also failed
for the expected missing native module.

```text
$ cargo test -p omics-alignment neon_ -- --nocapture
error[E0433]: cannot find `aarch64` in `super`
error[E0433]: cannot find module or crate `aarch64` in this scope
error: could not compile `omics-alignment` (lib test) due to 8 previous errors
```

I also added a vector-fixture coverage test before changing the fixture. It
failed because the existing eight-lane case contained only mismatches along
the vector anti-diagonal.

```text
$ cargo test -p omics-alignment \
    test_vector_differential_cases_include_mixed_substitutions -- --nocapture
assertion failed: comparisons.clone().any(|(reference, query)| reference == query)
test ...::test_vector_differential_cases_include_mixed_substitutions ... FAILED
```

The fixture now contains both equal and unequal byte pairs in a complete
eight-lane vector, so the NEON substitution mask path exercises both arms.

## GREEN evidence

The final native differential run directly exercises `NeonI16` and
`NeonI32`, compares complete global and local outcomes against
`engine::{global, local}`, and checks the safe dispatch wrappers.

```text
$ cargo test -p omics-alignment neon_ -- --nocapture
running 2 tests
test algorithm::simd::tests::neon_dispatch_matches_scalar ... ok
test algorithm::simd::wavefront::tests::neon_matches_scalar ... ok
test result: ok. 2 passed; 0 failed
```

`NeonI16` receives exhaustive binary inputs through length five plus explicit
eight-lane vector, vector-plus-tail, and two-vector cases. `NeonI32` receives
the exhaustive sweep, which reaches its four-lane vector path and scalar tail.
The dispatch test covers eligible `i16` and `i32` score bounds.

```text
$ cargo test -p omics-alignment --all-features
unit tests:        58 passed; 0 failed
integration tests:  6 passed; 0 failed
doc tests:          4 passed; 0 failed
```

Final validation also passed.

```text
cargo +nightly fmt -- --check
cargo clippy --all-features -- --deny warnings
cargo check -p omics-alignment --all-features --target x86_64-apple-darwin
cargo test --all-features
cargo test --all-features --examples
cargo deny --workspace check
(cd omics && cargo msrv verify --output-format minimal --all-features)
```

`cargo fmt -- --check` exited successfully. Stable `rustfmt` emitted the
repository's existing notices about nightly-only formatting options.
`cargo doc` also exited successfully; it retains an unrelated existing
`omics-molecule` warning for the unresolved `Nucleotide(s)` link.

## Target-feature and memory safety

`aarch64.rs` only compiles under
`all(target_arch = "aarch64", target_os = "macos")`. Each pointer-accessing
native method carries `#[target_feature(enable = "neon")]`, and every unsafe
block documents the NEON guarantee and the `Kernel`-contract bounds. Apple
Silicon AArch64 exposes NEON as a baseline feature; the host reports
`target_feature="neon"` and the runtime feature probe returns `true`.

`NeonI16` loads and stores exactly eight `i16` lanes. Its substitution path
loads eight reference bytes, loads the eight bytes preceding `query_end`,
reverses the query lanes with `vrev64_u8`, compares with `vceq_u8`, sign-widens
the byte masks to eight all-bit `u16` masks, and selects the match or mismatch
score with `vbslq_s16`.

`NeonI32` loads and stores exactly four `i32` lanes. Its substitution path
first materializes four bounded reference and reverse-query bytes in local
`i32` staging arrays, then loads, compares, and selects with four `u32` masks.
The comparison results are explicitly bit-reinterpreted at the signed-vector
`Kernel` boundary and again for `vbslq_*`, preserving all-one and all-zero mask
semantics.

The dispatch wrappers use `wavefront::lane_width`, which shares the checked
full-path score bound with the generic execution guard. `LaneWidth::Scalar`
calls `engine::{global, local}` directly, and each generic wavefront entry
independently rechecks eligibility before invoking a native kernel.

## Files changed

- `omics-alignment/src/algorithm/simd/aarch64.rs` adds `NeonI16`,
  `NeonI32`, NEON intrinsics, and safe lane-width dispatch wrappers.
- `omics-alignment/src/algorithm/simd.rs` conditionally compiles the Apple
  Silicon backend and routes the public byte APIs to it.
- `omics-alignment/src/algorithm/simd/wavefront.rs` exposes the private
  lane-width selector to the sibling backend, scopes the former dead-code
  expectation to non-Apple targets, and adds forced-native differential and
  mixed-substitution coverage.

## Self-review and concerns

I reviewed signed comparison masks, query-lane order, vector load/store lane
bounds, the conservative score-range dispatch, target gating, and the scalar
fallback. A separate read-only review found no significant issues. There are
no Task 4 concerns. The stable-format and documentation warnings noted above
pre-date this change and do not affect the NEON implementation.
