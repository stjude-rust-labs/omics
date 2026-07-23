# Task 5 Report: Linux AVX2

Task 5 adds the private Linux x86_64 AVX2 backend for the byte-oriented SIMD
alignment API. `Avx2I16` processes sixteen `i16` lanes and `Avx2I32`
processes eight `i32` lanes. The public API selects AVX2 only after a runtime
feature check and otherwise calls the scalar engine.

## Interface reconciliation

The reviewed private `Kernel` contract supports the required AVX2 operations
without modification. The implementation therefore adds target-specific
`Kernel` implementations and safe dispatch wrappers rather than changing the
generic wavefront interface or public API.

## RED evidence

I added the x86_64 Linux-gated native selector and
`avx2_matches_scalar_when_available` before adding the x86_64 module. The
cross-target test compilation failed for the expected missing backend:

```text
$ rustup target add x86_64-unknown-linux-gnu && \
    cargo check -p omics-alignment --tests --target x86_64-unknown-linux-gnu
error[E0433]: cannot find `x86_64` in `super`
  --> omics-alignment/src/algorithm/simd/wavefront.rs:1537:53
...
error: could not compile `omics-alignment` (lib test) due to 8 previous errors
```

The initial implementation exposed an unfulfilled `dead_code` expectation on
the new target. `cargo check -p omics-alignment --tests --target
x86_64-unknown-linux-gnu` reproduced the warning. The expectation had excluded
the existing Apple backend but not the new Linux backend. I extended its target
condition, then reran the command without warnings.

Rust `1.81` also required explicit unsafe blocks for AVX2 intrinsics, while
the current compiler treats several of those calls as safe and warns about
unnecessary blocks. I kept the explicit, documented blocks for the MSRV and
scoped `#[allow(unused_unsafe)]` to each AVX2 implementation with a comment
explaining the compiler-version distinction. This keeps the current lint gate
and the MSRV compile check clean.

## GREEN evidence

The final cross-target checks succeeded:

```text
$ cargo check -p omics-alignment --all-features --tests \
    --target x86_64-unknown-linux-gnu
Finished `dev` profile ... target(s) in 0.57s

$ cargo clippy -p omics-alignment --all-features --lib \
    --target x86_64-unknown-linux-gnu -- --deny warnings
Finished `dev` profile ... target(s) in 0.20s

$ docker run --rm --platform linux/amd64 \
    --mount type=bind,source="$PWD",target=/workspace,readonly \
    --workdir /workspace \
    -e PATH=/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    -e CARGO_TARGET_DIR=/cargo-target -e TMPDIR=/cargo-tmp \
    rust:1.81-bookworm \
    sh -c 'mkdir -p /cargo-target /cargo-tmp && RUSTFLAGS=-Dwarnings cargo check -p omics-alignment --lib'
Finished `dev` profile ... target(s) in 11.30s
```

The controller is Apple Silicon macOS. Docker provided a safe local
`linux/amd64` execution environment, so the following runtime evidence is
containerized x86_64 Linux evidence rather than evidence from a physical Linux
host. The container reported `avx2 detected`, so the feature-gated test
executed its native AVX2 branch:

```text
$ docker run --rm --platform linux/amd64 \
    --mount type=bind,source="$PWD",target=/workspace,readonly \
    --workdir /workspace \
    -e PATH=/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    -e CARGO_TARGET_DIR=/cargo-target -e TMPDIR=/cargo-tmp \
    rust:1.97-bookworm \
    sh -c 'mkdir -p /cargo-target /cargo-tmp && grep -qw avx2 /proc/cpuinfo && echo "avx2 detected" && cargo test -p omics-alignment avx2_matches_scalar_when_available -- --nocapture && cargo test -p omics-alignment avx2_dispatch_matches_scalar -- --nocapture'
avx2 detected
test algorithm::simd::wavefront::tests::avx2_matches_scalar_when_available ... ok
test result: ok. 1 passed; 0 failed; 0 ignored
test algorithm::simd::x86_64_tests::avx2_dispatch_matches_scalar ... ok
test result: ok. 1 passed; 0 failed; 0 ignored
```

The forced-native differential test compares global and local outcomes with
the scalar oracle for both lane widths. Its complete-vector fixture uses
17-byte inputs, which creates a 16-cell interior anti-diagonal for
`Avx2I16`; it includes eight equal and eight unequal substitutions. The
existing nine-byte fixture reaches a complete eight-lane `Avx2I32` vector.
The dispatch test exercises both `i16` and `i32` lane selection against the
scalar engine.

The final native Apple crate suite passed:

```text
$ cargo test -p omics-alignment --all-features
unit tests: 58 passed
integration tests: 6 passed
doc tests: 4 passed
```

I also ran the repository CI commands that apply to this change:

```text
$ cargo +nightly fmt -- --check
exit 0

$ cargo clippy --all-features -- --deny warnings
exit 0

$ cargo test --all-features
exit 0

$ cargo test --all-features --examples
test result: ok. 0 passed; 0 failed

$ cargo doc
exit 0

$ cargo deny --workspace check
advisories ok, bans ok, licenses ok, sources ok

$ (cd omics && cargo msrv verify --output-format minimal --all-features)
true
```

`cargo doc` retains an existing warning in `omics-molecule` for the unresolved
`Nucleotide(s)` link. It does not involve the Task 5 files.

## Safety reasoning

`x86_64.rs` compiles only for `target_arch = "x86_64"` and
`target_os = "linux"`. Both safe dispatch wrappers call
`is_x86_feature_detected!("avx2")` before entering private
`#[target_feature(enable = "avx2")]` functions. Unsupported processors return
directly through `engine::{global, local}`. The native differential test makes
the same runtime check before directly instantiating a private AVX2 `Kernel`.

Every target-feature helper is private and unsafe. Every unsafe block states
the feature or pointer invariant that justifies it. Score loads and stores use
unaligned AVX2 operations but access exactly the `Kernel` contract's sixteen
`i16` or eight `i32` values.

`Avx2I16` reads exactly sixteen reference bytes and the sixteen bytes
immediately before `query_end`, reverses the query block with
`_mm_shuffle_epi8`, compares bytes, sign-widens the all-zero or all-one byte
mask into sixteen `i16` masks, and blends match or mismatch scores.
`Avx2I32` first copies exactly eight reference and reverse-query bytes into
local `i32` arrays, then loads those bounded arrays, compares `i32` lanes, and
blends scores. Signed comparison instructions and byte-wise blending preserve
the generic wavefront's selection and sentinel semantics.

## Files changed

- `omics-alignment/src/algorithm/simd/x86_64.rs` adds the AVX2 kernels,
  target-feature helpers, and checked safe dispatch.
- `omics-alignment/src/algorithm/simd.rs` selects the Linux backend only on
  x86_64 Linux and otherwise retains the existing scalar fallback.
- `omics-alignment/src/algorithm/simd/wavefront.rs` adds Linux-gated
  forced-native differential coverage and adjusts the target-specific
  dead-code expectation.
- `.superpowers/sdd/task-5-report.md` records this work.

I did not edit, stage, or commit `docs/superpowers/plans/`.

## Self-review and concerns

I reviewed target gating, feature detection, lane widths, signed comparison
masks, query reversal, bounded input and score accesses, scalar fallback, and
the full-vector differential fixture. The checked dispatch and AVX2 paths have
Linux x86_64 runtime evidence. The non-AVX2 scalar-fallback branch has
cross-target compile coverage and direct code review but no runtime evidence
because the available Linux x86_64 container reports AVX2.

`rust:1.81-bookworm` cannot run this crate's test target because its Cargo
parser rejects the locked `clap_lex` manifest's `edition2024` feature before
building `omics-alignment`. The warning-free Rust `1.81` library compilation
above verifies the production AVX2 source at the workspace MSRV; Rust `1.97`
executes the Linux x86_64 test target.

An additional non-CI diagnostic,
`cargo clippy -p omics-alignment --all-features --tests --target
x86_64-unknown-linux-gnu -- --deny warnings`, flags the unchanged
`TestI16::splat` loop as `clippy::manual_slice_fill` under the current Clippy.
The exact repository CI lint command does not lint test targets and passes.
