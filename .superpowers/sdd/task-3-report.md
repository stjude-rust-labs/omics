# Task 3 Report: Generic Wavefront and Traceback

Implemented the private, target-independent affine-gap wavefront engine with a
generic `Kernel`, narrow `i16` and `i32` score lanes, packed traceback fields,
and an eight-lane scalar-array test backend. The public byte API remains on its
documented scalar fallback because this task intentionally adds no NEON or AVX2
backend.

## RED evidence

I added `TestI16`, exhaustive binary differential helpers, and the required
global and local tests before adding the `Kernel` or wavefront implementation.

```text
$ cargo test -p omics-alignment test_kernel_matches_scalar -- --nocapture
error[E0405]: cannot find trait `Kernel` in this scope
error[E0425]: cannot find function `global` in this scope
error[E0425]: cannot find function `local` in this scope
error: could not compile `omics-alignment` (lib test) due to 5 previous errors
```

The initial GREEN implementation made both required exhaustive tests pass.
During code review, I found that generic kernels could narrow an otherwise
valid scalar alignment outside the selected lane range. I added
`test_kernel_falls_back_outside_the_narrow_range` before correcting it. Its RED
run failed as expected:

```text
$ cargo test -p omics-alignment test_kernel_falls_back_outside_the_narrow_range -- --nocapture
Error: ScoreOverflow
test ...::test_kernel_falls_back_outside_the_narrow_range ... FAILED
```

The fix shares the checked full-path bound with `lane_width` and falls back to
the scalar engine before narrow conversion when the selected lane cannot prove
the range. The regression test then passed.

## GREEN evidence

```text
$ cargo test -p omics-alignment test_kernel_matches_scalar -- --nocapture
running 3 tests
test ...::test_kernel_matches_scalar_for_complete_vectors ... ok
test ...::test_kernel_matches_scalar_for_short_local_inputs ... ok
test ...::test_kernel_matches_scalar_for_short_global_inputs ... ok
test result: ok. 3 passed; 0 failed

$ cargo test -p omics-alignment test_kernel_falls_back_outside_the_narrow_range -- --nocapture
running 1 test
test ...::test_kernel_falls_back_outside_the_narrow_range ... ok
test result: ok. 1 passed; 0 failed
```

The required exhaustive tests enumerate every binary byte sequence through
length five under all four existing scoring fixtures and compare complete
successful outcomes plus the empty-global error against `engine`. An additional
focused test exercises a complete eight-lane vector diagonal, including local
and global traceback. The narrow-range regression covers both an individual
score too wide for `i16` and an in-range score whose full-path bound exceeds
`i16`.

## Verification

`cargo test -p omics-alignment --all-features` passed with 53 unit tests, 6
integration tests, and 4 documentation tests. The focused crate lint passed:

```text
$ cargo clippy -p omics-alignment --all-features -- --deny warnings
Finished `dev` profile ... successfully
```

I also ran the primary CI commands from `.github/workflows/CI.yml` successfully:

```text
cargo fmt -- --check
cargo clippy --all-features -- --deny warnings
cargo test --all-features
cargo test --all-features --examples
cargo doc
```

`cargo fmt -- --check` emitted the repository's existing stable-toolchain
warnings for nightly-only formatting settings but exited successfully. I ran
`cargo +nightly fmt` before the final test and lint passes. The remaining CI
checks also passed:

```text
cargo deny --workspace check
cargo msrv verify --output-format minimal --all-features
```

The MSRV command ran from `omics/`, matching its workflow working directory.

## Commit

`3516d7a` — feat: add generic `wavefront` alignment kernel

## Files changed

- `omics-alignment/src/algorithm/simd/wavefront.rs`
  - Adds `Narrow`, the private `Kernel` contract, generic global and local
    wavefront execution, rolling score diagonals, scalar tails, endpoint
    selection, packed traceback, and differential test kernel coverage.
- `omics-alignment/src/algorithm/engine.rs`
  - Exposes existing validation, dimension, and CIGAR helpers to sibling
    private modules so both engines preserve their construction behavior.
- `omics-alignment/src/algorithm/simd.rs`
  - Documents the intentional scalar fallback pending target-specific kernels.

## Self-review

The recurrence keeps `Narrow::MIN` unreachable and applies each candidate with
strict replacement in the scalar ordering `A,D,I`, `I,A,D`, or `D,A,I`.
Vector candidates retain unreachable predecessors before their native add,
whereas scalar tails use checked scalar addition. Local aligned scores reset
only when non-positive, with no predecessor field. Endpoint selection uses
`(score, Reverse(row), Reverse(column))`, rather than anti-diagonal traversal
order. Traceback reads the relevant two-bit field, recomputes `=` versus `X`,
and reuses `engine::path_to_cigar`.

A focused code review identified the missing lane-eligibility guard. The final
implementation uses the same conservative bound as `lane_width` and falls back
to `engine::{global, local}` whenever a kernel's narrow maximum cannot prove
the bound. This preserves scalar scores and `ScoreOverflow` behavior outside a
kernel's valid lane range.

## Concerns

No target-specific SIMD backend exists yet, so public `algorithm::simd`
functions deliberately retain scalar performance. The generic wavefront stores
one traceback byte per matrix cell and allocates nine rolling score arrays plus
small predecessor staging buffers. This task establishes exactness and backend
contracts; it does not include NEON or AVX2 benchmarking.

## Fix/test

`omics-alignment/src/algorithm/simd/wavefront.rs` now replaces the three
ad hoc complete-vector fixtures with deterministic `TestI16` differential
cases whose lengths prove which path the eight-lane loop must take.
`test_vector_differential_cases_cover_test_i16_paths` first shows why the
existing exhaustive binary sweep through length five never reaches the vector
loop—`min(5, 5) / 8 == 0`—then proves that the new `(8, 8)`, `(9, 9)`, and
`(16, 16)` inputs force one full vector, one full vector plus one scalar-tail
cell, and two full vector chunks. `test_kernel_matches_scalar_for_global_vector_paths`
and `test_kernel_matches_scalar_for_local_vector_paths` run those cases through
every `exhaustive_scorings()` fixture and compare `global::<TestI16>` and
`local::<TestI16>` against `engine::global` and `engine::local`.

I verified the harness with:

```text
cargo test -p omics-alignment test_vector_differential_cases_cover_test_i16_paths -- --nocapture
cargo test -p omics-alignment test_kernel_matches_scalar -- --nocapture
cargo test -p omics-alignment --all-features
```
