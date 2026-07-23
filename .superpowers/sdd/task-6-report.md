# Task 6 Report: Edge Cases and Public Differential Coverage

Task 6 hardens `omics-alignment` exactness with deterministic public regression cases and fixed-seed randomized differential coverage for the native SIMD backends. The work adds forced-backend regression helpers in `omics-alignment/src/algorithm/simd.rs` and public differential coverage in `omics-alignment/tests/algorithm.rs`.

## RED evidence

I added the new public regression and randomized tests before the supporting helpers. The first focused run failed for the expected reason: the new tests referenced missing coverage helpers.

```text
$ cargo test -p omics-alignment simd_public_api_matches_scalar_for_edge_case_regressions --all-features
error[E0425]: cannot find function `public_simd_regression_cases` in this scope
error[E0425]: cannot find function `assert_public_simd_case` in this scope
error[E0425]: cannot find function `backend_regression_cases` in this scope
error[E0425]: cannot find function `assert_random_public_simd_global_cases` in this scope
error[E0425]: cannot find function `assert_random_public_simd_local_cases` in this scope
```

A later focused run also caught an insufficient `NeonI32` fixture. The test failed because the first mixed-substitution attempt did not actually place equal and unequal anti-diagonal substitutions in the intended shape.

```text
$ cargo test -p omics-alignment neon_i32_mixed_substitutions_reach_two_full_vectors --all-features
FAILED
thread 'algorithm::simd::tests::neon_i32_mixed_substitutions_reach_two_full_vectors' panicked at ...
neon-i32-mixed-substitutions should include equal anti-diagonal substitutions
```

Those failures drove the final helpers, per-vector mixed-lane assertions, the corrected equal-local-maxima regression, and explicit native-backend availability assertions for the 1,000-case randomized runs.

## GREEN evidence

The focused Apple Silicon runs passed after the fixes.

```text
$ cargo test -p omics-alignment neon_i32_mixed_substitutions_reach_two_full_vectors --all-features
test algorithm::simd::tests::neon_i32_mixed_substitutions_reach_two_full_vectors ... ok

$ cargo test -p omics-alignment simd_public_api_matches_scalar_for_edge_case_regressions --all-features
test simd_public_api_matches_scalar_for_edge_case_regressions ... ok

$ cargo test -p omics-alignment simd_public_api_matches_scalar_for_1000_random_global_cases --all-features
test simd_public_api_matches_scalar_for_1000_random_global_cases ... ok

$ cargo test -p omics-alignment simd_public_api_matches_scalar_for_1000_random_local_cases --all-features
test simd_public_api_matches_scalar_for_1000_random_local_cases ... ok
```

The required crate suite also passed on the Apple Silicon host.

```text
$ cargo test -p omics-alignment --all-features
unit tests: 59 passed
integration tests: 9 passed
doc tests: 4 passed
```

The Linux `amd64` Docker run exercised the AVX2 path through the existing safe container approach. The container reported `avx2 detected`, then passed the AVX2 dispatch, the AVX2 mixed-substitution `i32` case, the public regression test, and the 1,000 global plus 1,000 local randomized public differential tests.

```text
$ docker run --rm --platform linux/amd64 \
    --mount type=bind,source="$PWD",target=/workspace,readonly \
    --workdir /workspace \
    -e PATH=/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    -e CARGO_TARGET_DIR=/cargo-target -e TMPDIR=/cargo-tmp \
    rust:1.97-bookworm \
    sh -c 'mkdir -p /cargo-target /cargo-tmp && grep -qw avx2 /proc/cpuinfo && echo "avx2 detected" && cargo test -p omics-alignment avx2_dispatch_matches_scalar --all-features && cargo test -p omics-alignment avx2_i32_mixed_substitutions_reach_two_full_vectors --all-features && cargo test -p omics-alignment simd_public_api_matches_scalar_for_edge_case_regressions --all-features && cargo test -p omics-alignment simd_public_api_matches_scalar_for_1000_random_global_cases --all-features && cargo test -p omics-alignment simd_public_api_matches_scalar_for_1000_random_local_cases --all-features'
avx2 detected
test algorithm::simd::x86_64_tests::avx2_dispatch_matches_scalar ... ok
test algorithm::simd::x86_64_tests::avx2_i32_mixed_substitutions_reach_two_full_vectors ... ok
test simd_public_api_matches_scalar_for_edge_case_regressions ... ok
test simd_public_api_matches_scalar_for_1000_random_global_cases ... ok
test simd_public_api_matches_scalar_for_1000_random_local_cases ... ok
```

I also reran the relevant repository CI commands from `.github/workflows/CI.yml`. They all passed after the Task 6 changes.

```text
$ cargo fmt -- --check
exit 0

$ cargo clippy --all-features -- --deny warnings
exit 0

$ cargo test --all-features
exit 0

$ cargo test --all-features --examples
exit 0

$ cargo doc
exit 0

$ cargo deny --workspace check
advisories ok, bans ok, licenses ok, sources ok

$ (cd omics && cargo msrv verify --output-format minimal --all-features)
true
```

## Coverage summary

The public regression table now covers empty global and local inputs, equal local maxima on different diagonals, zero-cost gaps, lengths around the 4-, 8-, and 16-lane native widths, asymmetric `1 x 257` inputs, `Score::MIN`, and scores that force `i32` or scalar fallback. Each public case compares complete `Result` values, including error parity through string and debug-form comparisons when needed.

The randomized public coverage uses a private fixed-seed `Xorshift64` helper, no extra dependency, sequence lengths through `128`, and exactly 1,000 global plus 1,000 local native-backend cases per supported native environment. The helpers now assert native-backend availability instead of silently passing without executing the required case count.

The forced-backend coverage adds backend-specific mixed-substitution multi-vector fixtures for `NeonI32` and `Avx2I32`. Each fixture proves the `i32` lane shape and now asserts mixed equal and unequal anti-diagonal substitutions inside every full vector chunk, not just somewhere in the whole input.

## Files changed

`omics-alignment/src/algorithm/simd.rs` adds backend regression fixtures, full-result comparison helpers, and backend-specific mixed-substitution vector-shape assertions. `omics-alignment/tests/algorithm.rs` adds the public regression table, the private fixed-seed `Xorshift64` helper, native-backend randomized differential tests, and the explicit equal-local-maxima diagonal check. `.superpowers/sdd/task-6-report.md` records this work.

I did not edit, stage, or commit `docs/superpowers/plans/`.

## Self-review

I requested a read-only review before finishing. The reviewer correctly flagged that my first equal-local-maxima case landed on one diagonal, that my first mixed-substitution fixtures only mixed equality across whole vectors rather than within each vector chunk, and that the randomized helpers should not silently succeed when no native backend runs. I fixed all three points, reran the focused Apple and Linux checks, and reran `cargo test -p omics-alignment --all-features`.

I also verified that the randomized tests stay deterministic, that the public checks compare complete `Result` values, and that the new backend fixtures reach the intended `NeonI32` and `Avx2I32` lane counts.

## Concerns

`cargo fmt -- --check` succeeds, but stable `rustfmt` still prints the repository's existing unstable-configuration warnings. `cargo doc` also retains the existing unresolved `Nucleotide(s)` link warning in `omics-molecule`. Neither warning comes from the Task 6 files.
