# Task 7 Report: Benchmarks, CI, and Documentation

Task 7 adds end-to-end scalar versus SIMD alignment benchmarks, native Apple
Silicon CI coverage, public dispatch documentation, and a focused hot-path
optimization required to satisfy the Linux AVX2 performance gate.

## Benchmark design

`omics-alignment/benches/algorithms.rs` now measures scalar
`algorithm::{global, local}` and byte-dispatched `algorithm::simd::{global,
local}` separately. Every iteration passes both input slices and `Scoring`
through `black_box`, then consumes the complete returned `Outcome`.

The suite registers `16`, `64`, `256`, and `1024` base reference lengths for
each algorithm, implementation, and fixture. The fixtures are `balanced`
equal-length inputs with substitutions, `asymmetric` inputs with a half-length
query, and `gap-heavy` inputs that insert a query-only four-byte region in
each 24-byte reference block. Throughput uses the actual matrix dimensions for
each fixture.

## Performance investigation and optimization

The first Linux `amd64` Docker measurement did not meet the gate. Its balanced
Criterion medians were as follows.

| Algorithm | Length | Scalar median | SIMD median | SIMD/scalar |
| --- | ---: | ---: | ---: | ---: |
| global | 256 | 1.3033 ms | 1.4901 ms | 1.1433 |
| global | 1024 | 22.5744 ms | 22.8204 ms | 1.0109 |
| local | 256 | 1.6305 ms | 1.6998 ms | 1.0425 |
| local | 1024 | 27.0529 ms | 26.7097 ms | 0.9873 |

The Docker CPU reported `vendor_id: VirtualApple` and advertised `avx2`.
Disassembly showed that the generic vector recurrence made repeated calls to
the target-feature AVX2 wrappers inside its innermost loop, including 15
static `splat_avx2` call sites. Repeated broadcasts of invariant scores and
state codes caused the fixed overhead that dominated the translated AVX2
environment.

I added a focused RED test named
`vector_values_broadcast_recurrence_constants`. It initially failed with
`cannot find type VectorValues in this scope`. The green implementation adds
`VectorValues<K>`, which broadcasts unreachable, zero, gap, and predecessor
constants once before the diagonal loop. `choose_vector` and `fill_vector`
reuse those vectors without changing recurrence order or tie behavior. The
focused test passed, then the full Apple and Linux SIMD suites passed.

The following final tables use Criterion
`new/estimates.json` `median.point_estimate`, converted to milliseconds.
`SIMD/scalar` below one is faster. The balanced rows are the acceptance
fixture, and the additional rows show that the new asymmetric and gap-heavy
coverage also improves on both targets.

### Apple Silicon NEON

| Fixture | Algorithm | Length | Scalar median | SIMD median | SIMD/scalar | Speedup |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| balanced | global | 256 | 0.8147 ms | 0.1041 ms | 0.1278 | 7.826x |
| balanced | global | 1024 | 13.2001 ms | 1.3567 ms | 0.1028 | 9.730x |
| balanced | local | 256 | 0.8849 ms | 0.1469 ms | 0.1661 | 6.022x |
| balanced | local | 1024 | 14.3010 ms | 1.9836 ms | 0.1387 | 7.210x |
| asymmetric | global | 256 | 0.4016 ms | 0.0532 ms | 0.1326 | 7.543x |
| asymmetric | global | 1024 | 6.4793 ms | 0.6369 ms | 0.0983 | 10.174x |
| asymmetric | local | 256 | 0.4294 ms | 0.0756 ms | 0.1761 | 5.677x |
| asymmetric | local | 1024 | 6.9617 ms | 0.9634 ms | 0.1384 | 7.226x |
| gap-heavy | global | 256 | 0.9395 ms | 0.1147 ms | 0.1221 | 8.188x |
| gap-heavy | global | 1024 | 15.1684 ms | 1.4641 ms | 0.0965 | 10.360x |
| gap-heavy | local | 256 | 1.0782 ms | 0.1644 ms | 0.1524 | 6.560x |
| gap-heavy | local | 1024 | 16.3925 ms | 2.2849 ms | 0.1394 | 7.174x |

### Linux `amd64` AVX2 Docker

| Fixture | Algorithm | Length | Scalar median | SIMD median | SIMD/scalar | Speedup |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| balanced | global | 256 | 1.3596 ms | 1.2343 ms | 0.9078 | 1.102x |
| balanced | global | 1024 | 22.0603 ms | 18.6228 ms | 0.8442 | 1.185x |
| balanced | local | 256 | 1.6466 ms | 1.3456 ms | 0.8172 | 1.224x |
| balanced | local | 1024 | 26.8574 ms | 20.5140 ms | 0.7638 | 1.309x |
| asymmetric | global | 256 | 0.6814 ms | 0.6188 ms | 0.9081 | 1.101x |
| asymmetric | global | 1024 | 10.6937 ms | 9.2520 ms | 0.8652 | 1.156x |
| asymmetric | local | 256 | 0.7936 ms | 0.6780 ms | 0.8544 | 1.170x |
| asymmetric | local | 1024 | 13.5564 ms | 10.3347 ms | 0.7623 | 1.312x |
| gap-heavy | global | 256 | 1.5600 ms | 1.4250 ms | 0.9135 | 1.095x |
| gap-heavy | global | 1024 | 25.5891 ms | 21.6860 ms | 0.8475 | 1.180x |
| gap-heavy | local | 256 | 1.8719 ms | 1.5894 ms | 0.8491 | 1.178x |
| gap-heavy | local | 1024 | 31.4778 ms | 24.2766 ms | 0.7712 | 1.297x |

The acceptance gate passes. On both targets, balanced global and local SIMD
medians improve at `1024`; neither result regresses at `256`. All additional
fixture rows also improve.

## Commands and environments

Apple Silicon host:

```bash
cargo bench -p omics-alignment --bench algorithms
cargo test -p omics-alignment --all-features
```

The host reports `Darwin ... arm64`; native NEON tests ran. The focused suite
passed 60 unit tests, 9 integration tests, and 5 documentation tests.

Linux AVX2 used a read-only source mount and a writable ignored target
directory. The command requires `/proc/cpuinfo` to advertise `avx2` before it
runs tests or benchmarks.

```bash
docker run --rm --platform linux/amd64 \
  --mount type=bind,source="$PWD",target=/workspace,readonly \
  --mount type=bind,source="$PWD/target/task-7-linux",target=/cargo-target \
  --workdir /workspace \
  -e PATH=/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e CARGO_TARGET_DIR=/cargo-target \
  -e TMPDIR=/cargo-tmp \
  rust:1.97-bookworm \
  sh -c 'mkdir -p /cargo-tmp && grep -qw avx2 /proc/cpuinfo && echo "avx2 detected" && cargo bench -p omics-alignment --bench algorithms'

docker run --rm --platform linux/amd64 \
  --mount type=bind,source="$PWD",target=/workspace,readonly \
  --mount type=bind,source="$PWD/target/task-7-linux",target=/cargo-target \
  --workdir /workspace \
  -e PATH=/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e CARGO_TARGET_DIR=/cargo-target \
  -e TMPDIR=/cargo-tmp \
  rust:1.97-bookworm \
  sh -c 'mkdir -p /cargo-tmp && grep -qw avx2 /proc/cpuinfo && echo "avx2 detected" && cargo test -p omics-alignment --all-features'
```

Both Docker commands printed `avx2 detected`. The full Linux focused suite
passed 60 unit tests, including AVX2 dispatch and wavefront tests, 9
integration tests including both 1,000-case public differential tests, and 5
documentation tests.

## CI and documentation

GitHub's authoritative
[GitHub-hosted runners reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners)
lists `macos-14` as an arm64 M1 GitHub-hosted runner. The `test` job now
matrices `ubuntu-latest` and `macos-14`; it preserves the existing stable
toolchain setup and literal `cargo test --all-features` command. The Linux
format, lint, example, documentation, deny, and MSRV jobs are unchanged.

The `algorithm` module now documents byte-slice dispatch with a compile-tested
scalar-parity example. It specifies NEON on macOS `aarch64`, runtime AVX2 on
Linux `x86_64`, scalar behavior on unsupported targets and processors, exact
result parity, `i16` then `i32` lane selection, and scalar fallback. Crate
documentation summarizes the same public behavior. The Unreleased changelog
records NEON and AVX2 global and local alignment support.

## Verification

The required repository commands completed successfully after final
formatting.

```bash
cargo fmt -- --check
cargo clippy --all-features -- --deny warnings
cargo test --all-features
cargo test --all-features --examples
cargo doc
```

`cargo clippy` completed with `--deny warnings`. The workspace test command
completed all listed unit, integration, and documentation suites successfully.
The examples command completed successfully. I also ran the unchanged CI
checks:

```bash
cargo +nightly fmt -- --check
cargo deny --workspace check
(cd omics && cargo msrv verify --output-format minimal --all-features)
git diff --check
```

Nightly formatting passed, `cargo deny` reported advisories, bans, licenses,
and sources as `ok`, the MSRV command printed `true`, and `git diff --check`
reported no whitespace errors. A dedicated read-only review found no
high-confidence issues.

## Files changed

* `.github/workflows/CI.yml`
* `omics-alignment/CHANGELOG.md`
* `omics-alignment/benches/algorithms.rs`
* `omics-alignment/src/algorithm.rs`
* `omics-alignment/src/algorithm/simd.rs`
* `omics-alignment/src/algorithm/simd/wavefront.rs`
* `omics-alignment/src/lib.rs`
* `omics-alignment/tests/algorithm.rs`
* `.superpowers/sdd/task-7-report.md`

`simd.rs` and `tests/algorithm.rs` contain formatter-only changes required for
the existing nightly format job. I did not edit, stage, or commit
`docs/superpowers/plans/`; its pre-existing untracked directory remains in
`git status --short`.

## Self-review and concerns

I checked the final diff for whitespace and scope, verified that the Criterion
IDs separate scalar and SIMD paths, confirmed the fixture inputs remain
positive for local alignment, checked that the CI matrix leaves existing
commands intact, and reviewed the constant-cache change against the original
recurrence order. The review found no correctness, CI, documentation, or MSRV
issue.

`cargo fmt -- --check` on the local stable toolchain emits existing warnings
for unstable `rustfmt.toml` settings, although it exits successfully. The
nightly CI-equivalent formatter check passes. `cargo doc` exits successfully
but retains one pre-existing `omics-molecule` warning for the unresolved
`Nucleotide(s)` link; Task 7 adds no rustdoc warnings. Finally, Docker reports
`VirtualApple` for the Linux CPU, so its AVX2 measurements validate the
required safe `linux/amd64` environment and runtime dispatch but are not a
substitute for measurements on physical AMD64 hardware.
