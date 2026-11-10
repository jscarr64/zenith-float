# Calibration (checklist item S0)

Plain facts, measured 2026-10-15 (evening, America/Chicago) on Pop (Intel Core i9-14900KF, 32 logical cores, cargo 1.99.1). Base commit 078533b plus S1 (3e89122).
Other agents were using the machine (load average 4 to 10 during the runs), so treat the numbers as upper-ish, not exact. Each timing was run 3 times; the range is given.
Wall seconds are the `time` of the already-built test binary (compile time excluded), unless stated.

## What could NOT be timed

`src/mag` does not exist at this head, so `cargo test --test mag` does not compile (8 errors, all `unresolved import zenith_bound::mag`). No Mag op, no property test and no `mag_ops_match_mpmath_goldens` could be timed. Nothing below is a time of a real Mag implementation.

## Model / golden families (model only; scratch stand-in for Mag)

The model/golden checks live in the `mag` target, which does not compile. To time them anyway, a scratch copy was built in /tmp/zb-harper-scratch/calib (outside this repo and outside /home/jsc64/MathWorks): this crate's tests/mag, tests/model and tests/goldens plus a stub `src/mag.rs` whose ops return Infinity (so `Mag`-dependent tests fail there; they were not timed). The filters below run only tests that do not call a Mag op.

| Family (filter) | Tests | Profile | Wall seconds (3 runs) |
|---|---|---|---|
| `model_selftest::` (tests/mag/model_selftest.rs) | 4 | debug | 0.001 to 0.002 |
| `model_selftest::` | 4 | release | 0.001 to 0.002 |
| `golden::model_` (model accepts every golden row; model rejects every perturbed row) over 20,954 rows of tests/goldens/mag_ops.txt (non-comment lines; count by `grep -vc`) | 2 | debug | 0.080 to 0.081 |
| `golden::model_` | 2 | release | 0.010 to 0.011 |
| `golden::golden_file_has` (all eight ops present) | 1 | debug | 0.008 to 0.011 |
| `golden::golden_file_has` | 1 | release | 0.003 |

Not timed: `golden::mag_ops_match_mpmath_goldens` (needs real Mag) and every other `mag` test.

## splitmix64 draw loop (throwaway test in the scratch crate, `std::hint::black_box` on the seed, count and result)

| Draws | Profile | Wall seconds (3 runs, in-process `Instant`) |
|---|---|---|
| 1,000,000 | debug | 0.0043 to 0.0061 |
| 1,000,000 | release | 0.0012 to 0.0014 |
| 100,000,000 | debug | 0.387 to 0.390 |
| 100,000,000 | release | 0.076 to 0.082 |

Arithmetic from the 100,000,000 rows (not separately measured): about 3.9 ns per draw in debug, about 0.8 ns per draw in release.

## tests/harness (S1 self-tests; was tests/common_selftest)

10 tests: debug 0.031 to 0.032 s wall (test binary), release 0.003 s wall. Not broken down per test.

## Reading for D7 (final case counts)

Only these facts are established: the generator is not the bottleneck (1,000,000 draws are about 5 ms in debug); the independent model checks 20,954 golden rows twice in about 0.08 s debug. Case counts for any property suite depend on the cost of the real Mag op plus the model check per case, which cannot be measured until src/mag exists. Re-run S0 (same commands, real Mag) when it does, before fixing the section 1 counts.
