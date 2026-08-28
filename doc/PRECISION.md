# Working precision and exponent growth

zenith-float is arbitrary-precision software floating-point: cost scales with the **bit width of the mantissa** being carried, not with hardware IEEE binary64. Several mechanisms bound how far internal precision can grow relative to the user-requested precision `p` and the operand exponent `e`.

## Correct-rounding retries

Transcendental and decimal conversion paths use a working precision `p_wrk` that starts near `p` (often `p + WORD_BIT_SIZE`) and may increase when a correct-rounding probe fails.

- **`MAX_PREC_RETRY`** (`zenith_float_num::MAX_PREC_RETRY`, default `256`): maximum number of word-sized retry steps beyond `p`.
- **`prec_retry_exhausted(p_wrk, p)`** returns true when `p_wrk > p + WORD_BIT_SIZE * MAX_PREC_RETRY`.
- **`bump_prec_retry`** increases `p_wrk` or returns `Error::InvalidArgument` instead of looping without bound.

At 64-bit limbs this caps extra work at roughly `256 * 64 = 16 384` bits above `p` per operation.

Public helpers:

- **`ziv_round(p, rm, compute)`** — same loop as the kernel: call `compute(p_wrk)`, then `try_set_precision(p, rm, p_wrk)` until the rounding is unique.
- **`Ball { mid, rad }`** — first-order interval arithmetic (`add` / `mul`) plus a rounding ulp, for proving a midpoint is the unique correctly rounded value.

## Stack inlining

Mantissa buffers of at most **`INLINE_WORDS` (2)** limbs (`128` bits at 64-bit `Word`) live on the stack. Larger precisions promote to the heap.

## Constant cache

`Consts` (`ConstCache`) keeps series state for π, e, ln 2, and ln 10 and extends it when more bits are requested. √2 and φ use a `CachedFBig`-style extra cache (`cache_info()` reports bit lengths).

## Benchmark history

Tab-separated logs (no JSON):

| File | Source |
| ------ | -------- |
| `doc/bench-baselines.tsv` | Criterion ids from `scripts/bench.sh` / `scripts/bench-compare.sh` |
| `doc/compare-results.tsv` | `scripts/compare-bench.sh` vs astro-float (132-bit core ops) |

Re-run `./scripts/compare-bench.sh --quick` before a crates.io publish to refresh the compare table. Dashu is optional (`--dashu`).

## `expr!` macro

The `expr!` macro (see `zenith-float-macro`) uses:

1. **Target precision** `p` and rounding mode from `Context`.
2. **Initial working precision** `p_rnd = p + WORD_BIT_SIZE`.
3. **Per-subexpression error slots** `errs[]` bumped when cancellation is detected; each retry uses `p_wrk = p_rnd + sum(errs)`.
4. **Exponent window** `emin` / `emax` from `Context`; results outside the window become `Inf` / `0` via `macro_util::check_exponent_range`.

Tighten `emin` / `emax` to the real exponent range of your computation so the macro does not reserve precision for unused magnitude.

## Operand exponent in transcendentals

Many elementary functions add terms proportional to `|e|` when the argument is far from the reduction interval:

| Pattern | Typical extra bits |
| -------- | ------------------- |
| Trig / hyperbolic arg reduction (`sin`, `cos`, `tan`, …) | `O(|e|)` for `rem_pi` / series region |
| `rem_pi` π cache | `min(|e|, 65536)` added to π precision |
| `asinh` large \|x\| | `2|e|` extra working bits |
| `log` / `pow` decimal scaling | bounded by `TEN_PWR_MAX_*` chunking in `conv.rs` |

Property tests in `ops/tests.rs` cap random exponents at **`TEST_EXP_BOUND = 1024`** to keep CI fast; MPFR oracles use wider ranges on a release gate.

## Parse / format

- Decimal output length is estimated from mantissa bit length (`~p * log10(2)`).
- General radix (2–36) conversion uses at most `O(p)` digits per conversion pass.
- Retry loops in `conv_to_dec` share the same **`MAX_PREC_RETRY`** cap as transcendentals.

## Practical guidance

1. Prefer **`RoundingMode::None`** for intermediate steps; round once at the end with `set_precision` / `round`.
2. Set **`Context` exponent limits** to the smallest interval that contains real results.
3. For huge arguments to trig functions, expect cost to grow with **`rem_pi`** precision (capped at `65536` extra bits for π).
4. If `Error::InvalidArgument` appears after many retries, increase `p` or reduce exponent spread; the library has hit the **`MAX_PREC_RETRY`** guard.

See also `doc/README.md` for ulp error bounds used by `expr!` cancellation heuristics.
