# `expr!` rounding contract

`expr!` evaluates a tree of `ExactNum` operations at extra working precision, then rounds **once** to the context precision `p` with the context rounding mode `rm`.

This is the Accumath numeric contract: leaves are treated as exact; the printed result is `set_precision(p, rm)` of the last value, after exponent-window clamping.

## Pipeline

1. **Context** supplies `p`, `rm`, `Consts`, and exponent limits `emin` / `emax`.
2. **Working precision** starts at `p_wrk = p + WORD_BIT_SIZE`, plus per-subexpression cancellation slots `errs[]` (add/sub/mul detect loss of significance and `continue` the retry loop with a larger `p_wrk`).
3. **Internal ops** use `RoundingMode::None` so intermediate values keep sticky bits.
4. **Final round:** `ret.set_precision(p, rm)`, then `check_exponent_range` (outside `[emin, emax]` becomes `0` or `Inf`).

`MAX_PREC_RETRY` still bounds kernel retries inside each leaf (sin, ln, …). The macro loop itself grows `errs[]` rather than that constant.

## Per-operation semantics

| Syntax | Kernel | Rounding |
| -------- | -------- | ---------- |
| `+` `-` `*` `/` `%` | `add` / `sub` / `mul` / `div` / `rem` | Working precision, `None`; add/sub bump `errs` on cancellation |
| `recip(x)` | `reciprocal` | Same |
| `sqrt` / `cbrt` / `root(x, n)` | `sqrt` / `cbrt` / `nth_root` | Kernel Ziv/`try_set_precision` to `p_wrk`, then final `set_precision` |
| `ln` `log2` `log10` `log` `log1p` | corresponding `ExactNum` methods | Need `Consts`; extra bits for argument reduction |
| `exp` `exp2` `exp10` `expm1` `pow` | `exp` / `pow` | Same |
| `rem_pi(x)` | `rem_pi` | Reduction into `(-2π, 2π)`, not MPFR `fmod` |
| `sin` `cos` `tan` `asin` `acos` `atan` `atan2` | trig | `rem_pi` / series; cost grows with `\|e\|` |
| `hypot` `fma` / `mul_add` | `hypot` / `fma` | `fma` is a single round of `a*b+c` (full product when magnitudes overlap) |
| `sinh` `cosh` `tanh` `asinh` `acosh` `atanh` | hyperbolic | Paired `sinh_cosh` is not a macro leaf; call the method |
| `erf` `erfc` `gamma` `ln_gamma` `bessel_j` | specials | Same pipeline; MPFR oracles are unit/compare coverage, not the full 1000-iter loop for Γ |

Integer literals and `exact!` / `fbig!` strings enter as exact `ExactNum` values.

## What `expr!` does **not** guarantee

- Bit-identity with MPFR for every leaf (trig/log use the same kernel as `ExactNum`, which is 1-ULP oracles for the compare suite, not a proof for composite trees).
- Correct rounding of the **entire** expression in the IEEE fused sense. Only the **root** of the tree is rounded to `p`; internal `None` rounding plus `errs[]` is the cancellation heuristic in `doc/README.md`.
- Unbounded exponent: set `emin`/`emax` to the real range or precision (and time) grow with unused magnitude.

## Practical rule

For Accumath terminals: put the formula in `expr!`, use a `Context` whose `p` is the output precision, `rm` is the required mode (usually `ToEven`), and `emin`/`emax` bound the formula class.
