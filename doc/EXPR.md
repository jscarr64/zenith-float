# `expr!` rounding contract

`expr!` evaluates a tree of `ExactNum` operations at extra working precision, then rounds **once** to the context precision `p` with the context rounding mode `rm`.

Leaves are treated as exact; the printed result is `set_precision(p, rm)` of the last value, after exponent-window clamping.

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
| `rem_pi(x)` | `rem_pi` | Reduction into `(-2π, 2π)`; identity + range; trig oracles on unreduced `x` |
| `sin` `cos` `tan` `asin` `acos` `atan` `atan2` | trig | `rem_pi` / series; cost grows with `\|e\|` |
| `hypot` `fma` / `mul_add` | `hypot` / `fma` | `fma` is a single round of `a*b+c` (full product when magnitudes overlap) |
| `sinh` `cosh` `tanh` `asinh` `acosh` `atanh` | hyperbolic | Paired `sinh_cosh` is the `ExactNum` method; expr uses `sinh`/`cosh` |
| `erf` `erfc` `gamma` `ln_gamma` `digamma` `gammainc` `gammainc_upper` `ei` `si` `ci` `li` `fresnel_s` `fresnel_c` `bessel_j` `bessel_j_nu` `bessel_y` `bessel_i` `bessel_k` `elliptic_k` `elliptic_e` `elliptic_e_inc` `elliptic_f` `elliptic_pi` `elliptic_pi_inc` `jacobi_am` `jacobi_sn` `jacobi_cn` `jacobi_dn` `jacobi_cd` `jacobi_ns` `jacobi_nc` `jacobi_nd` `jacobi_sc` `jacobi_sd` `jacobi_cs` `jacobi_ds` `jacobi_dc` `legendre_p` `legendre_p_assoc` `hypergeom_2f1` `betainc` `normal_pdf` `normal_cdf` `gamma_pdf` `beta_pdf` `poisson_pmf` `binomial_pmf` `chi_squared_cdf` `student_t_pdf` | specials | Same pipeline; MPFR 1-ULP oracles on bounded domains where MPFR has the function |
| `ldexp(x, n)` `scalb(x, n)` `logb(x)` | IEEE split | Integer `n`; `frexp`/`ilogb` are methods (tuple / `Option`) |
| `pi` `e` `ln_2` `ln_10` `sqrt2` `phi` `euler_gamma` | `Consts` | Cached at extra bits, then final `set_precision` |

Integer literals and `exact!` / `fbig!` strings enter as exact `ExactNum` values.

## `cexpr!`

Same pipeline on `ExactComplex`: working precision, **two** `errs[]` slots per add/sub (real and imaginary independently — `complex_cancel_bits` returns a pair), one final `set_precision` on re and im, then `check_complex_exponent_range`.

Complex `sin`/`cos`/`tan` use the `x+iy` identities. They do **not** call `rem_pi` on the complex value.

Full leaf list, branch cuts, and explicit exclusions: [LIBRARY.md](LIBRARY.md) §17.

Do not use `e` as a variable name (`e` is Euler’s number). Imaginary unit is `I`. No `atan2` or `rem_pi` in `cexpr!`.

## What `expr!` does **not** guarantee

- Bit-identity with MPFR for every leaf (trig/log use the same kernel as `ExactNum`, which is 1-ULP oracles for the compare suite, not a proof for composite trees).
- Correct rounding of the **entire** expression in the IEEE fused sense. Only the **root** of the tree is rounded to `p`; internal `None` rounding plus `errs[]` is the cancellation heuristic in `doc/README.md`.
- Unbounded exponent: set `emin`/`emax` to the real range or precision (and time) grow with unused magnitude.

## Practical rule

Put the formula in `expr!`. Use a `Context` whose `p` is the output precision, `rm` is the required mode (usually `ToEven`), and `emin`/`emax` bound the formula class.
