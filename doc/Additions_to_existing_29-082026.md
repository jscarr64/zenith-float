# Review of API/docs recommendations (2026-08-29)

The original note started at item **2** (item **1** was not in the pasted file). Items **2, 5, 6, 7, 9, 10** are in the crate. **Complex expressions use `cexpr!`** — this crate has macros; the Accumath “no macros” rule does not apply here. Items **4** and **8** remain load-bearing exclusions.

## Implemented (original 2, 5, 6, 7, 9, 10, and complex macros)

### Complex `cexpr!`

Same working-precision loop as `expr!`, `ExactComplex` leaves, cancellation on both parts. Imaginary unit `I`. This is a real expander, not a method-wrapping stub.

### 2. ExactComplex method coverage

`tan`, `sinh`/`cosh`/`tanh`, `sqrt`, `pow`, `asin`/`acos`/`atan`, `asinh`/`acosh`/`atanh` (principal branches). Call these directly. There is no complex `expr!`.

### 5. Public EFT primitives

`ExactNum::two_sum` and `two_product` (exact limb sum/product, then a `p`-bit high part). Reconstruct with `hi.add(&lo, p, rm)`.

### 6. `fused_sum` / `fused_dot`

Extra working bits, one final round. Length mismatch on `fused_dot` → NaN (`InvalidArgument`).

### 7. Horner `polyval`

`ExactNum::polyval(coeffs, x, p, rm)` — `a0, a1, …` via `fma`.

### 9. `Context::with_rounding_mode`

Restores the previous mode when the closure returns. Panic leaves the temporary mode.

### 10. Getting-started doc

[`GETTING_STARTED.md`](GETTING_STARTED.md), linked from the crate README.

---

## Load-bearing exclusions (original 4, 8)

### 3. Complex expressions — **done as `cexpr!`**

Not blocked by macros. `cexpr!` uses the same extra-precision loop and measures cancellation on both parts.

### 4. Bessel Y_n, I_n, K_n, fractional order — **will not add here**

Integer `J_n` is a factorial power series (`n ≤ 1024`). `Y_n` is a different numerical project (log + second series; recurrence is ill-conditioned in the wrong direction), with its own MPFR harness and cost model. Half-shipping it would undercut the rest of the library. Not a wrapper around `bessel_j`.

### 8. In-tree hardware IEEE conversion — **will not add**

A feature-gated module would force a CI exception, hedged docs, and an invariant of “no hardware float except.” Plotting/FFI belongs in a **separate crate** that depends on zenith-float and uses `frexp` / `ilogb`; the precision contract stays obvious.

---

## Already correctly excluded

- **`Hash` / total `Ord`:** NaN makes these wrong.
- **Assigning operators (`+=`):** explicit `(p, rm)` is the API.
