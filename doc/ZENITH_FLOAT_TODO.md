# zenith-float — Special Functions To-Do

This document lists the special functions that need to be added to `zenith-float` to support Accumath's build plan. Every item here blocks one or more Accumath capabilities. Nothing ships on either side until the dependency chain is complete.

All implementations follow the same rules as existing specials (`erf`, `gamma`, `bessel_j`, `besseli`):

- Software limb arithmetic only — no hardware float in any calculation
- Series or AGM with working precision `p_wrk = p + WORD_BIT_SIZE`; `MAX_PREC_RETRY` bounds retries
- MPFR oracle golds (`mpfr-tests` feature) on bounded domains
- Per-op precision at caller-chosen bits (not global `SOFT_PREC`) — this is the same fix needed for unary specials generally
- `expr!` leaf entry where applicable
- Derivative identity golded alongside the function

---

## Priority 1 — Integral result specials

These unblock Accumath item 2 (named non-elementary integration results) and items 5–6 (Laplace / Fourier numeric eval). **Done 2026-08-29:** `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c` on `ExactNum` with `expr!` leaves. SoftFloat calls them. Domain: `Ei`/`Ci` require `x>0`; `li` requires `x>1`. Large `|x|` uses the full factorial / auxiliary \(f,g\) expansions (not a one-term remainder). Golds: Accumath `eval.gold` (`si_zero`, `li(e)/Ei(1)`, Fresnel zeros), specials series identities, and two-precision large-argument checks. MPFR oracles for these leaves are still a follow-up.

| Function | Definition | Notes |
| --- | --- | --- |
| `Si(x)` — sine integral | `∫₀ˣ sin(t)/t dt` | Series for small `\|x\|`; asymptotic for large; `Si(-x) = -Si(x)` |
| `Ci(x)` — cosine integral | `γ + ln x + ∫₀ˣ (cos t - 1)/t dt` | `x > 0` only; `Ci(x)` for `x ≤ 0` is `NaN` |
| `li(x)` — logarithmic integral | `PV ∫₀ˣ dt/ln t` | `x > 0`, `x ≠ 1`; relation to `Ei`: `li(x) = Ei(ln x)` |
| `fresnel_s(x)` | `∫₀ˣ sin(πt²/2) dt` | Series for small `\|x\|`; auxiliary functions `f`, `g` for large |
| `fresnel_c(x)` | `∫₀ˣ cos(πt²/2) dt` | Same auxiliary function strategy as `fresnel_s` |

**Golds required for each:**
- Value at a known point vs MPFR oracle
- Derivative identity (e.g. `D(Si)(x) = sin(x)/x`)
- Asymptotic behavior gold (large `\|x\|`)
- `NaN` / domain guard gold

**`expr!` leaves:** `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`

---

## Priority 2 — Modified Bessel second kind

Unblocks Accumath item 8. `J_ν`, `Y_ν`, `I_ν` already exist; `K_ν` does not.

| Function | Definition | Notes |
| --- | --- | --- |
| `K_ν(x)` — modified Bessel second kind | `(π/2)(I_{-ν} - I_ν)/sin(νπ)` for non-integer ν; limit for integer ν | `x > 0` only; `K_{-n} = K_n` for integer `n`; exponential decay for large `x` |

**Algorithm:** For integer `n`, use the recurrence `K_{n+1} = (2n/x)K_n + K_{n-1}` starting from `K_0` and `K_1` computed by series. For large `x`, use the asymptotic expansion `K_ν(x) ~ √(π/2x) e^{-x} Σ`. For half-integer orders, use the closed form via `sinh`/`cosh`.

**Golds required:**
- `K_0(1)` vs MPFR oracle
- `K_1(1)` vs MPFR oracle
- `K_{1/2}(x) = √(π/2x) e^{-x}` — closed form gold
- `K_{-n} = K_n` symmetry gold
- `D(K_0) = -K_1` derivative identity gold
- `x ≤ 0` → `NaN` guard gold

**`expr!` leaf:** `bessel_k`

---

## Priority 3 — Elliptic integrals

Unblocks Accumath item 1. Nothing exists today — `elliptic_f` / `elliptic_e` appear only as unevaluated derivative skeletons.

| Function | Definition | Notes |
| --- | --- | --- |
| `K(k)` — complete elliptic first kind | `∫₀^{π/2} dθ/√(1-k²sin²θ)` | AGM algorithm; `\|k\| < 1`; `k = ±1` → `+∞` |
| `E(k)` — complete elliptic second kind | `∫₀^{π/2} √(1-k²sin²θ) dθ` | AGM-based; `\|k\| ≤ 1`; `E(0) = π/2`; `E(1) = 1` |
| `Π(n,k)` — complete elliptic third kind | `∫₀^{π/2} dθ/((1-n sin²θ)√(1-k²sin²θ))` | More complex AGM variant; `n < 1`, `\|k\| < 1` |
| `F(φ,k)` — incomplete elliptic first kind | `∫₀^φ dθ/√(1-k²sin²θ)` | `F(π/2, k) = K(k)` |
| `E(φ,k)` — incomplete elliptic second kind | `∫₀^φ √(1-k²sin²θ) dθ` | `E(π/2, k) = E(k)` |
| `Π(n,φ,k)` — incomplete elliptic third kind | `∫₀^φ dθ/((1-n sin²θ)√(1-k²sin²θ))` | `Π(n,π/2,k) = Π(n,k)` |

**Algorithm:** AGM (arithmetic-geometric mean) for complete forms. Descending Landen transformation for incomplete forms. Both converge quadratically in the number of AGM steps — well-suited to arbitrary precision.

**Golds required:**
- `K(0) = π/2` — exact gold
- `K(1/√2)` vs MPFR oracle
- `E(0) = π/2`, `E(1) = 1` — exact golds
- `F(π/4, 1/√2)` vs MPFR oracle
- Legendre relation `E(k)K'(k) + E'(k)K(k) - K(k)K'(k) = π/2` — identity gold
- `k ≥ 1` → `NaN` guard for `K` and incomplete forms
- Derivative identities: `dK/dk`, `dE/dk` in terms of `K` and `E`

**`expr!` leaves:** `elliptic_k`, `elliptic_e`, `elliptic_pi`, `elliptic_f`, `elliptic_e_inc`, `elliptic_pi_inc`

---

## Priority 4 — Legendre polynomials and associated functions

Unblocks Accumath items 9 and 10 (spherical harmonics depends on associated Legendre).

| Function | Definition | Notes |
| --- | --- | --- |
| `P_n(x)` — Legendre polynomial | Three-term recurrence `(n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}` | `n ≥ 0` integer; `\|x\| ≤ 1` for orthogonality; exact for integer `x` via recurrence |
| `P_n^m(x)` — associated Legendre | `P_n^m(x) = (-1)^m(1-x²)^{m/2} d^m/dx^m P_n(x)` | `0 ≤ m ≤ n`; Condon-Shortley phase `(-1)^m` included; `m > n` → 0 |

**Algorithm:** Three-term recurrence for `P_n`; forward recurrence in `m` for `P_n^m` starting from `P_m^m` and `P_{m+1}^m`.

**Named cap:** `LEGENDRE_N_MAX` — maximum degree `n`; return `NaN` above this. Set to a value where the recurrence remains numerically stable at the working precision.

**Golds required:**
- `P_0(x) = 1`, `P_1(x) = x`, `P_2(x) = (3x²-1)/2` — exact golds
- `P_5(0.5)` vs MPFR oracle
- `P_2^1(x) = -3x√(1-x²)` — exact gold
- Orthogonality: `∫₋₁¹ P_m P_n dx = 2/(2n+1) δ_{mn}` — numeric gold at 256 bits
- `D(P_n) = nP_{n-1} + xD(P_{n-1})` — derivative identity gold
- `m > n` → 0 gold
- `n > LEGENDRE_N_MAX` → `NaN` gold

**`expr!` leaves:** `legendre_p`, `legendre_p_assoc`

---

## Priority 5 — Gauss hypergeometric function

Unblocks Accumath item 11. This is the most complex item; series convergence is not guaranteed for all `(a,b,c,z)`.

| Function | Definition | Notes |
| --- | --- | --- |
| `₂F₁(a,b;c;z)` | `Σ_{n=0}^∞ (a)_n(b)_n/(c)_n · z^n/n!` | `(x)_n` is the Pochhammer symbol; `c` not a non-positive integer; `\|z\| < 1` for series |

**Algorithm:**
- Series for `\|z\| < 1` with convergence cap `HYPERGEOM_SERIES_MAX_TERMS`
- Euler transformation `₂F₁(a,b;c;z) = (1-z)^{c-a-b} ₂F₁(c-a,c-b;c;z)` to extend to `\|z\| < 1` from the other side
- Pfaff transformation `₂F₁(a,b;c;z) = (1-z)^{-a} ₂F₁(a,c-b;c;z/(z-1))` for `Re(z) < 1/2`
- Kummer transformation for `z = 1` when `Re(c-a-b) > 0`
- `Ok(None)` / `NaN` when no transformation brings `z` into the convergence region within the cap

**Named caps:** `HYPERGEOM_SERIES_MAX_TERMS` — maximum series terms before abandoning; `HYPERGEOM_TRANSFORM_MAX` — maximum transformation attempts.

**Golds required:**
- `₂F₁(1,1;2;z) = -ln(1-z)/z` — exact identity gold
- `₂F₁(1/2,1/2;1;k²) = (2/π)K(k)` — connection to elliptic `K` gold
- `₂F₁(a,b;c;0) = 1` — exact gold
- `₂F₁(a,b;c;1) = Γ(c)Γ(c-a-b)/(Γ(c-a)Γ(c-b))` when convergent — Gauss evaluation gold
- `c` a non-positive integer → `NaN` gold
- `\|z\| ≥ 1` with no applicable transformation → `NaN` / `Unsupported` gold

**`expr!` leaf:** `hypergeom_2f1`

---

## Summary table

| Function(s) | Priority | Blocks Accumath | Status |
| --- | --- | --- | --- |
| `Ei`, `Si`, `Ci`, `li`, `fresnel_s`, `fresnel_c` | 1 | Items 2, 5, 6 | ✅ 2026-08-29 |
| `K_ν` modified Bessel second kind | 2 | Item 8 | ⬜ Not implemented |
| `K(k)`, `E(k)`, `Π(n,k)` complete elliptic | 3 | Item 1 | ⬜ Not implemented |
| `F(φ,k)`, `E(φ,k)`, `Π(n,φ,k)` incomplete elliptic | 3 | Item 1 | ⬜ Not implemented |
| `P_n(x)` Legendre | 4 | Items 9, 10 | ⬜ Not implemented |
| `P_n^m(x)` associated Legendre | 4 | Items 9, 10 | ⬜ Not implemented |
| `₂F₁(a,b;c;z)` hypergeometric | 5 | Item 11 | ⬜ Not implemented |

---

## Cross-cutting requirements for all new specials

These apply to every item above without exception:

- No hardware float in implementation or tests
- Per-op precision at caller-chosen bits — not global `SOFT_PREC`
- `expr!` leaf added and documented in `EXPR.md` table
- `LIBRARY.md` §12 updated with the new method
- MPFR oracle gold in `mpfr-tests` feature
- Derivative identity golded
- Domain guard golded (`NaN` for out-of-domain input)
- `ZENITH_FLOAT_CAPABILITIES.md` §12 updated
- `scripts/ci_full.sh` green before marking done
