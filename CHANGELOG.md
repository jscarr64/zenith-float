# Changelog

## Unreleased

## 1.0.10 — 2026-10-05

### Gauss–Kronrod (7, 15) and Wynn ε

SoftFloat quadrature for `∫_a^b f(x) dx` on a finite interval. No `f64` path.

- `kronrod_pair` builds the (7, 15) nodes and both weight rows at `p + WORD_BIT_SIZE`. The recurrence is Laurie’s Jacobi–Kronrod matrix for the Legendre weight (*Math. Comp.* 66, 1997); weights are Golub–Welsch. Orders other than 7 return `None`.
- `gauss_kronrod_interval` returns the Kronrod sum `K`, the embedded Gauss sum `G`, and `abs_err_est = |K − G|`.
- `wynn_epsilon` is Wynn’s ε table on a sequence of partial sums.
- `integrate_adaptive_gk(f, a, b, tol_rel, max_subintervals, p, rm, cc)` bisects the panel with the largest `|K − G|` and, when that raw sum is still slow, accepts a stable ε value. The absolute stop is `max(tol_rel * max(1, |estimate|), 2^(1−p))`. Failure is [`QuadratureError`], not a panic and not a quiet wrong number.

The local error `|K − G|` and the ε acceleration for endpoint singularities are the QUADPACK design of Piessens, de Doncker-Kapenga, Überhuber, and Kahaner (*QUADPACK*, Springer, 1983; SLATEC’s port of that design is public domain). This tree recomputes the rule in SoftFloat `ExactNum`. It does not copy the ACM / SLATEC Fortran.

Still unsupported: infinite intervals (`QAGI`), oscillatory weights (`QAWO` / `QAWF`), algebraic endpoint weights (`QAWS`), Cauchy principal values (`QAWC`), Kronrod pairs other than (7, 15), and the full `QAGS` heap.

## 1.0.9 — 2026-10-04

The complex \(I_\nu\) power-series guard is sized to cancellation, not \(1.5\lvert z\rvert+16\) (that blanket is what 1.0.6 through 1.0.8 shipped; 1.0.7 and 1.0.8 were not this fix). When every term is positive (real \(z>0\), real \(\nu>-1\)) the guard is a fixed 32-bit pad. Otherwise it is that pad plus an upper bound on \(\lceil(\lvert z\rvert-\mathrm{Re}\,z)/\ln 2\rceil\), raised again if the summed peak still exceeds the result. Real \(I_\nu\) for large \(\lvert x\rvert\) now runs past the peak at \(k\approx\lvert x\rvert/2\). Real \(K_\nu\) uses the asymptotic series when the \(2.885\lvert x\rvert\) guard would exceed `BESSEL_GUARD_MAX`, instead of returning `InvalidArgument`. Complex \(K\) on the positive real axis uses that real kernel past the series switch. Clippy `-D warnings` cleanups in the same tree (manual `find`, `clamp`, `from_ref`, single-iteration test loops) do not change numeric results. The 1.85 `unused_unsafe` allow on `_addcarry_u64` / `_subborrow_u64` stays.

## 1.0.8 — 2026-10-03

README upgrade target corrected to 1.0.8; no math.

## 1.0.7 — 2026-10-02

1.0.7 keeps the 1.85-required unsafe wrappers on `_addcarry_u64`/`_subborrow_u64` and allows `unused_unsafe` so the crate builds on rustc 1.93+ (those intrinsics became safe in 1.93, not 1.87).

## 1.0.6 — 2026-10-02

Additive SoftFloat pack: Jacobi \(P_n^{(\alpha,\beta)}\), real large-\(|x|\) Bessel Hankel, complex Bessel cut/speed fixes, and the catalog specials Accumath was missing. No public signature was removed or changed — **1.0.6, not 2.0.0**. Jeff publishes crates.io — this tag is not published from the PR.

### Version path

Stayed on **1.0.6**. Every new entry is an added `ExactNum` / `ExactComplex` method (or a bugfix of an existing one). No existing signature changed, so a breaking 2.0.0 bump is not required.

### `ExactNum::jacobi_p`

- Jacobi \(P_n^{(\alpha,\beta)}(x)\) via the Bonnet three-term recurrence, same house style as `hermite_h` / `chebyshev_t` in `orthopoly.rs`. Cap `ORTHOPOLY_N_MAX`.
- Golds: \(\alpha=\beta=0\) reduces to Legendre \(P_n\); \(P_2^{(1,1)}(1/2)=3/16\); \(P_n^{(\alpha,\beta)}(-x)=(-1)^n P_n^{(\beta,\alpha)}(x)\).
- Independent mpmath 1.4.1 references for non-integer \((\alpha,\beta)\), including \(\alpha=1/3\), \(\beta=2/3\).

### Real Bessel large \(|x|\)

- `bessel_j` / `bessel_j_nu` / `bessel_y` use the Hankel expansion (DLMF 10.17, \(P,Q\) stopped at the smallest term) when \(|x|\ge\max(16,\,0.35(p+112))\) and \(|\nu|^2<2|x|\), and whenever the series cancellation guard would exceed `BESSEL_GUARD_MAX` under that order condition.
- Closes the 1.0.5 known issue: real \(J_0(20000)\), \(Y_0(20000)\) and the rest of that class no longer return `NaN(InvalidArgument)`.
- Negative \(x\) for integer order uses \(J_n(-x)=(-1)^n J_n(x)\). \(Y\) remains defined for \(x>0\).
- Regression tests vs mpmath at 128 and 256 bits, plus a Wronskian gold at \(x=20000\).

### Complex Bessel (negative real axis + \(K\) speed)

- On the cut \(\operatorname{Im} z=0\), \(\operatorname{Re} z<0\) with real \(\nu\), \(J,Y,I,K\) are assembled from the real kernels (DLMF 10.11 / 10.34) so algebraically zero parts stay exact zeros. Closes the 1.0.5 hang / `PrecisionRetryExhausted` on e.g. `bessel_i(-2, 1/2)`.
- \(I_\nu(z)\) for \(\mathrm{Re}\,z>0\) uses the well-conditioned \({}_0F_1\) series instead of a rotation through \(J\).
- \(K_\nu(z)\) for \(\mathrm{Re}\,z>0\) past the series/Hankel switch uses Temme \({}_2F_0\) (DLMF 10.32.10, modified Lentz) instead of \(J\pm iY\) with \(3|z|\) guard bits. The Hankel / series path remains the fallback. High-precision non-integer \(K\) just inside the series regime can still be slow; the Temme path covers the large-\(|z|\) cases that were the 74 s class.

### Catalog specials (additive)

New `ExactNum` methods, SoftFloat-only, with identity and mpmath golds:

- Scorer \(\mathrm{Gi},\mathrm{Hi}\) (`scorer_gi`, `scorer_hi`); \(\mathrm{Hi}=\mathrm{Bi}-\mathrm{Gi}\).
- Kelvin \(\mathrm{ber},\mathrm{bei},\mathrm{ker},\mathrm{kei}\) (order 0 and `_nu`).
- Struve \(\mathbf{H}_\nu\), Anger \(\mathbf{J}_\nu\), Weber \(\mathbf{E}_\nu\).
- Clausen \(\mathrm{Cl}_2,\mathrm{Cl}_3\); Barnes \(G\); Hurwitz / Riemann \(\zeta\) (integer \(s\ge 2\)); polygamma \(\psi^{(n)}\).
- Inverse Jacobi `jacobi_arcsn` / `arccn` / `arcdn`; \({}_0F_1\), \({}_1F_1\), \({}_pF_q\); Lambert \(W_0,W_{-1}\); \(\mathrm{Li}_n\) for \(|x|\le 1\).

`expr!` leaves for the one- and two-argument catalog names (`scorer_gi`, `kelvin_ber`, `struve_h`, `clausen_cl2`, `barnes_g`, `lambert_w0`, …).

### Meijer \(G\) / Fox \(H\) / Heun (additive, still 1.0.6)

Leftovers after the catalog pack — Accumath still had local SoftFloat clones for these. No existing signature changed.

- `ExactNum::meijer_g(an, ap, bm, bq)` — DLMF 16.17 residue + \({}_pF_q\). Type-1 or type-2 series by \((p,q,\lvert z\rvert)\). Coincident poles / numerator \(\Gamma\) poles → `NaN` (no logarithmic `hypercomb` limit). Negative \(z\) → `NaN` unless every used power is an integer.
- `ExactNum::fox_h` — pairs \((a,A)\), \((b,B)\). \(A=B=1\) is Meijer \(G\); positive rational scales use the Gauss multiplication formula (mpmath `foxh`). Non-rational or large scales → `NaN`.
- `ExactNum::heun_g(a,q,α,β,γ,δ)` — local Heun \(\mathrm{Hl}\), DLMF 31.3, \(\lvert z\rvert<1\).
- `ExactNum::heun_c(α,γ,δ,ε,q)` — confluent Heun, \(y(0)=1\) Frobenius solution of DLMF 31.12.1, \(\lvert z\rvert<1\). Not Maple `HeunC(α,β,γ,δ,η,z)`.
- Left out (no Accumath SoftFloat evidence, or no unique \(z=0\) series): Mathieu, Lamé, Appell \(F_{1..4}\), HeunB / HeunD / HeunT, complex \(G/H\)/Heun. `expr!` skipped (too many arguments, same as `hypergeom_pfq`).
- Golds: mpmath `meijerg` / `foxh` at 50 digits (including \(J_{1/4}\), series-2 \(\lvert z\rvert>1\), rational Fox scale); HeunG vs mpmath \({}_2F_1\); HeunC vs mpmath \({}_1F_1\).

## 1.0.5 — 2026-09-27

Bug-fix release: special-function accuracy. No public API changes (no new public items, no signature changes); results that were wrong now match mpmath to working precision. Every fixed case below is covered by a test with a 50-digit (or 2p+64-bit) mpmath reference.

### `hypergeom_2f1` (real)

- **What was wrong.** For \(1/2<z<1\) the Pfaff transform maps to \(\lvert w\rvert>1\), and the divergent series was summed there and silently truncated at the term cap. Example: \({}_2F_1(1,1;2;0.7)\) returned `6.06e175` (true `1.71996…` = \(-\ln(0.3)/0.7\)). Near \(z=-1\) the series was truncated too. For \(z\le-1\) the A&S 15.3.7 path returned NaN whenever \(b-a\) is an integer (e.g. \({}_2F_1(0.5,1.5;2;-3)\)). The series stop was absolute, not relative.
- **What changed.** New dispatch:
  - poles in \(c\) → NaN;
  - terminating (polynomial) series;
  - Gauss at \(z=1\) when \(c-a-b>0\);
  - \(z>1\) → NaN;
  - Pfaff for \(z<-1/2\);
  - direct series for \(\lvert z\rvert\le1/2\);
  - \(z\to1-z\) connection for \(1/2<z<1\) (DLMF 15.8.4), with the A&S 15.3.10–15.3.11 logarithmic case for integer \(c-a-b\).

  Each path measures cancellation and re-evaluates with that many guard bits. The series returns an error at its term cap instead of truncating.

### `betainc`

- **What was wrong.** The integer-\(b\) polynomial and the series were used at large \(a,b\), where they cancel catastrophically. Examples: \(I_{0.7}(30.5,2.25)\) returned `2.68e240` (true `3.0312…e-4`), \(I_{0.3}(1000,1000)\) returned `1.47e124`, and \(I_{0.5}(5000,5000)\) returned `1.8e2284`.
- **What changed.** A continued fraction (DLMF 8.17.22, modified Lentz), with the symmetry \(x>(a+1)/(a+b+2)\) applied. The polynomial is used only when its measured loss is ≤ 64 bits. \(I_{1/2}(a,a)=1/2\) is returned exactly.

### Other real special functions

- **`erf` / `erfc` / `normal_cdf`.** For every \(\lvert x\rvert\ge4\) (binary exponent ≥ 3), 1.0.2–1.0.4 summed the divergent asymptotic series past its smallest term, so the tail was garbage at every precision. At 256 bits: `erfc(4.2)` = `1.98e349` (true `2.855e-9`), `erfc(10)` = `4.99e17` (true `2.088e-45`), `erfc(12.5)` = `7.69e-56` (true `6.232e-70`, 1.2e14 times too large), `erf(4)` = `-4.40e366`, `erf(-4.2)` = `1.98e349`, `normal_cdf(-6)` = `1.53e433`. At 128 bits `erf(4)` = `-1.22e199` and `erfc(10)` = `1.3e-18`. 1.0.2 and 1.0.4 give identical wrong values (the code did not change between them). Now the asymptotic series is used only when \(x^2\ge0.7p+4\), where its smallest term is below \(2^{-p}\), and it stops at its smallest term. Otherwise the Maclaurin series runs with \(1.5x^2\) guard bits, and `erfc` for \(x>0\) carries another \(1.5x^2\) bits for the cancellation in \(1-\operatorname{erf}\). `normal_cdf` uses `erfc(-(x-μ)/(σ√2))/2`, so the lower tail keeps full relative accuracy (`normal_cdf(-40)` was garbage). A sweep of 561 points (x from −30 to 30, step 0.1 on \(3.5\le\lvert x\rvert\le30\)) for `erf`, `erfc` and `normal_cdf` at 64, 128, 192, 256, 512 and 1024 bits is within 2 bits of mpmath everywhere. 1.0.4 misses at 619 of these points at 128 bits.
- **`gamma` / `digamma`.** The Stirling start was fixed, which capped accuracy at about 274 bits (Γ) and 530 bits (ψ) at any precision. The shift is now precision-dependent. `gamma(11.5)` at 1024 bits had 276 correct bits.
- **Bernoulli numbers (`gamma`, `ln_gamma`, `digamma`, real and complex).** The Stirling coefficients came from the Akiyama–Tanigawa recurrence in floating point, which cancels catastrophically. With 64-bit limbs the Ziv loop usually hid this. With 32-bit limbs (`thumbv7em-none-eabihf`, i686) Γ at an internal 192 bits had about 101 correct bits, so `bessel_i(2.25, 1.5)` at 128 bits returned only 101 correct bits. The \(B_{2k}\) are now built from exact tangent numbers (Brent–Harvey) with one rounding each. They are computed once per evaluation instead of once per term (the old cost was cubic in the number of terms).
- **`ln_gamma`.** `ln_gamma(1e20)` and `ln_gamma(1e300)` returned `inf`; large arguments now use Stirling directly.
- **Euler–Mascheroni γ.** The constant (used by `ei`, `ci`, `li`, `bessel_y`, `bessel_k`, `digamma`) was correct to only about 383 bits at any precision. It is now computed with Brent–McMillan.
- **`ei`.** Series cancellation for \(x<0\) was unguarded: `ei(-50)` had 115 of 128 bits and `ei(-200)` was garbage. Guard bits are now sized from the measured loss.
- **`fresnel_s` / `fresnel_c`.** Series cancellation was unguarded (`fresnel_s(10)` had 36 bits, `fresnel_c(12)` was garbage), and the asymptotic switch now depends on precision.
- **`ai` / `bi` / `ai_prime` / `bi_prime`.**
  - The derivative asymptotic coefficients had the wrong sign/form, so `ai_prime(8)` and `bi_prime(-30)` had 7 correct bits.
  - The series/asymptotic switch now depends on precision: `ai(12)` had 85 bits, and `ai(-30)` had 320 bits at p=512.
- **Real Bessel.**
  - `bessel_j_nu` / `bessel_y` / the integer series had no cancellation guard, so J/Y at \(x=200\) were garbage. They now carry \(1.5\lvert x\rvert\) guard bits. For \(\lvert x\rvert\gtrsim10^4\) they return `InvalidArgument` rather than garbage (no large-argument Hankel expansion yet).
  - `bessel_k` used its asymptotic series below the precision where it converges (`bessel_k(20, 0)` had 62 bits).
  - Series stops are now relative: `bessel_j_nu(1e-10, 10.5)` had 276 of 512 bits.
- **Series stopping rule.** The absolute stopping rule was replaced by a relative one in `gammainc`, `si`, Fresnel and the Bessel series.
- **`elliptic_f(±1 | 1)`** returned `NaN(InvalidArgument)` even though \(\lvert x\rvert\le1\) is the documented domain. It now returns \(\pm\infty\), consistent with `elliptic_k(1)` = \(+\infty\).
- **Carlson \(R_C\), `elliptic_pi`, `elliptic_pi_complete`.** The \(x\approx y\) shortcut used an absolute threshold of \(2^{-(p/3+16)}\), which limited Π to about 340 bits at p=512. The threshold is now relative. The duplication cap scales with precision (`max(128, p/4+64)`), for real and complex.
- **Exactly representable results.** Ziv's loop could not certify these and returned `PrecisionRetryExhausted`: `log(100, 10)`, `log(2, 4)`, `log(0.25, 2)`, `nth_root(1e10, 5)`, `jacobi_cd(u, 1)`, `jacobi_dc(u, 1)`. They are now detected and returned exactly.

### Complex (`ExactComplex`)

- **`erf` / `erfc` / `fresnel_s` / `fresnel_c`.** The power series had no cancellation guard. At 512 bits, `erf(-15-5i)` was wrong from the 45th digit (true value is \(-1+4.4\times10^{-89}+2.3\times10^{-89}i\)) and `fresnel_s(-15-5i)` had 200 bits. The series now carries \(1.5\lvert u\rvert^2\) guard bits.
- **`ci`.** Wrong overall sign/branch: `ci(0.5+0.5i)` and `ci(-2+i)` had 0 correct bits. It is now principal and matches mpmath.
- **`si(i)`.** Returned NaN.
- **`ei`.** The \(i\pi\,\mathrm{sgn}(\operatorname{Im} z)\) term was missing from the asymptotic path, and the switch point was fixed rather than precision-dependent (`ei(20+20i)` had 22 bits).
- **`acosh` / `atanh` on their cuts.** `acosh(-2+0i)` had a negative real part and `acosh(-1+0i)` returned \(-i\pi\). `atanh(±1+0i)` was NaN; it is now \(\pm\infty+0i\).
- **`ai` / `bi`.**
  - The phases in the `bi` connection formula were swapped: `bi(0.001+10i)` had 0 bits.
  - `bi` on the Stokes line used only the dominant asymptotic series (`bi(20+20i)` had 110 bits).
  - The switch point now depends on precision (`ai(0.001+10i)` had 63 bits).
- **`bessel_*`.**
  - The Hankel P/Q sums ran past their smallest term, giving garbage at moderate \(\lvert z\rvert\).
  - `bessel_i` / `bessel_k` used the wrong rotation for \(\operatorname{Re} z<0\): `bessel_i(-2+i, 2+0.5i)` and `bessel_k(-2+i, 0)` were wrong. They now use DLMF 10.27.6 / 10.27.8.
  - The series now has cancellation guards.
  - `bessel_j_nu` / `bessel_y` with \(\operatorname{Re} z<0\) in the Hankel regime (\(\lvert z\rvert\gtrsim 0.35(p+112)\)) summed the large-argument expansion outside its sector \(\lvert\arg z\rvert<\pi\): `bessel_j_nu(-20000, 0)`, `bessel_j_nu(-500+i, 1)` and `bessel_y(-300+40i, 2.5)` had 0 correct bits (1.0.4 too). They now continue from \(-z\) (DLMF 10.11.1–10.11.2).
  - `bessel_i` / `bessel_k` at real \(z>0\) with real non-integer \(\nu\) never converged: the imaginary part is exactly 0 but was computed as a rounding residue that Ziv's loop cannot certify. 1.0.4 returned `PrecisionRetryExhausted` after about 3 s for `bessel_i(2, 0.5)`; the first 1.0.5 candidate took about 5 minutes and returned NaN. The imaginary part is now an exact 0 there.
  - `bessel_k` on its cut \((-\infty,0)\) with \(\operatorname{Im} z=0\) returned the value from below; it now returns the value from above, like `bessel_i`, `bessel_j_nu` and mpmath.
- **`gamma` / `ln_gamma` / `digamma`.** Precision was capped at about 275 bits. `ln_gamma` now returns the documented principal \(\ln(\Gamma(z))\): `ln_gamma(3-4i)` returned a different branch.
- **`hypergeom_2f1`.**
  - The documented 10 000-term cap was actually \(p+96\), and the series stop was absolute. `2F1(0.5,0.5;1;0.999999)` had 1 correct bit.
  - All-real arguments with \(z<1\) now use the real algorithm.
  - Integer \(m=c-a-b\) with \(z\) near 1 (\(\lvert 1-z\rvert\le 1/2\), or \(\lvert 1-z\rvert<1\) with \(\lvert z\rvert\ge1\)) returned NaN: the direct series converges too slowly there (or not at all for \(\lvert z\rvert\ge1\)) and the generic \(1-z\) transform has \(\Gamma(\pm m)\) poles. Such arguments now use the logarithmic connection (DLMF 15.8.10, A&S 15.3.10–15.3.11) for \(m\ge0\), after Euler's transformation \((1-z)^m F(c-a,c-b;c;z)\) for \(m<0\); cancellation is measured and the evaluation repeated with more bits. Example: \({}_2F_1(0.3+0.2i,\,0.7;\,2+0.2i;\,1.2+0.3i)\).
  - \(c-a-b\) within \(\delta\) of an integer (the generic \(1-z\) transform loses about \(\log_2(1/\delta)\) bits) now carries that many guard bits; \(z\) close to 1 inside the unit disc with non-integer \(c-a-b\) uses the transform instead of a slowly converging series.
  - A grid of 504 cases (9 \((a,b)\) pairs, \(m\in\{-3,\ldots,3\}\) plus near-integer \(c\), 9 values of \(z\) around 1 including \(\lvert z\rvert>1\)) agrees with mpmath to within 2 bits at 64, 128, 256 and 512 bits.

### `no_std`: built and tested

- The library already built without `std`. Now the reference tests also run without it.
- `nostd-tests/` is a `#![no_std]` crate (core + alloc only; not published; excluded from the workspace and the package) with 77 mpmath reference cases covering:
  - elementary functions, erf/erfc (including the tail), Γ/lnΓ/ψ at 128–320 bits, 2F1, betainc, Ei/Si/Ci, Fresnel, Airy, Bessel, elliptic, Jacobi, `normal_cdf`;
  - complex Γ/erf/Ci/Ai/J/I/2F1 (including the integer-\(c-a-b\) log case near \(z=1\), \(J\) with \(\operatorname{Re} z<0\) and \(I\) on the real axis);
  - quadrature, Brent and RK45.
- Two runners:
  - `nostd-host`: a `#![no_std]`/`#![no_main]` x86_64 Linux binary linking only libc;
  - `nostd-qemu`: bare-metal `thumbv7em-none-eabihf` (Cortex-M4F, 32-bit limbs) under `qemu-system-arm -machine mps2-an386` with semihosting.
- `scripts/ci_nostd.sh` builds `zenith-float` and `zenith-float-num` with `--no-default-features --target thumbv7em-none-eabihf` and runs both runners. Each passes 77/77.
- There is no separate `alloc` feature; an allocator is always required.

### Comparison, and the numeric methods that depend on it

- **`ExactNum::cmp` / `ExactNum::abs_cmp` return exactly `-1`, `0` or `1`.** The documented contract (positive / negative / zero) is unchanged.
  - For finite values with equal exponents, 1.0.4 returned the raw mantissa word difference: `3.cmp(2)` was `Some(2^62)`.
  - The crate's own numeric routines test `cmp(..) == Some(-1)` / `== Some(1)`, so in 1.0.4:
    - `gauss_legendre`, `tanh_sinh`, `bisect`, `illinois`, `brent`, `rk4`, `euler`, `rk45_adaptive`, `chebyshev_coeffs` and `chebyshev_eval` rejected (`None` / NaN) any interval whose endpoints share a binary exponent, e.g. `[2, 3]`;
    - `rk45_adaptive` with a non-zero `rtol` never accepted a step and returned `None`;
    - several convergence tests misfired.
- **`ExactNum::abs_cmp` compared signed values** for finite operands, contradicting its documentation: `(-2).abs_cmp(1)` was negative. It now compares magnitudes.
- **`brent` rewritten** as Brent's classic zero finder (bracket `[b, c]`, inverse quadratic / secant with a bisection fallback). The old loop had no bisection safeguard: it stalled like regula falsi and returned `None` with `root_default_tol` at 128 bits.
- **Resolution stop in `bisect` / `illinois` / `brent`.** When `tol` is below the working resolution (e.g. `root_default_tol` = \(2^{-256}\) at \(p<256\)), they now return the bracketed root instead of `None` after `ROOT_MAX_ITER`.

### Tests

- `ops::special` unit tests:
  - `test_hypergeom_2f1_mpmath`: 41 values across \(z\in[-10,1]\), including the log case, \(z=1\) and near-1 points, plus 4 error cases;
  - `test_betainc_mpmath`: 17 values;
  - `test_gamma_digamma_high_precision`: 1024 bits.
- New integration test `tests/special_audit_mpmath.rs`: 61 regression points (real and complex, 128/512/1024 bits). Each must agree with mpmath to within 2 ulp; 53 of them fail on 1.0.4.
- Reference generators: `refs/` in the maintainer workspace (mpmath 1.4.1, 50 digits or 2p+64 bits).
- New integration test `tests/numeric_methods_regress.rs`: 5 tests covering `cmp` normalization, `[2, 3]` intervals for quadrature, roots, ODE and Chebyshev, `rk45_adaptive` with `rtol`, and `brent` / `bisect` / `illinois` with `root_default_tol` at 128 bits. All 5 fail on 1.0.4.
- New integration test `tests/erf_tail_mpmath.rs`: 221 `erf` / `erfc` / `normal_cdf` cases with 50-digit mpmath references. They cover x in [3.5, 30] and negative x at 64/128/256 bits, `erfc` at 1024 bits, and the `normal_cdf` lower tail. A second test checks erf odd symmetry and erfc reflection. On 1.0.4, 109 of the 221 cases fail, and so does the symmetry test.
- New integration test `tests/edge_values.rs`: `elliptic_f(±1 | 1)` = ±∞.
- `no_std` harness `nostd-tests/` (77 cases) and `scripts/ci_nostd.sh` (see above).
- New integration test `tests/c2f1_log_mpmath.rs`: 54 complex \({}_2F_1\) cases near \(z=1\) with integer \(c-a-b\in\{-3,\ldots,3\}\) and near-integer \(c-a-b\), at 128 and 256 bits, each to \(p-4\) bits. The first 1.0.5 candidate returned NaN for 34 of them.
- New integration test `tests/complex_bessel_mpmath.rs`: 218 cases (J/Y with \(\operatorname{Re} z<0\) in the Hankel regime, I/K on the positive real axis with real non-integer \(\nu\), K on its cut), at 128 and 256 bits, each to \(p-4\) bits norm-wise.
- 32-bit targets (i686): two tests assumed 64-bit limbs; the library was correct. `ext::test_ext` asserted that a 128-bit value is stored inline, which holds only with 64-bit limbs (4 limbs > `INLINE_WORDS` = 2 on 32-bit); the assertion is now per limb width. The `add_commutes` property test built values at \(p=32\) (one 32-bit word), below the documented \(p\ge64\) minimum of `from_i64`, and got `NaN(InvalidArgument)`; it now uses \(p\) in multiples of 64 bits. `cargo test -p zenith-float-num --target i686-unknown-linux-musl` now passes.
- `cargo fmt --all` applied (1.0.4 was not `fmt --check` clean). Two clippy lints in the MPFR differential tests (`--all-features`) were fixed.

### Known issues (fix planned)

- **Complex `bessel_k` is slow for non-integer \(\nu\) at high precision when \(\lvert z\rvert\) is just inside the series regime** (\(\lvert z\rvert<0.35(p+112)\)). Still open; queued for 1.0.7.
- **Real `bessel_j_nu` / `bessel_y` return `NaN(InvalidArgument)` for \(\lvert x\rvert\gtrsim10^4\)**: closed in 1.0.6 (real Hankel expansion).
- **Complex Bessel on the negative real axis** (\(\operatorname{Im} z=0\), \(\operatorname{Re} z<0\)) with real \(\nu\): still open; queued for 1.0.7.

## 1.0.4 — 2026-09-20

Clippy debt clear under `cargo clippy --workspace --all-targets -- -D warnings`. SoftFloat IEEE inherent ops renamed to `soft_*` where they collided with `std::ops`; needless range loops / doc list / `?` tidy; workspace `clippy.toml` thresholds for SoftFloat kernels.

## 1.0.3 — 2026-09-19

Coordinated patch with latex-rust, hdf5-rust, and redb-view (pure-Rust FOSS family adjacent to Accumath; Accumath itself stays proprietary).

## 1.0.2 — 2026-09-04

This crate is a numeric library. CSV and the 16-byte binary record remain.

- Removed HDF5 I/O and the `hdf5` Cargo feature. This crate does not read or write that file format.
- Build sheets, checklists, and walk lists are gone. Capabilities, contributing, getting started, help, library, and the other user guides stay.

## 1.0.1 — 2026-09-03

Jacobi elliptic functions on the real line.

- `ExactNum` Jacobi family: `am`, `sn`, `cn`, `dn`, and the nine quotients `cd` `ns` `nc` `nd` `sc` `sd` `cs` `ds` `dc`. Parameter \(m=k^2\in[0,1]\). \(m=0\) is trigonometric; \(m=1\) is hyperbolic; otherwise AGM / descending Landen. Cap `JACOBI_AGM_MAX = 128`. Domain errors are `NaN`.
- Period reduction by \(4K(m)\) runs only when \(\lvert u\rvert\ge\pi\) and uses a single Carlson \(K\) pass (no nested Ziv). `RoundingMode::None` at 512 bits is golded.
- `expr!` leaves for every Jacobi method above.
- Elementwise wrappers on `ExactNumArray` and software IEEE arrays for the full Jacobi set (including quotients).
- Hex limb gold `jacobi_sn_1_half` at 64 bits.
- Golds: zeros at \(u=0\); \(m=0\) vs sin/cos; \(m=1\) vs tanh/sech; Pythagorean identities; `sn(K/2)`; `sn(u+4K)=sn(u)`; \(\partial_u sn=cn\,dn\); `ns·sn=1`; `F(sn(u|m)|m)=u`.

## 1.0.0 — 2026-08-31

First crates.io release.

- HDF5 array I/O (feature `hdf5`): crates.io `hdf5-rust` 1.0, no `libhdf5`. `Ieee64Array` / `Ieee32Array` / `ExactNumArray` `to_hdf5` / `from_hdf5`. Golds: 100×3 binary64 bits; 50×50 `ExactNum` at 256 bits; 1-D 1000 binary32 bits; nested `results/data`; append 10×3→20×3; wrong dataset name `Err`. CSV unchanged.
- Pre-publish gate: `scripts/zenith_prepublish.sh` (12 checks) and `scripts/ci_full.sh`. dashu 0.6.0 has no `euler_gamma` and no scoped rounding closure. `LIBRARY.md` `expr!` leaves include `gammainc_upper`, `ai`, `bi`.
- Hex limb CI: `golds/hex/reference.txt` (35 rows) bit-identical on i686 musl, wasm32-wasip1, and aarch64 musl vs x86_64 `to_bytes()`.
- Integer SIMD for software IEEE arrays: add/sub/mul/div/sqrt/fma (`IEEE_SIMD_LANE_WIDTH=4`). Gold: 1000-element `Ieee64Array` bit-identical to the scalar kernel and to `ExactNum` rounded to binary64.
- `lazy_static` uses `spin_no_std` so `thumbv7em-none-eabihf` builds with `default-features = false`.
- Locked `expr!`/`cexpr!` composite golds: `erf+erfc=1` at 256 bits; `J_0²+Y_0²` at working precision; complex `erf` at `p_wrk`.
- Criterion specials (64–1024 bits) and linalg (matmul/FFT/LU). `compare-bench.sh --quick` includes dashu. Baselines refreshed.
- `ExactComplex::mul` rounds each of the four real products at `(p, rm)` before the add/sub, matching the MPFR componentwise gold. `RoundingMode::None` series paths are unchanged.
- `proptest` (`PROPTEST_CASES=1000`): add commutes; directed round-then-coarser; `erf` odd; small-integer 2×2 matmul associativity; `ExactRational` `(a+b)-b=a`.
- `mpfr-tests`: GMP oracle for `ExactRational`; real-axis complex `erf`/`Γ`/`Ai`/`J_n` vs MPFR; `mpfr_gamma_inc`; identity golds where GNU MPFR/MPC have no function. GNU MPC 1.3 has no `erf`/`gamma`/Bessel and is not used.
- `HELP.md` rewritten as the user guide: precision model, rounding, `Consts`, macros vs methods, 30 recipes, mistakes, FAQ.
- `GETTING_STARTED.md` covers software IEEE, arrays, `(p, rm, cc)` specials, `ExactRational`/`ExactInt`, `cexpr!` cuts, and `Ball`.
- CSV for `Ieee64Array` / `ExactNumArray`: bit-pattern cells (empty → `NAN`); `Display@p=` for exact arrays. Golds: 100×3 bit round-trip; missing cell is `NAN`; extra column `Err`.
- Binary interchange: 16-byte big-endian inline record (`to_inline_bytes` / `write_inline_bytes`); heap `u32` limbs for wider mantissas; `ExactNumArray` packed records. Golds: finite / Inf / NaN flag / `p=256` heap / array shape / invalid `Err`; `u32::MAX+1=2^{32}` at `p=64`.
- Serde (`serde` feature): `ExactNum` / `ExactComplex` / `ExactRational` / `ExactInt` / arrays / `Ball`. Decimal strings carry `@p=`. IEEE arrays serialize integer bit patterns. Shape mismatch is `Err`.
- `IEEE_SIMD_LANE_WIDTH=4` public constant. `From<LayoutError>` → `MemoryAllocation`. LU/QR/SVD workspace `try_reserve_exact` returns `None` instead of aborting on reserve failure.
- SHA-256 / SHA-512 / HMAC-SHA-256 / `constant_time_eq`. Golds: FIPS empty and `abc`; RFC 4231 HMAC TC1.
- Modular `ExactInt`: `mod_pow`, `mod_inv`, `miller_rabin`, `pollard_rho` (Brent). Golds: `2^{100} ≡ 976371285 (mod 10^9+7)`; `mod_inv(3,7)=5`; Miller–Rabin on `2^{31}−1`; `8051=83×97`.
- Window functions: symmetric Hann / Hamming / Blackman, Kaiser (`I_0`), rectangular. Golds: `hann(4)=[0,3/4,3/4,0]`; Hamming endpoints `0.08`; Kaiser `β=0` is rectangular.
- Discrete transforms: `dct`/`idct` (type II / III), `dst`/`idst`, `fft_real`/`ifft_real`. Golds: `idct(dct(x))=x`; constant DCT is DC only; cosine energy at bins `k` and `N-k`; Parseval. `DSP_MAX_POINTS=2048`.
- ODE solvers: `rk4`, `rk45_adaptive` (Dormand–Prince 5(4)), `euler`. Golds: `y'=-y` RK4 1000-step error `<10^{-12}`; RK45 meets `atol=10^{-12}`; Euler 1000 vs 2000 is `O(h)`.
- Root finding: `bisect`, `newton`, `brent`, `illinois`. Golds: `bisect(sin,[3,4])=π`; `newton(x²−2)=√2`; Brent fewer iterations than bisection; `bisect(sin,[0,1])=None`.
- Quadrature: `gauss_legendre`, `tanh_sinh`, `gauss_laguerre`, `gauss_hermite`. Golds: \(x^2\) on \([-1,1]\) is \(2/3\); 20-point \(x^{38}\) is \(2/39\); \(1/\sqrt{1-x^2}=\pi\); Laguerre \(x^2=\Gamma(3)=2\); Hermite \(1=\sqrt{\pi}\).
- Orthogonal polynomials: `hermite_he`/`hermite_h`, `laguerre`/`gen_laguerre`, `chebyshev_t`/`chebyshev_u`, `gegenbauer`. Golds: `He_4(0)=3`, `L_3(0)=1`, `T_5(cos(π/5))=-1`, `C_2^{(1)}=4x²−1=U_2`, `2C_2^{(1/2)}=3x²−1`, `T_6=2xT_5−T_4`. `ORTHOPOLY_N_MAX=256`.
- Chebyshev interpolation: `chebyshev_coeffs` / `chebyshev_eval` / `clenshaw` / `chebyshev_error_bound`. `exp` on `[-1,1]` with 20 terms error `< 10^{-15}`; Clenshaw of `[1,2,3]` at `1/2` is `1/2`; `CHEBYSHEV_MAX_DEGREE=256`.
- `ExactNumPoly` — dense univariate (low-to-high coeffs); `div_rem(x²−1, x−1)=(x+1, 0)`; `gcd=x−1`; `compose(x², x+1)=x²+2x+1`; `∂(x³)=3x²`; `roots_real(x²−2)=±√2` (companion closed form through `POLY_COMPANION_CLOSED_DEG=2`).
- RNG: `random_uniform`, `random_gaussian` (Box–Muller), `random_exponential`, `ExactNumArray::random_fill` / `RandomDist`. Existing `random_normal(p, exp_from, exp_to)` is unchanged.
- Distribution kernels: `normal_pdf`/`cdf`, `gamma_pdf`, `beta_pdf`, `poisson_pmf`, `binomial_pmf`, `chi_squared_cdf`, `student_t_pdf`. Golds: \(1/\sqrt{2\pi}\), \(1/2\), \(e^{-1}\), \(\chi^2_2(2\ln 20)=19/20\).
- `ziv_round_vec` — same Ziv loop on a slice; `hypot(3,4)=5`; `atan2(1,1)=π/4` at 256 bits.
- `parse_exact` / `format_exact` — `0.1` is exact `1/10`; dyadic `0.5` / `0.125`; `parse("1.5e3")=1500`.
- `ExactInt` — limb integer; `20!`, `gcd(48,18)=6`, `2^100`, `div_rem(17,5)=(3,2)`.
- `ExactRational` — exact `num/den` with integer GCD reduction; `1/3+1/6=1/2`; `2/4=1/2`; sign onto the numerator; 256-bit `1/3`.
- `doc/REPRODUCIBILITY.md` — what determines a result, how to replay tests, citation line.
- `# Precision` rustdoc on `ExactNum` and `ExactComplex` specials (algorithm, thresholds, Ziv / ULP, MPFR oracle).
- `ExactNumArray::fft` / `ifft` — radix-2 Cooley–Tukey; impulse / cosine / Parseval golds; `FFT_MAX_POINTS=4096`.
- `ExactNumArray::eigen_decomp` — symmetric QR; \(Av=\lambda v\), \(V\Lambda V^T=A\); \(\begin{pmatrix}2&1\\1&2\end{pmatrix}\to(3,1)\); `EIGEN_ITER_MAX=64`.
- `ExactNumArray::svd_decomp` — Golub–Reinsch; \(U\Sigma V^T=A\), \(U^\top U=V^\top V=I\); \(\operatorname{diag}(3,2)\to\sigma=(3,2)\); `SVD_ITER_MAX=64`.
- Living docs are the build plan, capabilities, and build checklist only. `ZENITH_FLOAT_TODO.md` is retired.
- `ExactNumArray::qr_decomp` — modified Gram–Schmidt; \(QR=A\), \(Q^\top Q=I\); rank-deficient zero diagonal.
- `ExactNumArray::lu_decomp` — partial pivoting; \(PA=LU\) gold; singular → `None`.
- `ComplexBall` disk enclosures: `add`/`mul`/`exp`/`ln`/`sin`/`cos`. Unit-disk `exp` and \(\sin^2+\cos^2=1\) golds.
- `Ball::{cos,ln,sqrt,erf,bessel_j0,bessel_j1}` certified enclosures (`BALL_TRANSCENDENTAL_ERROR_TERMS`). `ExactNumArray` `(2×3)` `sin` and `bessel_j_nu` golds; `signum` ufunc.
- `ExactComplex::hypergeom_2f1` — series / Euler / Pfaff / Kummer in \(\mathbb{C}\). `cexpr!` leaf `hypergeom_2f1`. Cut on \([1,+\infty)\); non-positive integer \(c\) → NaN.
- `ExactComplex` elliptic \(K,E,\Pi\) (complete and incomplete) via Carlson \(R_F,R_C,R_D,R_J\) in \(\mathbb{C}\). `cexpr!` leaves `elliptic_k` / `elliptic_e` / `elliptic_f` / `elliptic_e_inc` / `elliptic_pi` / `elliptic_pi_inc`. \(K(1)=+\infty\); cut of \(K\) on \([1,+\infty)\).
- `cexpr!` — per-part `errs[]` for real vs imaginary cancellation; principal branch cuts documented; leaves match the complex-capable `expr!` set (`cbrt`/`root`, logs/exps, `hypot`/`fma`, `abs`/`arg`/`conj`, `ldexp`/`scalb`/`logb`). No `atan2` or `rem_pi` (complex trig uses `x+iy` identities).
- `ExactComplex::{log2, log10, log, log1p, exp2, exp10, expm1, ldexp, scalb, logb, cbrt, nth_root, hypot, fma, mul_add}`.
- Getting-started guide (`doc/GETTING_STARTED.md`) and help tutorial (`doc/HELP.md`).
- Maintainer build checklist is not packaged on crates.io.
- `Context::with_rounding_mode`.
- `ExactNum::{two_sum, two_product, fused_sum, fused_dot, polyval}`.
- `ExactComplex`: `tan`, `sinh`, `cosh`, `tanh`, `sqrt`, `pow`, `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`.
- `ExactComplex` rectangular complex type (`+ − × ÷`, `exp`, `ln`, `sin`, `cos`, `abs`, `arg`).
- `ExactNum::sin_cos`; `mul_add` alias of `fma`.
- Special functions: `erf` / `erfc`, `gamma` / `ln_gamma`, integer-order `bessel_j`; `expr!` leaves.
- MPFR compare coverage for `fma`, paired `sin_cos` / `sinh_cosh`.
- MPFR compare coverage for `exp2`, `exp10`, and `rem_pi`; release gate on Linux x86_64 in `scripts/ci.sh`.
- MPFR fuzz: `add` / `mul` / `sqrt` bit-exact under all five IEEE rounding modes (`tests/mpfr/fuzz_round_modes.rs`).
- Progressive `ConstCache` (`cache_info`, √2, φ) and `CachedFBig` extra-precision wrapper.
- Public `Ball` enclosures and `ziv_round` correct-rounding retry helper.
- Stack-inlined mantissas (`INLINE_WORDS` = 2 limbs) before heap promotion.
- `doc/EXPR.md` — `expr!` per-op rounding contract (working precision, cancellation slots, final `set_precision`).
- `fma` skips the full-width product when `|a*b|` and `c` differ by more than `p + 2` words.
- MPFR fuzz covers `nth_root` for n=4,5 at 1-ULP (`mpfr_rootn_ui`).
- Optional nightly bench gate: `CI_BENCH=1 ./scripts/ci.sh`.
- `Consts::euler_gamma` (Euler–Mascheroni γ) on the progressive constant cache.
- `ExactNum::{frexp, ldexp, scalb, logb, ilogb}` — IEEE-style exponent split without hardware floats.
- `LowerExp` / `UpperExp` (`{:e}` / `{:E}`) and `LowerHex` formatting.
- `SharedConsts` — mutex-wrapped constant cache for parallel batch evaluation (`std`).
- MPFR oracles for `erf`/`erfc`, `gamma`/`ln_gamma`, `bessel_j`, extra constants, and `ExactComplex` add.
- `expr!` constants `sqrt2`, `phi`, `euler_gamma`.
- `expr!` leaves `ldexp` / `scalb` / `logb`.
- `serde`: JSON round-trip for `ExactNum` (decimal string / integer, Inf/NaN) and `ExactComplex`; feature implies `std`.
- Structured radix 2–36 parse/format round-trip stress test (`tests/radix_roundtrip.rs`).
- Native `+`, `-`, `*`, `/` operators on `ExactNum` (default precision, round-to-even).
- `exact!` / `fbig!` compile-time decimal float literals.
- Parse/format for radices 2–36 (`Radix::try_new`); `_e` exponent separator for bases > 10.
- `RadixFloat` wrapper and `ExactNum::with_radix` for radix-tagged I/O.
- `ExactNum::nth_root` / `expr!` `root(x, n)` — general n-th root (`sqrt`/`cbrt` delegation, composite factors, Newton for primes).
- `ExactNum::sinh_cosh` — paired hyperbolic evaluation with a single `exp(|x|)` path.
- `ExactNum::copysign`, `ExactNum::next_after` — software IEEE sign and successor.
- MPFR bit-oracle for `copysign` at high precision (`compare_copysign_test`); sign is applied before rounding so directed modes match MPFR.
- Seeded test RNG (`reseed_random`, `ZENITH_TEST_SEED`; panic hook prints the seed).
- `scripts/ci.sh` enforces debug / MPFR wall-time budgets (10 min / 30 min).
- Allocation failure at huge precision returns `NaN` (`Error::MemoryAllocation`), with a large-precision add/mul smoke test.
- `ExactNum::fma` / `expr!` `fma(a, b, c)` — fused multiply-add with single final rounding.
- Criterion benchmarks in `zenith-float-num/benches/` (arithmetic, transcendentals, composite); `./scripts/bench.sh`.
- FFT-scale multiply benchmark tier (`arithmetic/mul_fft`, 346k and 524k bits).
- TSV benchmark baselines (`doc/bench-baselines.tsv`) and `scripts/bench-compare.sh` for regression checks without JSON.
- Cross-library compare harness (`zenith-float-compare/`, `scripts/compare-bench.sh`) for release vs astro/dashu.
- Newton reciprocal at three words and up.
- `hypot`, `atan2`, `log1p`, and `expm1` on `ExactNum` and in `expr!`.
- `exp2`, `exp10`, and public `rem_pi` on `ExactNum` and in `expr!`.

## 0.1.0

First public release.

- Software-limb `ExactNum` (no hardware floating-point in calculations).
- Elementary functions, constants cache, radix parse/format, `expr!`.
- Optional `std`, `random`, and `serde` features.
- Dual license MIT OR Apache-2.0.
