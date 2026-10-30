# zenith-float audit — 2026-10-30 (read-only)

Source audited: `git archive origin/main` @ dc6bc87 (workspace 2.0.0), plus the published crates.io tarballs zenith-float 1.0.12 / zenith-float-num 1.0.12 (published 2026-10-06 11:03 CDT; no git tag for it). All paths below are relative to `zenith-float-num/src/` on main unless prefixed. Nothing in ~/MathWorks was modified; no cargo build or test was run. **Coverage figures are grep heuristics, not test runs.**

## Blunt verdict on "a $2 calculator"
That claim is false. The reachable public surface is 663 items (inventory.csv). The crate has 523 real implementations: elementary, trig and hyperbolic functions with inverses, gamma, digamma, erf, Bessel, Airy, zeta, polylog, elliptic, Jacobi, orthogonal polynomials, complex, rational, interval and dist. There are 0 `todo!`/`unimplemented!` stubs. **The real problem is correctness, not breadth.** Published 1.0.12 ships dozens of wrong-answer bugs that were fixed only on unpublished main. Many functions have no proven accuracy, and the release checklist's MPFR sweeps are still open.

## Counts
- inventory.csv: 1159 rows. 663 reachable from the crate root, 392 internal, 92 macro-impl, 9 proc-macro.
  - Status: real 523 (+51 trait impls), delegating wrapper 411 (+131 trait, +32 macro-generated), trait decl 11, stub 0, always-error 0.
  - Documented accuracy (reachable): correctly-rounded 177, unspecified 29, working-precision-32-bit 10, enclosure 2, exact 1, **no statement 443**.
- Published 1.0.12: 627 reachable items, **zero** accuracy docs. 45 reachable items exist only on main (operator impls Add/Sub/Mul/Div, div_exact, pow2, sqr, beta, elliptic_*_phi, …).
- mpfr_gap.csv: 349 MPFR 4.x API names (source: the official MPFR manual function index, https://www.mpfr.org/mpfr-current/mpfr.html, fetched 2026-10-30). Present 192, MISSING 139, n/a (C memory/build API) 18. Most of the 139 are flags, printf, fits_* and C type I/O. The **math** functions that are missing are agm, const_catalan, rec_sqrt, rootn_si, powr, compound_si, dim, fmma/fmms, sec/csc/cot, sech/csch/coth, sinpi/cospi/tanpi, asinpi/acospi/atanpi/atan2pi, the *u family, exp2m1, exp10m1, log2p1, log10p1, fac_ui, remainder, remquo, fmodquo, modf, round-half-away, get_d/set_d on ExactNum, and get_q (float→rational).
- edge_cases.csv: 142 ExactNum math functions × 17 classes = 2414 rows. Heuristic result: covered 1074, **UNCOVERED 980**, n/a 360.
  - `jacobi_arccn` has no test hit at all.
  - Worst covered: lambert_wm1 (1/16), hypergeom_0f1, jacobi_arcdn, chebyshev_t/u, gegenbauer, gen_laguerre, hypergeom_pfq, struve_h, binomial_pmf, assoc_legendre_p.
  - Even core functions miss classes. acos misses 8, including branch boundary, all rounding modes and ties. ceil misses 8, including negative, huge/tiny exponent and exact. asin, ln_gamma: 7. exp10 misses negative and all rounding modes.
  - The MPFR random harness (tests/mpfr, random p 1..157 plus one p=1563 run) is behind `--features mpfr-tests` (Cargo.toml:59-60, tests/integration_tests.rs:6). Per RELEASE-2.0-CHECKLIST.md:133 it has not been run on the release SHA.
- 108 commits on main since 1.0.12 mention fix, wrong or known_wrong (fixes_after_1.0.12.txt).

## Top-20 gaps, ranked by risk of a wrong answer reaching users of the published crate
1. **Published 1.0.12 returns wrong special values from pow**: pow(±0,−Inf), x^±Inf, (−Inf)^negative-odd, NaN^0 and 1^NaN. Fixed only on main in 9923729, c953ad8, f5c40b5, c5dabeb and ad462f0. 1.0.13 was cancelled (RELEASE-2.0-CHECKLIST.md:22), so users still get the wrong values.
2. **1.0.12 riemann_zeta is wrong for real s**: zeta(1) is not invalid, and zeta(−2k) is not +0 (fixed 87789cb). On main it is still "Accuracy: unspecified" (catalog.rs:1116).
3. **1.0.12 gamma/ln_gamma/digamma**: gamma(−0) sign (507c685); ln_gamma poles and negative x (36e388e); digamma (26af74a, 7307f44); gamma timeouts (9a817e2).
4. **1.0.12 core arithmetic edge bugs**: full-precision add/sub/mul subnormal (0cfee92); finite/±Inf division sign (87ae1fb); fma exponent wrap (c45708b); fma far addend and subnormal reciprocal (RELEASE-2.0-CHECKLIST.md:20); recip (2b931d2).
5. **1.0.12 elementary functions**: expm1 (4acddc5), log1p overflow (2f466af), log10 on exact powers (45a6243), log2(1) (1a4327e), log with base (468ad72), pow underflow (da74a79, d783ed3, f78f4d6, 79ca61c), acosh guard bits (cab00e8), atan2 huge/tiny returned Inf (0099044), ldexp (825e93e), nth_root (a068c0c, de2ded0), negative nth_root with odd n≥5 (checklist:20).
6. **1.0.12 constants**: pi, e, ln2 and ln10 were not rounded at exactly p (c04ee8f). Every function that uses them inherits this.
7. **1.0.12 erf/erfc**: directed rounding modes were wrong (a219bee), and erfc underflow returned −0 (50a81d2). On main, erfc still has no detected test for overflow, ties or high precision (edge_cases.csv).
8. **1.0.12 complex mul/div** were wrong (0cdb297, bc9df25, f97212a, de02bb2).
9. **round_exact makes unproven functions look correctly rounded.** It rounds f(round_p(p)+64, ToZero) to odd and then to p (ext.rs, round_exact just before round_odd_path). That is only correct if f is accurate to p+64 bits. The 29 catalog functions marked "Accuracy: unspecified" (catalog.rs:575-1581: hurwitz_zeta :1067, riemann_zeta :1116, lambert_w0 :1470, polylog :1582, polygamma …) and 4 in heun_fox.rs (:576, :640, :771, :967) carry no such guarantee. The output therefore looks exact when it may not be. 1.0.12 also had a bug where round_exact's fallback returned an unrounded value (388afc6).
10. **dist.rs functions can silently return wrong digits.** gamma_pdf :156, beta_pdf :193, poisson_pmf :239, binomial_pmf :266, chi_squared_cdf :304 and student_t_pdf :330 compute at p+64 with RoundingMode::None and have no Ziv check. They are self-described as "MPFR oracle: no". Cancellation (for example in the cdf tails) gives garbage that is presented as a value.
11. **Operators silently compute at 128 bits.** `a+b`, `a-b`, `a*b` and `a/b` (ext.rs:4610-4615 and siblings) use DEFAULT_P=128 (defs.rs:113) with ToEven, whatever the operand precision. A user with 1024-bit operands gets a 128-bit result and no error. From<i8..u128> (ext.rs:4529-4548) and Default (ext.rs:4493) also use 128.
12. **The MPFR sweeps the release relies on were not done**: full suite on the release SHA (RELEASE-2.0-CHECKLIST.md:133); exponent extremes across all 5 modes (:137); high-precision comparison (:138); retry budget at p≥8192 (:139). The random harness draws exponents only from the i32 band (tests/mpfr/common.rs:226-231), so 62-bit exponent extremes are not randomly tested.
13. **Ziv retry exhaustion.** MAX_PREC_RETRY=256 (common/util.rs:80) returns Error::PrecisionRetryExhausted (defs.rs:165). That is acceptable as a "no answer". However, it has not been validated at p≥8192 (checklist:139), and how hard-to-round cases near the budget behave is unknown.
14. **Trig at huge arguments.** sin (ops/sin.rs:52) reduces modulo 2π at mantissa+extra bits (ops/util.rs reduce_trig_arg). There is no detected test of cos at a huge exponent (edge_cases.csv: cos "huge exponent" UNCOVERED), and the exponent-extremes sweep is open. rem_pi (ops/util.rs:13) ignores its rounding-mode argument.
15. **sin_cos and sinh_cosh do not go through round_exact.** sin_cos (ext.rs:2819) just calls sin and cos, and sinh_cosh is at ext.rs:2887. If their component paths differ from sin and cos, the results can diverge. They have no independent accuracy statement.
16. **443 reachable items have no accuracy statement at all**, and 10 are documented as only "working precision, 32 bits" (inventory.csv documented_accuracy). Users cannot tell which functions are guaranteed.
17. **Missing functions that users will build by composition and get wrong**: sinpi/cospi (sin(π·x) loses exact zeros at integers), rec_sqrt (1/sqrt double-rounds), sec/csc/cot, exp2m1/log2p1 and relatives (catastrophic cancellation near 0), IEEE remainder/remquo (rem is fmod only, ext.rs:723), agm (only a private agm_small, ops/jacobi.rs:305), and const_catalan.
18. **ln_gamma differs from MPFR**: it returns ln|Γ| with no sign output, so the sign of Γ is lost (mpfr_gap.csv: lngamma/lgamma). Code ported from MPFR's lngamma (NaN for Γ<0) will get a silent value instead.
19. **Integer-argument limits**: bessel_j takes n: usize (ext.rs:3314), so jn(−n) is impossible. nth_root takes usize (ext.rs:2679), so there is no rootn_si. powi takes usize and powsi takes i64, with no big-integer exponent.
20. **Weakly tested families.** Orthogonal polynomials, hypergeometric functions, Jacobi inverses, struve_h and lambert_wm1 each cover only 1–2 of their applicable edge classes (edge_cases.csv). jacobi_arccn has no test reference at all.

## Places that return a value for input they should refuse
- **FromStr** (ext.rs:4500-4507): malformed text returns `Ok(NaN(InvalidArgument))` instead of `Err`. `"abc".parse::<ExactNum>()` succeeds.
- **Eq** (ext.rs:4455): `impl Eq for ExactNum {}`, while PartialEq (ext.rs:4441) makes NaN≠NaN. This breaks Eq's reflexivity contract, so HashMap/HashSet/dedup behave incorrectly for NaN.
- **Operators** (ext.rs:4610-4615): they return a 128-bit-rounded result for higher-precision operands instead of refusing or using the operand precision (see #11).
- **dist.rs** (:156-:330): returns values with no error bound under cancellation (see #10).
- **"Unspecified" catalog functions** (catalog.rs:575-1581, heun_fox.rs:576/640/771/967): they return correctly-rounded-looking values with no accuracy proof (see #9).
- **ln(±0)** = −Inf (gen_wrapper_log, ext.rs:2510) matches MPFR/IEEE (divide-by-zero), so it is acceptable. However, there is no flag API to tell the user that a pole was hit (mpfr_gap.csv: divby0_p MISSING).
- **Documentation contradiction**: root lib.rs says "Precision is always rounded up to the nearest word", but round_exact rounds to the exact p. One of the two is wrong for users.
- Fixed on main but **live in 1.0.12**: p=0 was not rejected (ceb15e7, 6c85925, d20f094), min_positive(0) panicked (7e1b7d3), and zeta(1) was not invalid (87789cb).

## Method / caveats
- inventory.csv comes from extract.py (static parse; strips cfg(test); expands gen_wrapper macros). The status heuristic is: ≤3-line body = wrapper.
- edge_cases.csv uses a ±12-line keyword window around each call site in the test files. A "covered" result means evidence exists, **not proof**. "UNCOVERED" means no keyword evidence was found. Functions exercised only by the feature-gated MPFR random harness may show precision-1 or rounding-mode classes as UNCOVERED, even though that harness draws p from 1..157.
