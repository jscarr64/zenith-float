# zenith-bound plan review: suggested additions

Review of `zenith-bound_PLAN_2026-10-14.md` (revision 2026-10-14). Prepared 2026-10-14.

The design is sound. The additions below fall into five groups: inconsistencies in the plan, missing rules that affect soundness, missing scope, testing, and open-source readiness.

## 1. Inconsistencies in the plan

- **Part G total.** The 10-14 column sums to 175–**315**, but the headline says 175–335. Either a row is off by 20 or the headline is stale. The "5 → 11 items" for zf 2.1 also doesn't match the 13 items (P1–P13) listed in C1.
- **Directed rounding of transcendentals is unaudited.** K4 only audited add, sub, mul, div and sqrt. Part H says to make Up/Down MPFR coverage "for every kernel zenith-bound uses" a C1 gate, but P1–P13 has no such item. Add a **P14**. Directed-mode Ziv loops also have a known trap: exact cases like `exp(0)` never terminate unless special-cased. The alternative is to use `ToEven` with a proven ≤1 ulp bound.
- **Three exponent conventions.** zf uses 0.m×2^e, Part I says value = m·2^e, and `Mag` uses mant·2^(exp−63). The `zf-ball-enc` layout never says which one `exp i64` means. Pin one in `CONTRACT.md`.
- **The encoding has a NaN tag.** `tag 3 nan` contradicts "never store NaN". The decoder should reject it.
- **Negative zero.** K2 says exact cancellation under `Down` gives −0. Normalize −0 to +0 at the ball boundary, or byte-equality breaks.
- **BLAKE3.** It needs the `pure` feature. Otherwise it pulls in SIMD, asm and CPU detection, which conflicts with Parts B and J. The `op_ver` of "2.2.0" in the Part I example also doesn't match the 0.x versioning.

## 2. Missing rules that affect soundness

1. **Zero `Mag`.** Zero is `Finite{0,0}`, but `add_up` and `mul_up` align by exponent, so zero would be treated as 2^0. Every `Mag` function needs a zero early-return, and zero belongs in every test set. Also decide `0·∞` (recommended: `Infinity`).
2. **`Mag` underflow.** Only exponent overflow is specified. Underflow below `i64::MIN` must round **up** to the smallest positive `Mag`, never to zero and never by wrapping.
3. **`Mag` exponents can fall outside zf's range.** C6 itself notes `exp(mid) − prec − 1` can go below `i32::MIN`, so `to_exact_num` can fail for valid radii. That means `contains`, `is_unique_at` and the sqrt domain check can't go through it. They need dedicated functions such as `Mag::cmp_exact_num` that never materialize the radius.
4. **Exact add can blow up memory.** `add_full_prec` on operands with a huge exponent gap (up to \~2^32 bits) needs a mantissa of hundreds of MB. Only do exact work when the gap is at most `prec + guard`. Otherwise use the sticky/directed path or return `Resource`. The sqrt "exact subtraction" is unnecessary anyway, because `m − r > 0` is just `m > r`. Add this as an adversarial set.
5. **No division formula.** The contract gives add and mul rules only. With `q` the midpoint quotient and `|q|` bounded from above: `rad = (r_a + |q|·r_b) / (|m_b| − r_b) + ε`, with the denominator rounded **Down**. If the ball touches 0, return `Domain`.
6. **exp radius.** `exp(r) − 1` cancels badly for small `r`. Use `r·exp_up(r)`, which is always valid because `e^r − 1 ≤ r·e^r`.
7. **One domain policy.** sqrt clips to the in-domain part, while log and division return errors. A single `DomainPolicy { Strict, Clip }` applied uniformly would be cleaner.

## 3. Missing scope

- **L0 ops.** Add `abs`, `sqr` (a ball straddling 0 should give `[0, (|m|+r)²]`, not `m² ± …`), `powi`, `hypot`, `sum`/`dot`, and `floor`/`ceil`/`sgn` returning `Tri` or `Indeterminate`. `union` and `intersect` are listed in gap 3 but absent from the L0 file plan. Avoid `PartialEq`/`PartialOrd` on `Ball`; use `identical()` and `Tri` instead.
- **L1 functions.** Add asin, acos, atan2, sinh, cosh, tanh and their inverses, log1p, expm1, log2, log10 and exp2. All have simple derivative-based radius rules.
- **Decimal parse and print.** Part F depends on it but the file plan has no `decimal.rs`. "0.1" must become an enclosing ball, and printing must round outward in integer arithmetic. The printer rule also needs a fallback for a ball like `[0.9999, 1.0001]`, where no common digit string exists.

## 4. Testing

- **Circular oracle.** Tests that take "exact truth" from zf's own `*_full_prec` would share zf's bugs. Use a small independent dev-only rational/bigint model, or Python-generated truth.
- **Mutation testing.** Run `cargo-mutants` alongside the "shrink the radius by 1 ulp" check.
- **7B-model safeguard.** Have a stronger model or a human review every `mag/*` file after 2026-10-30. One off-by-one there is unsound, and the test sets may not catch it.

## 5. Open-source readiness

- **Clean-room rule.** The plan cites Arb routines (`acb_calc_integrate`, `acb_elliptic_rf`, `arb_printd`). Arb/FLINT and MPFR (via `rug`) are LGPL. Implement from the papers, never from their source, and say so in `CONTRIBUTING`. Note in `NOTICE` that the oracles and the golden data derived from them are dev-only.
- **Claims policy.** Put a per-function status table (Guaranteed / Heuristic / Refused) in the README, and only say "rigorous" for functions whose contract tests pass.
- **Project files.** Add `SECURITY.md` treating containment failures as security-grade, plus a CONTRIBUTING covering the no-float rule, red-first tests, the B3 ban and the oracle setup.
- **Checks and licensing hygiene.** Run `reuse lint` on all files, not only the puremp headers. Add `cargo-deny`, `cargo-semver-checks`, and a promise that the encoding is frozen. Check that the crates.io name is free.
- **Rust version and `unsafe`.** `rust-version = 1.99` forces that toolchain on every consumer, so decide on an MSRV policy. The `fesetround` test needs `unsafe`, so put it in a separate dev crate to keep "no unsafe" true for the library.


---

## Dispositions (recorded in docs/PLAN.md Part R, commit ee5dbd5)


### R1. Inconsistencies
| # | Item | Disposition | Where |
|---|---|---|---|
| 1.1 | Part G totals didn't match the headline or the column sums | **Accepted.** The 10-14 column summed to 173–331 while the table said 175–335. After the review's scope changes the column sums to **195–345**; the headline, Summary and table total all say that. Part J's own effort line (16–30, plus weeks) now matches Part G's 13–22, and the leftover week counts in Parts I and J are gone | Summary, Part G, I, J |
| 1.2 | "5 → 11" vs C1 listing P1–P13 | **Accepted.** Now "5 → 14 (P1–P14)" in C1 and Part G | C1, G |
| 1.3 | Add P14: Up/Down MPFR coverage for every transcendental kernel; Ziv exact-case trap; ToEven + proven ≤ 1 ulp alternative; move Part H's gate into the P list | **Accepted.** P14 added to C1 and gates L1; Part H now points to P14. The stray "Option A in C1" (no such option exists) is removed | C1, L1, H |
| 1.4 | One exponent convention for zf `0.m·2^e`, Part I `m·2^e` and Mag `mant·2^(exp−63)` | **Accepted (convention), changed (Mag).** CONTRACT.md clause 4a pins zf's `0.m·2^e`. Mag keeps its fields (`u64` mantissa with the top bit set, `i64` exponent), but `exp` now means `0.mant·2^exp`, so the mul_up exponent formulas and the half-ulp boundary set were redone. The encoding and traces use the same convention. **Jeff should confirm, because C6 was locked** | B, C6, L 4a, I, J, M |
| 1.5 | Encoding tag 3 = nan contradicts the no-stored-NaN rule | **Accepted.** Tags 2 and 3 are reserved and the decoder rejects them; the whole line is mid +0 with `rad = Infinity`. Trace `nan`/`inf` mids are gone too | L 1, I, J |
| 1.6 | −0 from exact cancellation under Down | **Accepted.** −0 is normalized to +0 at the ball boundary; the encoding rejects sign = 1 on zero | K2, L 1, J |
| 1.7 | BLAKE3 pure feature | **Accepted.** `default-features = false, features = ["pure"]`: no SIMD, no asm, no CPU detection | I |
| 1.8 | op_ver example should be 0.x | **Accepted.** `"0.1.0"`; op_ver is defined as the zenith-bound version | I, J |

### R2. Soundness
| # | Item | Disposition | Where |
|---|---|---|---|
| 2.1 | Zero Mag early-returns everywhere; zero in every test set; 0·∞ = Infinity | **Accepted** | C6, L 6–7, L0 tests |
| 2.2 | Underflow below `i64::MIN` rounds UP to the smallest positive Mag | **Accepted.** `checked_*` everywhere; the low side gives `mant = 2^63, exp = i64::MIN` | C6, L 7 |
| 2.3 | `Mag::cmp_exact_num` etc.; contains / is_unique_at / sqrt domain never build the radius in zf | **Accepted.** `cmp_exact_num` and `cmp_exact_diff` added in `src/mag/cmp.rs`; `is_unique_at` is restated as a comparison against the rounding-cell boundaries | C6, L0 file plan |
| 2.4 | Exact add only within the gap ≤ prec+guard, otherwise the sticky/directed path or Resource; sqrt domain is the comparison m > r; adversarial huge-gap set | **Accepted.** guard = 64 bits | C6, L0 tests |
| 2.5 | Division radius `(r_a + |q|_up·r_b)/(|m_b|−r_b)_down + ε`; a divisor touching 0 gives Domain | **Accepted.** Domain under both policies (the old "Indeterminate when whole line" option is gone) | C6, L 6–7 |
| 2.6 | exp radius `r·exp_up(r)` | **Accepted.** `rad' = exp_up(m)·r·exp_up(r) + ulp`; the same pattern covers exp2/expm1/sinh/cosh | L1 |
| 2.7 | One DomainPolicy {Strict, Clip}, applied uniformly | **Accepted.** Default Strict; replaces the sqrt "strict domain gets Indeterminate" wording | C5, C6, L 7, L0 `policy.rs`, L1 |

### R3. Scope
| # | Item | Disposition | Where |
|---|---|---|---|
| 3.1 | L0: abs, sqr (straddle → `[0,(|m|+r)^2]`), powi, hypot, sum/dot, floor/ceil/sgn → Tri or Indeterminate, union/intersect in the file plan; no PartialEq/PartialOrd, use identical() and Tri | **Accepted** | C5, L0 file plan |
| 3.2 | L1: asin, acos, atan2, sinh, cosh, tanh and their inverses, log1p, expm1, log2, log10, exp2 | **Accepted.** Each has a radius rule, domain policy and adversarial cases; all need P14 | L1, C1 P14 |
| 3.3 | `decimal.rs`: parse → enclosing ball; print rounds outward in integers, with a fallback when there's no common digit string | **Accepted.** In the L0 file plan; the Part F printer points to it; the fallback is `mid ± r` or the interval form | L0, F |
| 3.4 | Adjust the estimates (relative runs only) | **Accepted.** zf 2.1 14–24, scaffold/OSS 3–5, B 5–9, L0 22–34, L1 28–46, new truth-model/mutants row 4–8; total 195–345 | G |

### R4. Testing
| # | Item | Disposition | Where |
|---|---|---|---|
| 4.1 | No truth from zf `*_full_prec` (circular) | **Accepted.** A dev-only rational/bigint model crate, or Python `Fraction`/mpmath/python-flint goldens | E, L0 tests, C6 |
| 4.2 | cargo-mutants plus the shrink-radius-by-1-ulp check | **Accepted.** Shrinking is a `cfg(zb_shrink_radius)` switch; cargo-mutants runs nightly on mag/ball/encode | E |
| 4.3 | Every `mag/*` file gets an extra stronger-reviewer review (Paul or Cassie) | **Accepted** | L0 review rule, CONTRIBUTING.md |

### R5. Open-source readiness
| # | Item | Disposition | Where |
|---|---|---|---|
| 5.1 | Clean-room rule for LGPL Arb/FLINT and MPFR in CONTRIBUTING; NOTICE says oracles and derived goldens are dev-only | **Accepted.** CONTRIBUTING.md has the rule; a new top-level `NOTICE` and a README section say it; `NOTICE-puremp` is unchanged | CONTRIBUTING.md, NOTICE, README.md, H |
| 5.2 | README per-function status table (Guaranteed / Heuristic / Refused); "rigorous" only after the contract tests pass | **Accepted.** Every row starts as "Not yet (pre-release)" | README.md |
| 5.3 | SECURITY.md: containment failures are security-grade | **Accepted** | SECURITY.md |
| 5.4 | CONTRIBUTING: no-float rule, red-first, the single-exact-input-test ban, oracle setup | **Accepted** | CONTRIBUTING.md |
| 5.5 | Reuse lint, cargo-deny, cargo-semver-checks, encoding-freeze promise | **Accepted** | C3 |
| 5.6 | Check whether `zenith-bound` is free on crates.io | **Accepted.** On 2026-10-14 the API returned HTTP 404, "crate `zenith-bound` does not exist", so the name was free. Nothing was reserved (`publish = false`) | C3 |
| 5.7 | MSRV policy: the frozen Rust 1.99 pin | **Accepted** | C3, CONTRIBUTING.md |
| 5.8 | The unsafe fesetround test goes in a separate dev crate | **Accepted.** `fpu-check/` (`publish = false`); zenith-bound keeps `forbid(unsafe_code)` | B.4 |


**Update 2026-10-14:** Jeff approved item 1.4. The Mag exponent convention `0.mant·2^exp`, with the redone `mul_up` and half-ulp boundaries from `ee5dbd59`, is LOCKED (docs/CONTRACT.md clause 1).
