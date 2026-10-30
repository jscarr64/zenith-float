# Proposal: split Stage B so L1 can start on a closed set of zenith-float functions

**Note (2026-10-16 8:28 AM CT, Jeff via Paul):** full scope, no scope cuts, no corner cutting; the ship date and the scope are Jeff's concern, not a team constraint. This split is an ordering for parallelism and quality only. It contains no schedule projection and no date comparison; Rusty will provide measured numbers after the Ziv-loop conversion commits land.

Status: **PROPOSAL, with Rusty's answers folded in as provisional recommendations.** Everything labelled "Rusty's recommendation 2026-10-16, Paul to confirm" is provisional and not frozen. Written 2026-10-16 (CT) on branch `spence/zb-stageb-docs`, docs only. It is lever 4 of `docs/REPLAN-2026-10-16.md` section 6, and follows the staged pin Paul ruled (DECIDED 2026-10-16, Paul under Jeff's delegation; PLAN C1 and the L0/L1 gates say so).

**Marks.** **[F: file:line]** read from a file or git. **[I]** my inference. **[R]** Rusty's recommendation 2026-10-16, Paul to confirm (provisional). zenith-float (zf) was read only, through `git show` on `/home/jsc64/MathWorks/zenith-float` (`wip/option-a` @ `fe5c15d`, local, 2026-10-16 07:56; `main` @ `abb2d5a` = 1.0.12); the checkout itself was not touched. zf line numbers are `zenith-float-num/src/ext.rs` at `fe5c15d`.

## 1. What is proposed

Stage B today (REPLAN section 1) is P11, P14 and Option A layer 1a for "the functions L1 calls". The proposal makes that list **explicit and closed** and orders it first:
- **B1** = the Ziv/transcendental family and the `round_exact` pivot for the functions in section 2, plus P11 and P14. **[R]**
- **B2** = Option A layers 2-8 and the final `Clone` removal. B2 must not change any B1 public item. **[R]**
L1 starts on B1. Every B1 function still gets its P14 finding and a green gate on real kernels; nothing is dropped.

## 2. The function list (derived from PLAN L1)

Source: PLAN:342 ("Functions: sqrt (from L0), exp, exp2, expm1, log, log2, log10, log1p, sin, cos, tan, asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, pow; constants pi, e, ln2, ln10") [F]; radius rules PLAN:343-349; P14 scope PLAN:195 (the same list and the constants; `sqrt` is not in it) [F]. A zf function is on the list if an L1 v1 kernel needs it for a midpoint or an enclosed constant. Radius arithmetic in `Mag` (`div_up`, `sub_lower`, Mag `sqrt`, `min(r, 2)`) is zenith-bound code and needs no zf function [F: PLAN:344-348]. **The list is agreed [R]: 22 functions plus the constants; `ln` is the natural log; `sin_cos`/`sinh_cosh` are not needed; zf `powi` is not used; `fma` is not needed.**

Signatures are from `wip/option-a` @ `fe5c15d` [F]. Single-argument functions are `(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum`. **Rusty's recommendation 2026-10-16, Paul to confirm: every function below keeps its current public signature; consuming `self` is not needed.**

| # | zf function (signature today = signature at B1 **[R]**) | Called by L1 kernel | Stage today | Effect of layer 1a / P14 / P11 |
|---|---|---|---|---|
| 1 | `sqrt(&self, p, rm)` (:2406; no `cc`, no retry loop) | `sqrt` | B | Implementation loses the `exact_operand` copy; signature unchanged **[R]**; not in P14 scope [F] |
| 2 | `exp` (`gen_wrapper_ovf_rm_cc`, :2244) | `exp` (mid and the `exp_up(m)` bound, PLAN:344) | B | Internals only **[R]**; P14 finding; accuracy contract changes only if the kernel is found not rigorous (`ToEven` + proven 1 ulp bound, PLAN:195, 349) |
| 3 | `exp2` (same macro) | `exp2`; exact case `2^n` | B | as #2 |
| 4 | `expm1` (:1347) | `expm1`; exact case `expm1(0)` | B | as #2 |
| 5-6 | `sinh`, `cosh` (`gen_wrapper_ovf_rm_cc`) | `sinh`, `cosh` | B | as #2 |
| 7-9 | `ln`, `log2`, `log10` (`gen_wrapper_log`, :2269) | `log` (= natural log, PLAN:346), `log2`, `log10`; exact cases `ln(1)`, `log2(2^n)` | B | as #2 |
| 10 | `log1p` (:1323) | `log1p` | B | as #2 |
| 11-15 | `sin`, `cos`, `tan`, `asin`, `acos` (`gen_wrapper_arg_rm_cc`, :2209) | same names; `sin(0)` exact case | B | as #2 |
| 16 | `atan` (:1209) | `atan`; also uses the constant π in the radius (`min(r, π)`) | B | as #2 |
| 17 | `atan2(&self, x: &Self, p, rm, cc)` (:1228) | `atan2` (quadrant logic and the (0,0) `Domain` are zenith-bound code) | B; W1 value fixes are Stage A [F] | as #2; second operand already borrowed |
| 18 | `tanh` (:1372) | `tanh` | B | as #2 |
| 19-21 | `asinh`, `acosh`, `atanh` (`gen_wrapper_arg_rm_cc`) | same names | B | as #2 |
| 22 | `pow(&self, n: &Self, p, rm, cc)` (:761) | `pow` with a non-integer or ball exponent; integer exponents are L0b `powi` built from Ball mul/sqr, not zf `powi` [I] | B | as #2; exact-case trap `pow(x, 0)` (PLAN:195); zf `pow(0,0)` is exactly 1, `inexact()` false on `wip/option-a` and 1.0.12 [F: `/tmp/zf-logs/zf-pow-zero-zero.md`] |
| 23-26 | `Consts::pi`, `e`, `ln_2`, `ln_10` (`&mut self, p, rm) -> ExactNum`, `ops/consts/mod.rs:213-240`) and the new `Consts::get_exact_prec` | constants; π also for `atan` radius, `tan` pole test, argument reduction; ln2/ln10 for exp-family reduction and the `1/ln2`, `1/ln10` scaling (PLAN:348, 350) | B (P11) | `get_exact_prec` is additive (P11) **[R]**. **Open: Jeff's decision 5** (borrowed `&ExactNum` second operand vs rebuild from retained integer series state) decides whether the existing getters' signatures change [F: `option-a-evidence.md` §3.6; `zf-pin-items-extract.md` P11 row]; Jeff has not ruled |
| (L0b) | `hypot(&self, other: &Self, p, rm)` (:1287) | L0b `hypot` midpoint, **not** an L1 kernel (PLAN:302) | B (hand-written, not Ziv) | **On the B1 list as signature-unchanged [R]**; no conversion needed beyond what the shared helper change gives it |

Not on the list (no L1 v1 kernel calls them): `fma`, `sin_cos`, `sinh_cosh`, `exp10`, `powi`, `powsi`, `cbrt`, `nth_root`, `reciprocal`, `sqr`, `log(x, base)` (zf's two-argument base form, :904), `rem`, `fused_*`, and every special function.

## 3. The mechanism behind "signatures unchanged" (Rusty's recommendation)

- Jeff approved at 8:25 AM CT Rusty's helper option: **Ziv-loop functions keep a borrowed operand and take one fresh copy per pass through the single audited helper `ziv_pass_operand(&self, p)`**; the approved scan list is pinned **[R, approved by Jeff]**. This replaces the earlier open "X family" question (`&ExactNum` kernels vs input-reconstruction contract) in `/tmp/zf-logs/ziv-loop-shapes.md` and `option-a-evidence.md` §2.2 [F] for the L1-called functions.
- **Layer 1a changes internals only**: the shared `round_exact`/`exact_operand` copy is replaced by the helper per pass; plus the accuracy contract if P14 finds a kernel not rigorous. Neither changes a signature **[R]**.
- **Midpoint evaluated twice (earlier question Q2, RESOLVED [R]).** L1 evaluates the same midpoint in `Down` and `Up` (PLAN:349) and bounds `f'` at it. With the wrappers staying `&self`, L1 does this with plain borrowed calls: zero copies, no bracketing API, no second owned midpoint. A paired `f_bounds(&self, p, cc) -> (lo, hi)` is a possible later **additive** speedup, **not in B1**.

## 4. What stays out of B1 (B2)

From `/tmp/zf-logs/option-a-evidence.md` §4.1 and `zf-pin-items-extract.md` §0 [F]: Option A covers about 150 public methods in layers 0a, 0b, 1b, **1a (ExactNum core)**, 2 arrays/dsp, 3 poly/numerics, 4 ball, 5 complex, 6 specials, 7 macro (`expr!`/`cexpr!`), 8 compare crate, then removal of `impl Clone`, `exact_copy*`, `exact_operand`. Done on `wip/option-a`: neg, abs, frexp, ldexp/scalb, int/fract/ceil/floor, round, add/sub/mul/div (Rule R), sqr, layer 0a clone removals. Not done: layer 1a and layers 2-8.
- **B1:** the section 2 functions' internals (layer 1a for them), P11, P14 for #2-22 (21 kernels, including pow and atan2) and the constants (sqrt is outside P14).
- **B2:** layers 2-8 for all other methods, `impl Clone` removal and the extended clone-scan gate, `fma`/`cbrt`/`nth_root`/`reciprocal` and the specials, and any audit of zf specials (REPLAN 3c: none planned; L2d needs one or its own proofs). P8 (`ulp`/`next_up`/`next_down` on `ExactNum`) is not a blocker; if wanted it is an additive read-only accessor in B2, unless L0b uses it **[R]**.
- **Not part of the split:** Stage A (P1, P2, P3, W1 fixes, P9/P12/P15). P13 (all-modes MPFR differential) is largely in place and is re-verified in the Stage A gate **[R]**. P5's named `MIN_NORMAL_EXPONENT` is dropped, because CONTRACT 9.4 puts the subnormal clamp in the L1 caller **[R]**. `wip/option-a` already has P1 `leading_bit_exponent` (`23ce0c6`), P2 `top_u64_mantissa_and_sticky` (`2f0c807`), P3 `from_mag_parts`/`pow2` (`e876fe9`), 07:51-07:54 today [F: git log]; they are not yet a pinned gate-green SHA.
- **Versioning [R]:** B1 is a **git pin on a gate-green SHA**; the published 2.x release must contain the same signatures. B2 must not change any B1 public item (needed because Accumath pins both crates exactly, PLAN:208 [F]).

## 5. Risks and checks, each way

**If B is not split:** L1 and everything needing L1 (L2a complex, L2d, L3, expression evaluation in L4) wait for all of layers 1a-8 before any pin exists; L1 code can still be written against a stub or mpmath mids.

**If B is split:**
1. *B1 is most of the pivot.* The L1 wrappers sit in macros (`gen_wrapper_ovf_rm_cc`, `_log`, `_arg_rm_cc`, `_exact_rm` for sqrt) over the shared `round_exact`/`exact_operand` machinery [F: ext.rs:2173-2300]; layer 1a is "the long pole" [F: evidence §4.3]. What B2 defers is layers 2-8 and specials, not the pivot.
2. *P14 order.* A written audit of kernels that layer 1a rewrites is "likely to be invalidated" [F]; within B1 the order is 1a for the listed functions, then P14.
3. *More than one pin SHA* (A, B1, final). Pin changes can change grid digests (TEST_PLAN 15 end) [F]; each pin change bumps the op-version per PLAN Part J.
4. *A later layer wants a different shape for a B1 function.* Mitigated by **[R]**: B2 may not change B1 public items, and B1 keeps today's signatures.
5. *L1 stays red on real kernels until B1's P14 findings and MPFR sets are green* (TEST_PLAN 15) [F]. Stub-backed L1 code is useful for parallelism and early API feedback only.
6. *Accuracy-contract change on a frozen signature.* If P14 finds a kernel not rigorous, its contract becomes `ToEven` plus a proven 1 ulp bound (PLAN:195, 349), a documented change that does not alter the signature **[R]**; zenith-bound's radius adds 1 ulp for it.

## 6. Order inside B1 (dependencies only)

1. Jeff's decision 5 on constants (open); everything else in section 2 is not blocked by a decision, subject to Paul's confirmation of the **[R]** items.
2. `ziv_pass_operand` and layer 1a for the section 2 functions (internals), with the approved scan list.
3. P14 findings and MPFR Up/Down sets (256 and 512 bits, every mode, exact-case trap) for #2-22 and the constants; P11 `get_exact_prec`.
4. Gate-green SHA recorded as the B1 git pin; the published 2.x carries the same signatures.
Stage A is independent and goes first.

## 7. Schedule

None. **No measured estimate exists; Rusty will give a real estimate after the Ziv-loop conversion commits land.** The earlier guessed share of Stage B inside B1 (50/60/70%) and every date derived from it are withdrawn and must not be re-added.

## 8. Still to confirm (Paul, with Rusty)

All **[R]** items above (signatures unchanged; list; hypot on B1 signature-unchanged; B1/B2 boundary; git-pin versioning; P5, P8, P13 handling), and Jeff's decision 5 on constants. Questions Q1 (X family: answered by the `ziv_pass_operand` approval) and Q2 (midpoint twice: resolved, section 3) are closed; the call-site-count question is answered provisionally: no L1 call site changes **[R]**.
