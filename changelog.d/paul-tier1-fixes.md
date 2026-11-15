### Fixed (Tier 1, zenith-float-num)

- Subnormals: frexp/ldexp/logb/ilogb use the true exponent of a subnormal (ZF-01, ZF-03); ieee_soft add/fma/from_exact
  round a subnormal result once with the sticky bits jammed before the subtract, div normalizes subnormal significands,
  and fma with an exact zero product and an opposite-sign zero gives +0 (ZF-12..15).
- next_after steps toward and from ±Inf, shrinks the magnitude for a negative value going up, renormalizes below a power
  of two and rounds an unrepresentable self toward the target (ZF-02).
- log with base 0 follows ln x / ln 0 with IEEE division (ZF-04); complex atan on the axes follows C99 catan.
- Complex mul forms the four products exactly and rounds each part once; when an exact product of nonzero operands
  underflows, the directed modes keep the per-product rounding so the least subnormal is reached. Complex exp/sin/cos
  round each part through the Ziv loop (ZF-18).
- Gamma at p ≥ 2048 scales the Stirling term count with p so it stays under the shift cap; trig and Clausen carry the
  full argument exponent in π up to 2^22 bits (ZF-16); atan/atan2 invert directed rounding before negating π (ZF-17, ZF-36).
- Jacobi: K reduction carries the argument exponent bits and am keeps its 2π quasi-period (ZF-60, ZF-23); arccn for
  x < 0 (ZF-62); m = 1 far out of range (ZF-59).
- Special values: distribution zero scale/shape is NaN (ZF-45); bessel_y at x = 0 (ZF-26); elliptic_pi at x = ±1
  (ZF-34); gammainc_upper at s = 0 is E1(x) (ZF-61); barnes_g at non-positive integers and tiny |z| (ZF-77, ZF-74);
  bi, bi_prime, digamma, two_sum, two_product map a finite overflow by the rounding mode (ZF-27, ZF-29).
- Rational round uses exact integer division, ties away from zero (ZF-09); from_i128/from_u128 below 128 bits round
  once to p (ZF-07); rem takes the exponent gap as abs_diff (ZF-55); precision-limit overflows return errors (C2).
- Build: the `random` feature needs std (E2).
- Grid digests complex_mul_div, trig_hyperbolic and jacobi re-pinned for these output changes (jacobi held count 24 -> 0).
- No new `.clone()`: the Tier 1 sites build the value fresh or use a borrowing op, so the clone baseline does not rise.
