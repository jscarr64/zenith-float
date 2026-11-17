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
- Gamma at p ≥ 2048 scales the Stirling term count with p so it stays under the shift cap (BROKEN: 361 of 385 r.gamma
  rows at p ≥ 2048 pass; 24 still return NaN(PrecisionRetryExhausted), and the Γ family has 1,047 of 9,057 rows still
  failing, deferred to Tier 2); sin/cos/tan carry the full argument exponent in π with no cap, the width limited only by
  representability and allocation (ZF-16 BROKEN: 222 of 588 rows still fail, in rem_pi and Clausen, deferred);
  atan/atan2 invert directed rounding before negating π (ZF-17, ZF-36).
- Jacobi: K reduction carries the argument exponent bits with no cap (ZF-60 BROKEN: 317 of 590 rows still fail, in cd,
  dc, sn, ns, dn; deferred) and am keeps its 2π quasi-period (ZF-23); arccn for
  x < 0 (ZF-62); m = 1 far out of range (ZF-59).
- Special values: distribution zero scale/shape is NaN (ZF-45); bessel_y at x = 0 (ZF-26); elliptic_pi at x = ±1
  (ZF-34); gammainc_upper at s = 0 is E1(x) (ZF-61); barnes_g at non-positive integers and tiny |z| (ZF-77, ZF-74);
  bi, bi_prime, digamma, two_sum, two_product map a finite overflow by the rounding mode (ZF-27, ZF-29).
- Rational round uses exact integer division, ties away from zero (ZF-09); from_i128/from_u128 below 128 bits round
  once to p (ZF-07); rem takes the exponent gap as abs_diff and its modular power keeps bit 63 of the gap, so rem(±MAX, minpos) is an exact ±0
  (ZF-55, 5 of 5 rows); precision-limit overflows return errors (C2).
- Build: the `random` feature needs std (E2).
- Grid digests complex_mul_div, trig_hyperbolic and jacobi re-pinned for these output changes (jacobi held count 24 -> 0).
- No new `.clone()`: the Tier 1 sites build the value fresh or use a borrowing op, so the clone baseline does not rise.
- 2F1 outside the unit disc via DLMF 15.8.8; z = 1 with no finite limit is ±Inf (ZF-44, ZF-67). Terminating 2F1 is summed
  exactly by binary splitting at any degree and rounded once (ZF-70); size limit removed per Jeff's order 2026-11-15;
  the 8 TiB (2^46-bit) guard was Paul's addition, mislabelled as a ruling, removed 2026-11-16. A real allocation
  failure is still Err(MemoryAllocation) from try_reserve. Exact complex parts for exp2/exp10/log2/log10/nth_root/cbrt/sqrt on axes (ZF-70).
- ln_gamma: Stirling at any argument width for z >= 2^47 (fixes the sign of zeta(−(2^62+1)), ZF-11); tiny x, x < 1/2
  and x near 1 or 2 without overflow or absolute-error loss; ln_gamma(1) = ln_gamma(2) = +0 (ZF-30).
- gammainc_upper: DLMF 8.11.2 asymptotic for x > s, measured cancellation otherwise (ZF-25). Si/Ci series carry their
  cancellation bits (ZF-24). Fresnel angle uses exact x² mod 4 (ZF-32). Lambert W0/W−1 correctly rounded with a
  certified bracket; x < −1/e is NaN; W−1(±0) is NaN(InvalidArgument). Jeff wrote "Lambert W staying NaN at zero";
  6227d48's −∞ was a mislabelled 2026-11-15 ruling and is reverted (ZF-31).
- Binary format: wide-exponent heap and inline records (ZF-08); a value the 16-byte inline record cannot hold is
  InvalidArgument (format limit) (Jeff ruling 2026-11-16 8:18 PM; previously Paul's call, mislabelled as a 2026-11-15 ruling); to_bytes encodes every finite value.
- ZF-51: gk_adapt reports success only when the whole error estimate, final rounding included, meets tol, so Ok only when
  the bound is within tol (Jeff ruling 2026-11-16 8:18 PM; previously Paul's call, mislabelled as a 2026-11-15 ruling).
- ZF-78 (partial): r.format of max_value in radix 3/36 no longer hangs. conv_to_general
  estimates the radix exponent and scales with powsi instead of dividing once per unit
  (~10^18 iterations at EXPONENT_MAX). convert_from uses the same scale; format retries
  a short correctly rounded spelling so a directed parse round-trips. Test
  zf_78_format_max_value_general_radix.
- ZF-78 (partial): Bessel/Kelvin and complex gamma/digamma at p=2048 no longer hang.
  Kelvin is Re/Im of complex J/K. Those (and the later gamma/digamma cluster) hung
  because complex Stirling used a fixed 64-term tail, forcing ~2^{p/127} shifts
  (about a million at p=2048). The term count now scales with p as the real kernel
  already did. Integer-order complex K uses the K_0 series instead of J±iY
  (ker/kei of 700). Real Y_n uses Hankel Y_0/Y_1 plus recurrence when |x| is in
  the Hankel range but |ν|² ≥ 2|x| (Y_100(1000)). Test
  zf_78_bessel_kelvin_high_prec. I_{2+i}(-2±0i) and J/K/Y of the same cut
  were NaN(PrecisionRetryExhausted): complex pow of a real or imaginary
  axis base with a non-real exponent went through exp(w Log z) and could
  not certify an exact-axis result ((-1)^{2+i}=e^{-π}, i^{2+i}=-e^{-π/2}).
  pow_real_axis now uses the exact argument θ=mπ/2 and forces exact zero
  parts. kei_100(2^{-30}) was NaN(InvalidArgument) because
  BESSEL_INTEGER_MAX=64 refused a computable order (P9); integer orders
  that fit in i32 are recognized, and the K recurrence shares
  BESSEL_K_RECUR_MAX. ker/kei of −700 hung because integer-order K of a
  left-half argument skipped the K_0 series; Re z<0 now uses DLMF 10.34.2
  to reduce to K and I of −z. Tests zf_78_bessel_cut_and_large_integer_order
  and zf_78_bessel_kelvin_high_prec. BESSEL_K_RECUR_MAX=1<<22 is removed:
  the K recurrence scales with ν (K_{2^{22}+1}(1) at p=64 in ~3.4s).
  Test zf_78_k_recur_scales_with_nu. On the cut, J, Y, I, and K follow
  C99: Im z = +0 is the limit from above (Arg = +π), Im z = −0 from
  below (Arg = −π). FLINT's principal value is the +0 side; −0 is
  judged against the DLMF 10.11 / 10.27 / 10.34 continuation. Test
  zf_78_bessel_cut_and_large_integer_order. The 2F1 cluster remains open.
