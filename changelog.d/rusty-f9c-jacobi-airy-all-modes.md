### Jacobi and Airy in all seven modes (Spence spot-check at cc85211, items 2 and 3)

The new `tests/mpfr/all_modes_jacobi_airy.rs` checks every mode. The digests in ext.rs still pin only ToEven and
Down; those old pins do not move.
- **Ai:** checked strictly against MPFR `mpfr_ai`.
- **All thirteen Jacobi functions, Ai, Bi, Ai' and Bi':** checked against `bracket_expect`. It takes the function at
  rp(p, 64) + 128 bits in Down and Up, then lets MPFR round the bracket's midpoint (or the exact value) in each mode.
  The helper has its own tests against MPFR (exact, inexact, an exact tie at p, a non-word p).
- **Inputs:** the inputs are the digest's operands, which avoid Bi and Bi' overflow. Under None, ToZero and
  truncating Up/Down, those overflows are the ZF-M27 F-2 rows.
- **Held:** the tiny Jacobi operand 21 is HELD for new findings N-6 (directed modes give NaN(PrecisionRetryExhausted))
  and N-7 (None gives 1 above values below 1, and the flag is wrong).
