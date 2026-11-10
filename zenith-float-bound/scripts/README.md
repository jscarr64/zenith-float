# scripts

## gen_goldens.py

Generates `tests/goldens/mag_ops.txt`, the mpmath golden rows for the integer-only `Mag` operations (`add_up`, `mul_up`, `mul_lower`, `sub_lower`, `div_up`, `div_lower`, `sqrt_up`, `sqrt_lower`). The rows come from `mpmath.libmp` directed rounding to 64 bits, an oracle that is independent of zenith-float and of the Rust model. Every row is also checked with exact rationals (`fractions.Fraction`) before it is written, when its exponents are close enough to compute; the script prints how many rows were checked and how many were skipped.

Run it with the **team env** (since 2026-10-30, Jeff: never `~/.local/oracle-venv`; goldlib refuses it and any path that is not `~/team-env/<name>` passing `~/team-env/verify_env.sh <name>`, and records the env name in the log and in each output dir's `PROVENANCE`) `~/team-env/oracle-2026-10-30/bin/python` (every generator has it as its shebang):

    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_goldens.py

Oracle venv (since 2026-10-23, Jeff's decision): mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0, the pins of zenith-float `scripts/oracle-requirements.lock` (recorded in `requirements.lock`). `goldlib.require_oracle` runs when any generator imports `goldlib` and exits nonzero unless all four versions match exactly, `sys.prefix` is the venv and no `~/.local` path outside the venv is on `sys.path` (the `flint` in `~/.local` is an unrelated tool). `/usr/bin/python3` (Ubuntu mpmath 1.2.1) is refused. Goldens regenerated in the venv were byte-identical to the 1.2.1 ones apart from the provenance header line.

The output is deterministic (fixed operand sets, a fixed xorshift seed, order-preserving de-duplication), so regenerating must leave `tests/goldens/mag_ops.txt` byte-identical (a gate step that regenerates and diffs it is planned, shared-harness item S8; it does not exist yet).

Operand sets: binary rows use four mantissas and eleven exponents (including `i64::MIN`, `i64::MAX` and their neighbours) plus a wide block of twelve mantissas at the exponent extremes; unary rows use twelve mantissas (the same list as `tests/model/sets.rs`) and perfect or near-perfect squares. Round-down overflow rows expect `Finite{i64::MAX, u64::MAX}`; round-up overflow rows expect `inf`.

## gen_review_goldens.py

Generates `tests/goldens/mag_ops_review.txt`: boundary rows added after the L0a review (sub_lower at exponent gap 191, division with remainder exactly 1 and mb - 1, add_up and sub_lower at the alignment gaps around 64, 128 and 192, mul carry rows at the exponent limits, `ulp_mag` and `half_ulp` at the `i32` exponent limits and large precisions). The oracle is a big-integer script inside the generator (no mpmath, no Rust code), but it still imports `goldlib`, so it runs only in the oracle venv like the others:

    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_review_goldens.py > tests/goldens/mag_ops_review.txt

The output is deterministic. `tests/mag/review_rows.rs` runs every row against the crate.

## gen_ball_goldens.py

Generates `tests/goldens/ball_exact_ops.txt` (L0b-7): rows `OP PREC RM IPREC A [B] R INEXACT` for `add`, `sub`, `mul`, `div` and `sqrt` of EXACT inputs, rounded to `PREC` bits under each of the six rounding modes, at the PS-fast precisions 1, 2, 3, 4, 8, 24, 31, 32, 33, 53, 63, 64, 65, 127, 128, 129, 256 and 1024. The row format (value fields `SIGN EXP MANT_HEX`, zero-sign rules, tie rule, sort order) is documented in the file header. Independent of any zenith-float pin: the model is an unbounded exponent range and the results are plain rounded exact values.

Method: every row is computed with exact `Fraction`s and `goldlib.round_fraction`, and cross-checked against `mpmath.libmp` (`mpf_add`, `mpf_sub`, `mpf_mul`, `mpf_div`, `mpf_sqrt`) for the five modes mpmath has (ToEven, Up, Down, ToZero, FromZero, by value and inexact flag); ToOdd is cross-checked through mpmath's down, up and nearest neighbours; rows with an exact zero result are checked against `goldlib` only (mpmath has no signed zero). The script exits 3 on any disagreement and 4 on a row whose inputs are not exactly representable at their stated precision (`IPREC`). It prints the row, cross-check and byte counts.

Operand sets come from `goldlib.SplitMix64` streams of suite `l0b_ball` (same BASE and seeds table as the Rust side); the exponent bands, mantissa patterns and tie patterns mirror `tests/common/generators.rs` in spirit, NOT bit for bit. Structured sets (zeros, carries, exact ties at `p`, 81/64 at p = 3, cancellation, huge exponent gaps, perfect squares) plus random sets are all used under all six modes.

    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_ball_goldens.py            # writes tests/goldens/ball_exact_ops.txt
    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_ball_goldens.py --check DIR # regenerates and compares with DIR/ball_exact_ops.txt
    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_ball_goldens.py --selftest

The Rust side (`tests/common/ballrows.rs`, `tests/harness/ball_goldens_selftest.rs`) re-verifies every row with an independent exact-integer checker: full recomputation for `PREC <= 64` (u128), structural checks only above.

The gate wrappers (`tools/full_gate.sh`, `tools/goldens_regen.sh`) export `PYTHONNOUSERSITE=1` and call the oracle venv (`GOLDENS_PYTHON` may override it, but the startup check still applies). `tools/goldens_regen.sh` also regenerates `mag_ops_review.txt` and compares it byte for byte.

## gen_ball_enclosure.py

L0b-7 ball enclosure rows (TEST_PLAN section 9) at the full schedule PS: truth interval rounded outward at W = 2*PREC + 64 and R_ref, for
containment checks only (tau and lambda are D2, not asserted). Exact Fractions, cross-checked with mpmath.iv.

    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_ball_enclosure.py            # writes tests/goldens/ball_enclosure_rows.txt
    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_ball_enclosure.py --check DIR # regenerates and compares with DIR/ball_enclosure_rows.txt
    ~/team-env/oracle-2026-10-30/bin/python scripts/gen_ball_enclosure.py --selftest
