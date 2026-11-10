#!/usr/bin/env bash
# Self-test of the format-2 LA gold columns and the Arb cross-check of scripts/gen_la_goldens.py (Harper 2026-10-29, Spence's spec
# from Cassie's review). Cases: the python-flint/gmpy2 pins are enforced; a clean ln|det| enclosure passes the Arb cross-check; a
# wrong Arb bracket (GEN_LA_ARB_FAULT=shift selftest hook) MUST fail with ArbMismatch; scripts/la_cond_cols.py gives the exact
# cond_row, k_row, scale_row, containment_only/reason and rref on hand-checked matrices.
# Needs the oracle venv and a popslot scope; run via poprun/popslot.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap gen_la_cond_selftest "$@"   # run log (tools/runlog.sh)
set -u
here=$(cd "$(dirname "$0")" && pwd)
fail() { echo "gen_la_cond selftest FAIL: $*"; exit 1; }
. "$here/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/gen-la-cond-selftest.XXXXXX") || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
py=${ZB_ORACLE_PY:-$HOME/.local/oracle-venv/bin/python}; [[ -x $py ]] || fail "oracle python missing: $py"
s="$here/../scripts"
pyc() { (cd "$t" && "$py" -c "import sys; sys.path.insert(0, '$s'); $1") > "$t/out" 2>&1; echo $?; }
G="import gen_la_goldens as g; F = g.F"
# 1. pins: flint and gmpy2 are pinned in the generator like sympy/mpmath
ck; [ "$(pyc "$G; assert g.FLINT_PIN == '0.9.0' and g.GMPY2_PIN == '2.3.2'; print(g.flint.__version__)")" = 0 ] || fail "pins: $(tail -3 "$t/out")"
# 2. clean enclosures pass the Arb cross-check (256 and 4096 bits, a negative and a huge det)
ck; [ "$(pyc "$G; [g.enclose_log(F(d), b) for d in (F(3, 7), F(-5, 2), F(10) ** 300 + 1) for b in (256, 4096)]")" = 0 ] || fail "clean cross-check failed: $(tail -3 "$t/out")"
# 3. mismatch MUST fail: a shifted Arb bracket raises ArbMismatch (never a silent pass)
ck; [ "$(GEN_LA_ARB_FAULT=shift pyc "$G
try:
    g.enclose_log(F(3, 7), 256)
except g.ArbMismatch as e:
    print('ArbMismatch', e); sys.exit(0)
sys.exit(9)")" = 0 ] || fail "arb mismatch did not fail: $(tail -3 "$t/out")"
ck; grep -q '^ArbMismatch arb bracket not inside' "$t/out" || fail "mismatch reason not named: $(tail -3 "$t/out")"
C="from la_cond_cols import *; from gmpy2 import mpq as F
def lu(A):
    import gen_la_goldens as g; return g.lu_pp(A)"
# 4. identity 2x2: cond 1, k_row = 4 + 1 + 0 = 5, tight, logdet rref = n*tr = 4, scale 1
ck; [ "$(pyc "$C
A = [[F(1), F(0)], [F(0), F(1)]]; p, L, U = lu(A)
t = cond_tokens(la_cond_cols(A, p, L, U, scale_row_ln(), 'logdet'))
assert t == ['1/1', '5', '1/1', '0', 'tight', '4/1'], t")" = 0 ] || fail "identity cols: $(tail -3 "$t/out")"
# 5. [[2,1],[1,1]] solve with x = [1,-1]: A^-1 = [[1,-1],[-1,2]], L = [[1,0],[1/2,1]], U = [[2,1],[0,1/2]], |L||U| = [[2,1],[1,1]],
#    |A^-1| G = [[3,2],[4,3]] -> cond 7, k_row = 4 + 1 + 3 = 8; scale ||x|| = 1; rref = 2 * max(5, 7) = 14
ck; [ "$(pyc "$C
A = [[F(2), F(1)], [F(1), F(1)]]; x = [F(1), F(-1)]; p, L, U = lu(A)
t = cond_tokens(la_cond_cols(A, p, L, U, scale_row_solve(x), 'solve', x))
assert t == ['7/1', '8', '1/1', '0', 'tight', '14/1'], t")" = 0 ] || fail "2x2 solve cols: $(tail -3 "$t/out")"
# 6. near-singular [[1,1],[1,1+2^-200]]: k_row + 8 >= 192 -> containment_only 1 k_row_too_large
ck; [ "$(pyc "$C
A = [[F(1), F(1)], [F(1), 1 + F(1, 2 ** 200)]]; p, L, U = lu(A)
c = la_cond_cols(A, p, L, U, scale_row_ln(), 'logdet')
assert c.k_row == 4 + 1 + ceil_log2(c.cond_row) and c.k_row >= 184 and (c.containment_only, c.reason) == (1, 'k_row_too_large'), c[1:5]")" = 0 ] || fail "near-singular cols: $(tail -3 "$t/out")"
# 7. ceil_log2 exact at and around powers of two; lstsq scale wrapper is ||x||_inf
ck; [ "$(pyc "$C
assert [ceil_log2(F(v)) for v in (1, 8, 9, F(1, 3), F(1, 4), F(7, 2))] == [0, 3, 4, -1, -2, 2]
assert scale_row_lstsq([F(-3, 2), F(1)]) == F(3, 2)")" = 0 ] || fail "ceil_log2/lstsq: $(tail -3 "$t/out")"
# 8. a singular input is refused (no cond columns on ERR rows)
ck; [ "$(pyc "$C
A = [[F(1), F(2)], [F(2), F(4)]]; p, L, U = lu(A)
try:
    la_cond_cols(A, p, L, U, scale_row_ln(), 'logdet')
except AssertionError:
    sys.exit(0)
sys.exit(9)")" = 0 ] || fail "singular accepted: $(tail -3 "$t/out")"
# 10. lstsq (16.LA.5, Rusty): A = [[1],[1]]. Residual b = [1,3]: x = 2, r = [-1,1], |A^+||A| = 1, |(A^TA)^-1||A^T||r| = 1,
#     cond = (1*2 + 1)/2 = 3/2, k_row = 4 + 1 + 1 = 6, rref = 2*(2 + 1) = 6. Exact fit b = [2,2] (r = 0) reduces to the first term:
#     cond 1, k_row 5, rref 2*2 = 4. x = 0 (A = [[1],[-1]], b = [1,1], r = b != 0): lstsq_x_zero, cond from the first term (1).
ck; [ "$(pyc "from la_cond_cols import *; from gmpy2 import mpq as F
A = [[F(1)], [F(1)]]
t = cond_tokens(la_cond_cols_lstsq(A, [F(1), F(3)], [F(2)])); assert t == ['3/2', '6', '2/1', '0', 'tight', '6/1'], t
t = cond_tokens(la_cond_cols_lstsq(A, [F(2), F(2)], [F(2)])); assert t == ['1/1', '5', '2/1', '0', 'tight', '4/1'], t
t = cond_tokens(la_cond_cols_lstsq([[F(1)], [F(-1)]], [F(1), F(1)], [F(0)])); assert t == ['1/1', '5', '0/1', '2', 'lstsq_x_zero', '2/1'], t
try:
    la_cond_cols_lstsq(A, [F(1), F(3)], [F(1)])
except AssertionError:
    sys.exit(0)
sys.exit(9)")" = 0 ] || fail "lstsq cond: $(tail -3 "$t/out")"
st_done gen_la_cond 10   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "gen_la_cond selftest: ok"
