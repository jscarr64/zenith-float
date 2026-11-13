#!/usr/bin/env bash
# Self-test of the parallel (jobs > 1) cond-column path of scripts/la_cond_cols.py and scripts/gen_la_goldens.py (Paul/Jeff
# 2026-11-11, P15 use all cores: the heavy logdet `main` row of cosmo_spread.p2e1000000_lowrank64 ran on one core for hours).
#   la_cond_cols(..., jobs=N) and lu_pp(A, jobs=N) fork N workers over independent rows/columns/entries (|L||U| rows, A^-1
#   columns, |A^-1| G rows with their row sums, lu_pp's L U == P A check rows) and reassemble in index order. The arithmetic is
#   exact (gmpy2.mpq), so every jobs value must give the jobs=1 bytes; jobs=1 (the default) is the old serial code.
# Cases: lu_pp factors and la_cond_cols tokens (logdet and solve kinds) equal across jobs 1, 2, 3, 7, n+5 on n=1, n=2,
# pivoting/permutation, tie, near-singular, spread and low-rank matrices; jobs < 1 is refused; lu_pp with jobs still rejects a
# wrong factorization; generator --jobs 4 rows equal the --jobs 1 rows and the committed la_logdet.txt row; --jobs 0 exits 2.
# Runs on a copy of scripts/ and tests/goldens/slices/ under $(zb_scratch). Needs the oracle venv and a popslot scope.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap gen_la_jobs_selftest "$@"   # run log (tools/runlog.sh)
set -u
here=$(cd "$(dirname "$0")" && pwd)
fail() { echo "gen_la_jobs selftest FAIL: $*"; exit 1; }
. "$here/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/gen-la-jobs-selftest.XXXXXX") || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
py=${ZB_ORACLE_PY:-$HOME/team-env/oracle-2026-10-30/bin/python}; [[ -x $py ]] || fail "oracle python missing: $py"
mkdir -p "$t/r/tests/goldens" && cp -r "$here/../scripts" "$t/r/" && cp -r "$here/../tests/goldens/slices" "$t/r/tests/goldens/" || fail "copy failed"
gen="$t/r/scripts/gen_la_goldens.py"; gold="$t/r/tests/goldens/slices/cholesky/la_logdet.txt"
cat > "$t/eq.py" <<'PY'
import sys, random
sys.path.insert(0, sys.argv[1])
import gen_la_goldens as g
import la_cond_cols as c
from gmpy2 import mpq as F
what = sys.argv[2]
g.rng = random.Random("gen_la_jobs_selftest")
def mats():
    yield "n1", [[F(3)]]
    yield "n1neg", [[F(-5, 7)]]
    yield "n2swap", [[F(0), F(1)], [F(1), F(0)]]
    yield "n2", [[F(1), F(2)], [F(3), F(4)]]
    yield "cyclic4", [[F(int((j - i) % 4 == 1)) + F(int(i == j), 3) for j in range(4)] for i in range(4)]
    yield "tie3", [[F(2), F(1), F(0)], [F(-2), F(3), F(1)], [F(2), F(0), F(5)]]
    yield "near_singular", [[F(1), F(1)], [F(1), F(1) + F(1, 10 ** 30)]]
    yield "spread8", g.ldlt(g.unit_lower(8), g.spread(8, 300, -300, 10))
    d = g.spread(9, 3000, -3000, 2); U = g.imat(9, 2)
    yield "lowrank9", g.lowrank_diag(d, U)
    yield "rand7", [[F(g.rng.randint(-9, 9), g.rng.randint(1, 5)) for _ in range(7)] for _ in range(7)]
bad = 0
if what == "eq":
    for name, A in mats():
        n = len(A); ref = g.lu_pp(A); x = [F(i + 1, 2) for i in range(n)]
        r_ld = c.cond_tokens(c.la_cond_cols(A, *ref, c.scale_row_ln(), "logdet"))
        r_so = c.cond_tokens(c.la_cond_cols(A, *ref, c.scale_row_solve(x), "solve", x))
        for j in (1, 2, 3, 7, n + 5):
            f = g.lu_pp(A, jobs=j)
            if f != ref: bad += 1; print(f"DIFF lu_pp {name} jobs={j}")
            if c.cond_tokens(c.la_cond_cols(A, *f, c.scale_row_ln(), "logdet", jobs=j)) != r_ld: bad += 1; print(f"DIFF logdet {name} jobs={j}")
            if c.cond_tokens(c.la_cond_cols(A, *f, c.scale_row_solve(x), "solve", x, jobs=j)) != r_so: bad += 1; print(f"DIFF solve {name} jobs={j}")
    print("eq: bad", bad)
elif what == "refuse":
    A = [[F(1), F(2)], [F(3), F(4)]]; ref = g.lu_pp(A)
    for j in (0, -1):
        for fn in (lambda: g.lu_pp(A, jobs=j), lambda: c.la_cond_cols(A, *ref, F(1), "logdet", jobs=j)):
            try: fn(); bad += 1; print(f"accepted jobs={j}")
            except ValueError: pass
    print("refuse: bad", bad)
elif what == "badfactor":
    # a wrong factorization must still fail the L U == P A check when the check runs in workers
    real = g.lu_pp.__globals__["mm"]
    g.mm = lambda X, Y: [[v + 1 for v in r] for r in real(X, Y)]
    try:
        g.lu_pp([[F(1), F(2)], [F(3), F(4)]], jobs=2); print("badfactor: accepted")
    except AssertionError: print("badfactor: refused")
sys.exit(1 if bad else 0)
PY
run() { "$py" "$t/eq.py" "$t/r/scripts" "$1" > "$t/out" 2>&1; echo $?; }
# 1. lu_pp and la_cond_cols (logdet + solve) are identical for every jobs value on every matrix
ck; [ "$(run eq)" = 0 ] || fail "jobs results differ: $(tail -5 "$t/out")"
ck; grep -q '^eq: bad 0$' "$t/out" || fail "eq summary missing: $(tail -3 "$t/out")"
# 2. jobs < 1 is refused (ValueError), by lu_pp and by la_cond_cols
ck; [ "$(run refuse)" = 0 ] || fail "jobs < 1 accepted: $(tail -5 "$t/out")"
# 3. the worker-side L U == P A check still refuses a wrong product
ck; run badfactor > /dev/null; grep -q '^badfactor: refused$' "$t/out" || fail "wrong factorization not refused with jobs=2: $(tail -3 "$t/out")"
# 4. generator: --jobs 4 unit rows equal the --jobs 1 rows and the committed serial row; --jobs 0 exits 2
for j in 1 4; do
  ck; "$py" "$gen" --logdet-unit cosmo_spread.e300_12:all --row-cache "$t/rc$j" --jobs $j > "$t/out" 2>&1 || fail "unit --jobs $j failed: $(tail -3 "$t/out")"
done
for p in main enc4096 enc16384; do ck; cmp -s "$t/rc1/logdet.cosmo_spread.e300_12.$p.row" "$t/rc4/logdet.cosmo_spread.e300_12.$p.row" || fail "--jobs 4 $p row differs from --jobs 1"; done
ck; cmp -s <(grep '^logdet\.cosmo_spread\.e300_12 ' "$gold") "$t/rc4/logdet.cosmo_spread.e300_12.main.row" || fail "--jobs 4 main row differs from $gold"
ck; "$py" "$gen" --logdet-unit cosmo_spread.e300_12:main --row-cache "$t/rc0" --jobs 0 > "$t/out" 2>&1; [ $? = 2 ] || fail "--jobs 0 not exit 2"
st_done gen_la_jobs 11
echo "gen_la_jobs selftest: ok (11 cases)"
