#!/usr/bin/env bash
# Self-test of the logdet row cache of scripts/gen_la_goldens.py (Paul 2026-11-11 05:11 CT: the heavy cosmo logdet cases run
# concurrently, one process per case, and the output must be byte-identical to the serial run).
#   --logdet-unit CASE:PART --row-cache DIR  computes one cosmo_* logdet case (PART all|main|enc4096|enc16384) exactly as the
#   serial generator does and writes its row text(s) to DIR; a later normal run given --row-cache DIR takes those rows instead
#   of recomputing them (each cached row must start with the case's own input fields, else it refuses).
# Cases: a unit run writes its 3 rows; the cached main row equals the committed la_logdet.txt row; a --check --skip-heavy run
# with the cache (gate + budget-tier cosmo cases cached) is check: OK and names its cache hits, so the bytes are identical to
# the serial run; a cache row with wrong input fields is refused; a cache row with a wrong value fails the check; bad units exit 2.
# Runs on a copy of scripts/ and tests/goldens/slices/ under $(zb_scratch). Needs the oracle venv and a popslot scope.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap gen_la_rowcache_selftest "$@"   # run log (tools/runlog.sh)
set -u
here=$(cd "$(dirname "$0")" && pwd)
fail() { echo "gen_la_rowcache selftest FAIL: $*"; exit 1; }
. "$here/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/gen-la-rowcache-selftest.XXXXXX") || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
py=${ZB_ORACLE_PY:-$HOME/team-env/oracle-2026-10-30/bin/python}; [[ -x $py ]] || fail "oracle python missing: $py"
mkdir -p "$t/r/tests/goldens" && cp -r "$here/../scripts" "$t/r/" && cp -r "$here/../tests/goldens/slices" "$t/r/tests/goldens/" || fail "copy failed"
gen="$t/r/scripts/gen_la_goldens.py"; rc="$t/rc"; gold="$t/r/tests/goldens/slices/cholesky/la_logdet.txt"
unit() { "$py" "$gen" --logdet-unit "$1" --row-cache "$rc" > "$t/out" 2>&1; echo $?; }
chk() { "$py" "$gen" --check "$t/c$1" --skip-heavy --row-cache "$rc" > "$t/out" 2>&1; echo $?; }
# 1. unit runs: a gate cosmo case and a budget-tier cosmo case, all three parts each
ck; [ "$(unit cosmo_cond.e50_06:all)" = 0 ] || fail "unit cosmo_cond.e50_06:all failed: $(tail -3 "$t/out")"
ck; [ "$(unit cosmo_spread.e300_08:all)" = 0 ] || fail "unit cosmo_spread.e300_08:all failed: $(tail -3 "$t/out")"
for p in main enc4096 enc16384; do ck; [ -s "$rc/logdet.cosmo_cond.e50_06.$p.row" ] || fail "no $p row written: $(ls "$rc" 2>&1)"; done
# 2. the cached main row is the committed serial row, byte for byte
ck; cmp -s <(grep '^logdet\.cosmo_cond\.e50_06 ' "$gold") "$rc/logdet.cosmo_cond.e50_06.main.row" || fail "cached main row differs from $gold"
# 3. a full --check --skip-heavy run through the cache: check: OK (gate and budget files identical) and the hits are named
ck; [ "$(chk 1)" = 0 ] || fail "check with row cache failed: $(grep -v '^progress' "$t/out" | tail -5)"
ck; grep -q '^row-cache: hit logdet.cosmo_cond.e50_06 main' "$t/out" || fail "cache hit not named: $(grep row-cache "$t/out")"
ck; grep -q '^row-cache: hit logdet.cosmo_spread.e300_08 main' "$t/out" || fail "budget-tier cache hit not named"
ck; grep -q '^check: OK cholesky/la_logdet.txt' "$t/out" || fail "la_logdet.txt not OK through the cache"
ck; grep -q '^check: OK heavy cholesky/la_logdet_budget_heavy.txt' "$t/out" || fail "budget file not OK through the cache: $(grep budget "$t/out")"
# 4. a cached row whose input fields are not the case's own is refused (exit nonzero, named)
cp "$rc/logdet.cosmo_cond.e50_06.main.row" "$t/main.orig"
sed 's/^logdet\.cosmo_cond\.e50_06 6 /logdet.cosmo_cond.e50_06 6 7 /' "$t/main.orig" > "$rc/logdet.cosmo_cond.e50_06.main.row"
ck; [ "$(chk 2)" != 0 ] || fail "row with foreign input fields accepted"
ck; grep -q '^gen_la_goldens: row-cache MISMATCH logdet.cosmo_cond.e50_06 main' "$t/out" || fail "mismatch not named: $(tail -3 "$t/out")"
# 5. a cached row with the right inputs but a wrong value is caught by the check (DIFF), never passed
awk '{ $NF = $NF "0"; print }' "$t/main.orig" > "$rc/logdet.cosmo_cond.e50_06.main.row"
ck; [ "$(chk 3)" = 1 ] || fail "row with a wrong value passed the check"
ck; grep -q '^check: DIFF cholesky/la_logdet.txt' "$t/out" || fail "wrong value not reported as DIFF: $(grep check: "$t/out" | tail -3)"
cp "$t/main.orig" "$rc/logdet.cosmo_cond.e50_06.main.row"
# 6. bad units exit 2: unknown case, unknown part, no --row-cache
ck; [ "$(unit cosmo_nope.x:all)" = 2 ] || fail "unknown case not exit 2"
ck; [ "$(unit cosmo_cond.e50_06:enc99)" = 2 ] || fail "unknown part not exit 2"
ck; "$py" "$gen" --logdet-unit cosmo_cond.e50_06:main > "$t/out" 2>&1; [ $? = 2 ] || fail "--logdet-unit without --row-cache not exit 2"
st_done gen_la_rowcache 18
echo "gen_la_rowcache selftest: ok (18 cases)"
