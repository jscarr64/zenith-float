#!/usr/bin/env bash
# Self-test of the regular-gate heavy check of scripts/gen_la_goldens.py (Spence, 2026-10-24): `--check DIR --skip-heavy` never
# recomputes heavy rows; it logs `check: heavy mode=pins` and verifies each slice's HEAVY.SHA256 by pin/sha only: every entry is a
# strict `<64 hex>  <name>` line and the entry names are exactly the heavy files the generator defines (none missing, none stale).
# `--check DIR` without --heavy is the same pins check (case 6); `--check DIR --heavy HDIR` logs `check: heavy mode=recompute`; tools/goldens_regen.sh runs it only at merge points.
# Runs on a copy of scripts/ and tests/goldens/slices/ under $(zb_scratch). Needs the oracle venv and a popslot scope; run via poprun/popslot.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap gen_la_pins_selftest "$@"   # run log (tools/runlog.sh)
set -u
here=$(cd "$(dirname "$0")" && pwd)
fail() { echo "gen_la_pins selftest FAIL: $*"; exit 1; }
. "$here/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/gen-la-pins-selftest.XXXXXX") || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
py=${ZB_ORACLE_PY:-$HOME/team-env/oracle-2026-10-30/bin/python}; [[ -x $py ]] || fail "oracle python missing: $py"
mkdir -p "$t/r/tests/goldens" && cp -r "$here/../scripts" "$t/r/" && cp -r "$here/../tests/goldens/slices" "$t/r/tests/goldens/" || fail "copy failed"
h="$t/r/tests/goldens/slices/gls/HEAVY.SHA256"; cp "$h" "$t/gls.orig"
run() { "$py" "$t/r/scripts/gen_la_goldens.py" --check "$t/c$1" --skip-heavy > "$t/out" 2>&1; echo $?; }
# 1. untouched pins: pass, mode logged, pins named, and no heavy row computed (no tier=heavy progress, no partial/heavy file)
ck; [ "$(run 1)" = 0 ] || fail "untouched pins failed: $(tail -5 "$t/out")"
ck; grep -qx 'check: heavy mode=pins (HEAVY.SHA256 pin/sha check only, no recompute)' "$t/out" || fail "pins mode not logged: $(head -3 "$t/out")"
ck; grep -q '^check: OK pin gls/la_gls_heavy.txt ' "$t/out" || fail "gls pin not checked: $(grep pin "$t/out")"
ck; grep -q 'tier=heavy\|mode=recompute' "$t/out" && fail "pins mode recomputed heavy rows"
# 2. a missing pin fails, naming it
grep -v la_gls_heavy "$t/gls.orig" > "$h"
ck; [ "$(run 2)" = 1 ] || fail "missing pin passed"
ck; grep -q '^check: HEAVY.SHA256 pin MISSING gls/la_gls_heavy.txt' "$t/out" || fail "missing pin not named: $(grep -i heavy "$t/out")"
# 3. a malformed pin fails
{ grep '^#' "$t/gls.orig"; echo 'e6b771a5  la_gls_heavy.txt'; } > "$h"
ck; [ "$(run 3)" = 1 ] || fail "malformed pin passed"
ck; grep -q '^check: HEAVY.SHA256 malformed line gls: ' "$t/out" || fail "malformed pin not named: $(grep -i heavy "$t/out")"
# 4. a stale pin (no such heavy file) fails
{ cat "$t/gls.orig"; printf '%064d  la_gls_old_heavy.txt\n' 0; } > "$h"
ck; [ "$(run 4)" = 1 ] || fail "stale pin passed"
ck; grep -q '^check: HEAVY.SHA256 pin STALE gls/la_gls_old_heavy.txt' "$t/out" || fail "stale pin not named: $(grep -i heavy "$t/out")"
# 5. write mode keeps foreign SHA256SUMS lines (e.g. Rusty's la_lstsq.txt) byte-identical and in place, rewrites only the
#    lines for files it generates, and inserts a missing owned line in sorted position (here: la_solve.txt after la_lstsq.txt)
s="$t/r/tests/goldens/slices/lu/SHA256SUMS"
{ grep '^#' "$s"; grep '  la_lu.txt$' "$s"; printf '%064d  la_lstsq.txt\n' 7; grep '  la_solve.txt$' "$s"; echo '# foreign trailer kept'; printf '%064d  zz_rusty.txt\n' 9; } > "$t/lu.want"
grep -v '  la_solve.txt$' "$t/lu.want" > "$s"
ck; "$py" "$t/r/scripts/gen_la_goldens.py" --skip-heavy > "$t/out" 2>&1 || fail "write mode failed: $(tail -5 "$t/out")"
ck; grep -q '  la_lstsq.txt$' "$s" || fail "foreign line la_lstsq.txt dropped: $(cat "$s")"
ck; cmp -s "$t/lu.want" "$s" || fail "SHA256SUMS not merged in place: $(diff "$t/lu.want" "$s")"
# 6. a plain `--check DIR` (no --skip-heavy, no --heavy) must be passable (Cassie, 2026-10-29: on d9cfacb it recomputed
#    logdet.cosmo_spread.p2e100000_lowrank32 for 2877 s against the 1800 s case budget): it checks heavy rows by HEAVY.SHA256
#    pin only; recompute needs an explicit --heavy DIR. --case-budget 30 makes a recompute fail fast (exit 5) instead of hanging.
cp "$t/gls.orig" "$h"
ck; "$py" "$t/r/scripts/gen_la_goldens.py" --check "$t/c6" --case-budget 30 > "$t/out" 2>&1; r=$?; [ "$r" = 0 ] || fail "plain --check exited $r: $(grep -m3 -i 'heavy\|budget' "$t/out")"
ck; grep -qx 'check: heavy mode=pins (HEAVY.SHA256 pin/sha check only, no recompute)' "$t/out" || fail "plain --check not in pins mode: $(head -1 "$t/out")"
ck; grep -q '^check: OK pin gls/la_gls_heavy.txt ' "$t/out" || fail "plain --check did not check the gls pin"
st_done gen_la_pins 16   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "gen_la_pins selftest: ok"
