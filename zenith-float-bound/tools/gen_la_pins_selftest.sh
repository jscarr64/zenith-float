#!/usr/bin/env bash
# Self-test of the regular-gate heavy check of scripts/gen_la_goldens.py (Spence, 2026-10-24): `--check DIR --skip-heavy` never
# recomputes heavy rows; it logs `check: heavy mode=pins` and verifies each slice's HEAVY.SHA256 by pin/sha only: every entry is a
# strict `<64 hex>  <name>` line and the entry names are exactly the heavy files the generator defines (none missing, none stale).
# A full `--check` (no --skip-heavy) logs `check: heavy mode=recompute`; tools/goldens_regen.sh runs it only at merge points.
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
st_done gen_la_pins 10   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "gen_la_pins selftest: ok"
