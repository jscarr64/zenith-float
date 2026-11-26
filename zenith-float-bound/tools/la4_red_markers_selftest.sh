#!/usr/bin/env bash
# Self-test of tools/la4_red_markers.sh on fixture trees: markers with a CONTRACT item pass, a marker without one fails, LA4 marked done
# with an LA4 marker left fails, and "done" elsewhere (another slice's section, or markers naming only other slices) does not trip it.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap la4_red_markers_selftest "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
here="$(cd "$(dirname "$0")" && pwd)"; . "$here/scratch.sh"
fail() { echo "la4_red_markers selftest FAIL: $*"; exit 1; }
t=$(mktemp -d "$(zb_scratch)/la4-red-selftest.XXXXXX") || fail "mktemp failed"; runlog_on_exit "rm -rf \"$t\""; . "$here/selftest_cases.sh"; st_init "$t"
[ -x "$here/la4_red_markers.sh" ] || fail "tools/la4_red_markers.sh missing or not executable"
# fx NAME LA4_STATUS_LINE OTHER_STATUS_LINE [MARKER_TEXT...]: a fixture tree with docs/CONTRACT.md and one test file per marker
fx() {
  local d="$t/$1" i=0; mkdir -p "$d/docs" "$d/tests/harness"
  printf '## Slices\n### Slice LA2\n- **LA2.3** lu rows.\n### Slice LA3 (`src/qr.rs`)\n%s\n### Slice LA4 (`src/gls.rs`)\n%s\n- **LA4.1** gls rows.\n- **LA4.2** logdet rows; see also LA4.7 (prose mention only).\n### Slice LA5\n_Empty._\n### Slice S2 (`src/bridge.rs`)\n- **16.S2.1 cmp_exact_num.** cmp rows.\n- **16.S2.2 cmp_exact_diff.** diff rows; see also 16.S2.7 (prose mention only).\n' "$3" "$2" > "$d/docs/CONTRACT.md"
  shift 3
  for m in "$@"; do i=$((i + 1)); printf '#[test]\n#[ignore = "%s"]\nfn t%s() {}\n' "$m" "$i" > "$d/tests/harness/t$i.rs"; done
}
run() { "$here/la4_red_markers.sh" "$t/$1" > "$t/out" 2>&1; echo $?; }
G1='red: not implemented: zb::gls (CONTRACT LA4.1: gls rows)'; G2='red: not implemented: zb logdet (CONTRACT LA4.2: logdet rows)'
fx clean '_Empty._' '_Empty._'; ck; [ "$(run clean)" = 0 ] || fail "clean tree failed: $(cat "$t/out")"
fx items '_Empty._' '_Empty._' "$G1" "$G2"; ck; [ "$(run items)" = 0 ] || fail "markers with items failed: $(cat "$t/out")"
ck; grep -q '2 LA4 red marker' "$t/out" || fail "marker count not reported: $(cat "$t/out")"
fx noitem '_Empty._' '_Empty._' 'red: not implemented: zb::gls'; ck; [ "$(run noitem)" = 1 ] || fail "marker without a CONTRACT item passed"
ck; grep -q 'without its CONTRACT' "$t/out" || fail "missing item not named: $(cat "$t/out")"
fx done 'Status: done (2026-11-01)' '_Empty._' "$G1"; ck; [ "$(run done)" = 1 ] || fail "LA4 done with a red marker passed"
ck; grep -q 'marked done' "$t/out" || fail "done violation not named: $(cat "$t/out")"
fx donebold '**Status: DONE**' '_Empty._' "$G2"; ck; [ "$(run donebold)" = 1 ] || fail "bold DONE with a red marker passed"
fx donegreen 'Status: done' '_Empty._'; ck; [ "$(run donegreen)" = 0 ] || fail "LA4 done with no markers failed: $(cat "$t/out")"
fx otherdone '_Empty._' 'Status: done' "$G1"; ck; [ "$(run otherdone)" = 0 ] || fail "LA3 done tripped the LA4 check: $(cat "$t/out")"
fx otherslice 'Status: done' '_Empty._' 'red: not implemented: zb lu (CONTRACT LA2.3: lu rows)'; ck; [ "$(run otherslice)" = 0 ] || fail "an LA2 marker blocked LA4 done: $(cat "$t/out")"
# Spence (2026-10-24): every CONTRACT id a marker names must exist as a **<id>** entry in docs/CONTRACT.md (a prose mention is not an entry)
fx noentry '_Empty._' '_Empty._' 'red: not implemented: zb::gls (CONTRACT LA4.9: gls rows)'; ck; [ "$(run noentry)" = 1 ] || fail "a marker naming a CONTRACT id with no entry passed"
ck; grep -q 'LA4.9 has no entry in docs/CONTRACT.md' "$t/out" || fail "missing entry not named: $(cat "$t/out")"
fx prose '_Empty._' '_Empty._' 'red: not implemented: zb::gls (CONTRACT LA4.7: gls rows)'; ck; [ "$(run prose)" = 1 ] || fail "a prose-only mention counted as a CONTRACT entry"
# Harper (2026-10-29): `red: CONTRACT <id>: ...` markers (slice clause ids 16.<SLICE>.<n>, e.g. 16.S2.1) must be well-formed and have a
# `- **<id> <title>.**` entry in docs/CONTRACT.md (the form the 16.x clauses use); an unknown or malformed id fails.
C1='red: CONTRACT 16.S2.1: cmp rows'; C2='red: CONTRACT 16.S2.2: diff rows'
fx cgood '_Empty._' '_Empty._' "$C1" "$C2" "$G1"; ck; [ "$(run cgood)" = 0 ] || fail "known 16.S2.x CONTRACT markers failed: $(cat "$t/out")"
ck; grep -q '1 LA4 red marker' "$t/out" || fail "LA4 count changed by CONTRACT markers: $(cat "$t/out")"
fx cunknown '_Empty._' '_Empty._' 'red: CONTRACT 16.S2.9: no such clause'; ck; [ "$(run cunknown)" = 1 ] || fail "a CONTRACT marker with an unknown id passed"
ck; grep -q '16.S2.9 has no entry in docs/CONTRACT.md' "$t/out" || fail "unknown id not named: $(cat "$t/out")"
fx cprose '_Empty._' '_Empty._' 'red: CONTRACT 16.S2.7: prose only'; ck; [ "$(run cprose)" = 1 ] || fail "a prose-only 16.S2.7 mention counted as an entry"
fx cmalformed '_Empty._' '_Empty._' 'red: CONTRACT S2.1: missing the 16. prefix'; ck; [ "$(run cmalformed)" = 1 ] || fail "a malformed CONTRACT id passed"
ck; grep -q 'malformed' "$t/out" || fail "malformed id not named: $(cat "$t/out")"
fx cnocolon '_Empty._' '_Empty._' 'red: CONTRACT 16.S2.1x rows'; ck; [ "$(run cnocolon)" = 1 ] || fail "a CONTRACT id without its colon passed"
mkdir -p "$t/nocontract/tests"; ck; [ "$(run nocontract)" = 2 ] || fail "missing CONTRACT.md not a setup error"
st_done la4_red_markers 23   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "la4_red_markers selftest: ok"
