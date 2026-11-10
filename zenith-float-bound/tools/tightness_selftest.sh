#!/usr/bin/env bash
# Self-test for tools/tightness_report.py and tools/tightness_check.sh (S12): builds fixtures under /tmp and proves the ratchet PASSES on an
# unchanged or improved result and FAILS on every failure mode (worse maximum, unrecorded cell, cell that disappeared, malformed results,
# malformed baseline), that --update never raises a baseline, and that the report prints the per-op/per-precision table.
# The RESULTS fixture below is byte-identical to the text pinned by tests/harness/tightness_selftest.rs (the Rust writer): format agreement.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap tightness_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
export LC_ALL=C PYTHONNOUSERSITE=1
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d "$(zb_scratch)"/tightness-selftest.XXXXXX)"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
tool=("/usr/bin/python3" "$here/tightness_report.py")
fail() { echo "tightness SELFTEST FAIL: $*"; exit 1; }
expect_rc() { # expect_rc <rc> <label> <cmd...>
  local want="$1" label="$2"; shift 2
  out=$("$@" 2>&1); rc=$?
  ck; [ "$rc" = "$want" ] || fail "$label: exit $rc, want $want: $out"
}
RESULTS_TEXT='# zenith-bound tightness results (S12), written by the tightness sets; ratio = rad / R_ref in units of 1/65536 (q16), rounded up
# columns: OP PREC COUNT MIN MEDIAN P99 MAX
add 24 3 65536 65536 98304 98304
add 53 2 0 0 131072 131072
mul 24 4 10 20 40 40'
BASE_TEXT='# tightness ratchet (S12): PROVISIONAL, see header of tools/tightness_baseline.txt
add 24 98304
add 53 131072
mul 24 40'
printf '%s\n' "$RESULTS_TEXT" > "$t/results.txt"
printf '%s\n' "$BASE_TEXT" > "$t/base.txt"

# 1. unchanged: pass
expect_rc 0 "unchanged results" "${tool[@]}" check "$t/results.txt" "$t/base.txt"
# 2. the report table
ck; out=$("${tool[@]}" report "$t/results.txt" 2>&1) || fail "report failed: $out"
ck; echo "$out" | grep -c '^add  *24  *3  *1.0000  *1.0000  *1.5000  *1.5000  *2  *ok' > /dev/null || fail "report lacks the add/24 row: $out"
ck; echo "$out" | grep -c '^add  *53  *2  *0.0000  *0.0000  *2.0000  *2.0000  *2  *ok' > /dev/null || fail "report lacks the add/53 row: $out"
ck; echo "$out" | grep -c '^mul  *24  *4 ' > /dev/null || fail "report lacks the mul/24 row: $out"
# 3. a worse maximum fails (one q16 unit above the baseline)
sed 's/^mul 24 4 10 20 40 40$/mul 24 4 10 20 30 41/' "$t/results.txt" > "$t/worse.txt"
expect_rc 1 "worse maximum" "${tool[@]}" check "$t/worse.txt" "$t/base.txt"
ck; echo "$out" > /dev/null; { "${tool[@]}" check "$t/worse.txt" "$t/base.txt" 2>&1 || true; } | grep -c 'WORSE' > /dev/null || fail "worse maximum not named WORSE"
# 4. an improved maximum passes, with a note
sed 's/^mul 24 4 10 20 40 40$/mul 24 4 10 20 30 39/' "$t/results.txt" > "$t/better.txt"
expect_rc 0 "improved maximum" "${tool[@]}" check "$t/better.txt" "$t/base.txt"
ck; { "${tool[@]}" check "$t/better.txt" "$t/base.txt" 2>&1 || true; } | grep -c 'improved' > /dev/null || fail "improvement not noted"
# 5. an unrecorded cell fails
{ cat "$t/results.txt"; echo "sqrt 24 5 1 2 3 4"; } > "$t/extra.txt"
expect_rc 1 "unrecorded cell" "${tool[@]}" check "$t/extra.txt" "$t/base.txt"
ck; { "${tool[@]}" check "$t/extra.txt" "$t/base.txt" 2>&1 || true; } | grep -c 'UNRECORDED' > /dev/null || fail "unrecorded cell not named"
# 6. a cell that disappeared fails
grep -v '^mul 24 ' "$t/results.txt" > "$t/gone.txt"
expect_rc 1 "missing cell" "${tool[@]}" check "$t/gone.txt" "$t/base.txt"
ck; { "${tool[@]}" check "$t/gone.txt" "$t/base.txt" 2>&1 || true; } | grep -c 'MISSING' > /dev/null || fail "missing cell not named"
# 7. no results file and an empty baseline pass; no results file and a non-empty baseline fail
printf '# empty\n' > "$t/empty.txt"
expect_rc 0 "absent results, empty baseline" "${tool[@]}" check "$t/nonexistent.txt" "$t/empty.txt"
expect_rc 1 "absent results, recorded baseline" "${tool[@]}" check "$t/nonexistent.txt" "$t/base.txt"
# 8. malformed results are an error (exit 2), never a pass
for bad in 'add 24 3 65536 65536 98304' 'add 24 3 70000 65536 98304 98304' 'add 24 0 0 0 0 0' 'add 024 3 1 1 1 1' 'add 24 3 1 1 1 1x' 'mul 24 4 10 20 40 40
add 24 3 65536 65536 98304 98304' 'add 24 3 1 1 1 1
add 24 3 1 1 1 1'; do
  { head -2 "$t/results.txt"; printf '%s\n' "$bad"; } > "$t/bad.txt"
  expect_rc 2 "malformed results: $bad" "${tool[@]}" check "$t/bad.txt" "$t/base.txt"
done
# 9. malformed baseline is an error (exit 2)
printf 'add 24\n' > "$t/badbase.txt"
expect_rc 2 "malformed baseline" "${tool[@]}" check "$t/results.txt" "$t/badbase.txt"
# 10. --update: lowers and adds, never raises, never drops
cp "$t/base.txt" "$t/b2.txt"
expect_rc 0 "update lowers" "${tool[@]}" update "$t/better.txt" "$t/b2.txt"
ck; grep -q '^mul 24 39$' "$t/b2.txt" || fail "update did not lower mul/24"
cp "$t/base.txt" "$t/b3.txt"
expect_rc 0 "update adds a cell" "${tool[@]}" update "$t/extra.txt" "$t/b3.txt"
ck; grep -q '^sqrt 24 4$' "$t/b3.txt" || fail "update did not add sqrt/24"
cp "$t/base.txt" "$t/b4.txt"
expect_rc 1 "update refuses to raise" "${tool[@]}" update "$t/worse.txt" "$t/b4.txt"
ck; cmp -s "$t/base.txt" "$t/b4.txt" || fail "a refused update changed the baseline"
cp "$t/base.txt" "$t/b5.txt"
expect_rc 1 "update refuses to drop a cell" "${tool[@]}" update "$t/gone.txt" "$t/b5.txt"
ck; cmp -s "$t/base.txt" "$t/b5.txt" || fail "a refused update changed the baseline (drop)"
# 11. the report marks a maximum above the provisional tau
{ head -2 "$t/results.txt"; echo 'add 24 3 1 1 1 131073'; } > "$t/over.txt"
ck; { "${tool[@]}" report "$t/over.txt" 2>&1 || true; } | grep -c 'OVER' > /dev/null || fail "report did not mark a maximum above tau"
# 12. the committed baseline and the wrapper agree with each other on an absent results file (no sets exist yet) or on real ones
"$here/tightness_check.sh" >/dev/null 2>&1; rc=$?
ck; [ "$rc" = 0 ] || fail "tools/tightness_check.sh fails on the committed tree (exit $rc)"
st_done tightness 33   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "tightness selftest: ok"
