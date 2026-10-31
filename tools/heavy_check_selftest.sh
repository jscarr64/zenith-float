#!/usr/bin/env bash
# Self-test of the heavy-tier check (ZB-P03, Paul 2026-10-30): tools/heavy_check.sh and the full_gate steps matrix-heavy and skip-scan.
# A heavy log PASSES only with exactly one `HEAVY s6_heavy_rows_replay_through_zb: rows run=30 passed=30` line and no SKIP line; a SKIP
# line in any gate step log is NOT-RUN (FAIL); the matrix-heavy step FAILS with ZB_MATRIX_HEAVY unset (never a skip) or set to a file
# that does not match HEAVY.SHA256. Fixtures only, under $(zb_scratch); no repo file is touched.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap heavy_check_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
export LC_ALL=C
root=$(cd "$(dirname "$0")/.." && pwd)
fail() { echo "heavy_check selftest FAIL: $*"; exit 1; }
t=$(mktemp -d "$(zb_scratch)"/zb-heavy-check-selftest.XXXXXX) || { echo "mktemp failed"; exit 2; }; runlog_on_exit "rm -rf \"$t\""
. "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
hc() { "$root/tools/heavy_check.sh" "$@" > "$t/hc.out" 2>&1; }
ok30='HEAVY s6_heavy_rows_replay_through_zb: rows run=30 passed=30 file=/x/matrix_rows_heavy.txt'
skip='SKIP s6_heavy_rows_replay_through_zb: ZB_MATRIX_HEAVY unset, 0 of 30 heavy rows run'
res='test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 361 filtered out; finished in 1.68s'
printf 'running 1 test\n%s\ntest slice_matrix::s6_heavy_rows_replay_through_zb ... ok\n%s\n' "$ok30" "$res" > "$t/good.log"
printf 'running 1 test\n%s\ntest slice_matrix::s6_heavy_rows_replay_through_zb ... ok\n%s\n' "$skip" "$res" > "$t/skip.log"
sed 's/run=30 passed=30/run=29 passed=29/' "$t/good.log" > "$t/rows29.log"
sed 's/run=30 passed=30/run=30 passed=29/' "$t/good.log" > "$t/pass29.log"
{ cat "$t/good.log"; printf '%s\n' "$skip"; } > "$t/goodskip.log"
{ cat "$t/good.log"; printf '%s\n' "$ok30"; } > "$t/twice.log"
# --heavy: one heavy-run log
ck; hc --heavy "$t/good.log" || fail "real heavy log (run=30 passed=30) failed: $(cat "$t/hc.out")"
ck; hc --heavy "$t/skip.log" && fail "heavy log with a SKIP line passed"
ck; grep -q 'NOT-RUN' "$t/hc.out" || fail "SKIP not reported as NOT-RUN: $(cat "$t/hc.out")"
ck; hc --heavy "$t/rows29.log" && fail "heavy log with run=29 passed"
ck; hc --heavy "$t/pass29.log" && fail "heavy log with passed=29 passed"
ck; hc --heavy "$t/goodskip.log" && fail "heavy log with run=30 AND a SKIP line passed"
ck; hc --heavy "$t/twice.log" && fail "heavy log with two HEAVY lines passed"
ck; hc --heavy "$t/missing.log" && fail "missing heavy log passed"
# --logs: every step log of one gate tag
mkdir -p "$t/g"; cp "$t/good.log" "$t/g/x-build-release.log"; printf 'clean\n' > "$t/g/x-fmt.log"
ck; hc --logs "$t/g" x skip-scan || fail "clean gate logs failed: $(cat "$t/hc.out")"
printf 'step output\n%s\n' "$skip" > "$t/g/x-build-debug.log"
ck; hc --logs "$t/g" x skip-scan && fail "a gate log with a SKIP line passed"
ck; hc --logs "$t/g" none skip-scan && fail "no gate logs at all passed (vacuous)"
# --pin: the heavy file must match HEAVY.SHA256
printf 'rows\n' > "$t/h.txt"; printf '# pin\n%s  matrix_rows_heavy.txt\n' "$(sha256sum "$t/h.txt" | awk '{ print $1 }')" > "$t/PIN"
ck; hc --pin "$t/PIN" "$t/h.txt" || fail "matching heavy file failed the pin: $(cat "$t/hc.out")"
printf 'other\n' > "$t/o.txt"
ck; hc --pin "$t/PIN" "$t/o.txt" && fail "unpinned heavy file passed the pin"
# full_gate: matrix-heavy never skips, skip-scan fails a SKIP line in any step log
fe=(env -u GATE_STEPS -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE RUNLOG_POPLOG=/nonexistent RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1)
rcof() { sed -n "s/^$2=//p" "$t/logs/$1.rc"; }
mkdir -p "$t/repo/tests/goldens/slices/matrix" "$t/logs"; cp "$t/PIN" "$t/repo/tests/goldens/slices/matrix/HEAVY.SHA256"
ck; "${fe[@]}" -u ZB_MATRIX_HEAVY GATE_STEPS=matrix-heavy T="$t/logs" REPO="$t/repo" "$root/tools/full_gate.sh" unset > /dev/null 2>&1 && fail "gate passed matrix-heavy with ZB_MATRIX_HEAVY unset"
ck; [ "$(rcof unset matrix-heavy)" != "" ] && [ "$(rcof unset matrix-heavy)" != 0 ] || fail "matrix-heavy did not run and fail with ZB_MATRIX_HEAVY unset (rc '$(rcof unset matrix-heavy)')"
ck; "${fe[@]}" ZB_MATRIX_HEAVY="$t/o.txt" GATE_STEPS=matrix-heavy T="$t/logs" REPO="$t/repo" "$root/tools/full_gate.sh" badpin > /dev/null 2>&1 && fail "gate passed matrix-heavy with an unpinned heavy file"
printf 'step output\n%s\n' "$skip" > "$t/logs/planted-build-debug.log"
ck; "${fe[@]}" GATE_STEPS=skip-scan T="$t/logs" REPO="$t/repo" "$root/tools/full_gate.sh" planted > /dev/null 2>&1 && fail "gate passed with a SKIP line in a step log"
ck; [ "$(rcof planted skip-scan)" = 1 ] || fail "skip-scan did not fail on the planted SKIP line (rc '$(rcof planted skip-scan)')"
st_done heavy_check 19
echo "heavy_check selftest: ok (19 cases)"
