#!/usr/bin/env bash
# Red test: tools/full_gate.sh runs the i686 lib tests as a blocking step `test-i686-lib`
# (cargo test -p zenith-float-num --lib --target i686-unknown-linux-musl, own target subdir), in the default
# step list right after build-release; GATE_STEPS selects steps by exact name; a red i686 run FAILs the gate
# at that step (non-zero exit, step named in the RESULT line). Stub cargo/rustup on PATH, throwaway POPLOG_ROOT.
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY CARGO_BUILD_JOBS RUST_TEST_THREADS \
  GATE_STEPS
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin"
cat > "$tmp/bin/cargo" <<'SH'
#!/usr/bin/env bash
echo "$*" >> "${0%/bin/cargo}/cargo.args"
case "$*" in audit\ --version) echo "cargo-audit 0.0"; exit 0;; esac
if [[ "$1" == test && "$*" == *--message-format=json* ]]; then
  b="$CARGO_TARGET_DIR/debug/deps/stub-0"; mkdir -p "${b%/*}"; : > "$b"
  printf '{"reason":"compiler-artifact","executable":"%s","profile":{"test":true}}\n' "$b"; exit 0
fi
if [[ "$1" == test && "$*" == *i686* && "${STUB_I686_RED:-0}" == 1 ]]; then
  echo "test result: FAILED. 0 passed; 1 failed; 0 ignored; 0 measured"; exit 101
fi
if [[ "$1" == test ]]; then echo "test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured"; fi
exit 0
SH
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"
chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
fail=0
# gate CASE [VAR=VALUE...]: run the stubbed gate in its own POPLOG_ROOT/T; sets rc, log, out
gate() {
  local c=$1; shift; rc=0; : > "$tmp/cargo.args"
  out=$(env "$@" HOME="$tmp" PATH="$tmp/bin:$PATH" POPLOG_ROOT="$tmp/pop-$c" T="$tmp/gate-$c" \
    bash "$root/tools/full_gate.sh" "i686-$c" 2>&1) || rc=$?
  log=$(find "$tmp/pop-$c" -name '*.log' 2>/dev/null | head -1 || true)
  [[ -n "$log" ]] || { echo "FAIL[$c]: no poplog"; fail=1; log=/dev/null; }
}
# 1. default list: the step exists, right after build-release, red i686 -> step FAIL, gate FAIL naming it
gate default STUB_I686_RED=1
grep -A3 'STEP END build-release' "$log" | grep -q 'STEP START test-i686-lib' || { echo "FAIL[default]: no test-i686-lib step right after build-release"; fail=1; }
grep -q 'STEP END test-i686-lib verdict=FAIL' "$log" || { echo "FAIL[default]: red i686 run is not a FAIL step"; fail=1; }
grep -q 'test -p zenith-float-num --lib --target i686-unknown-linux-musl --no-run --message-format=json' "$tmp/cargo.args" || { echo "FAIL[default]: no cargo test --lib --target i686-unknown-linux-musl via fresh_test_run.sh"; fail=1; }
[[ -f "$tmp/target/i686/.worktree-owner" ]] || { echo "FAIL[default]: test-i686-lib did not use its own target subdir \$CARGO_TARGET_DIR/i686"; fail=1; }
[[ $rc -ne 0 ]] || { echo "FAIL[default]: gate exited 0"; fail=1; }
# 2. GATE_STEPS=test-i686-lib, red: only that step runs, the gate FAILs at it (blocking)
gate red STUB_I686_RED=1 GATE_STEPS=test-i686-lib
[[ $(grep -c 'STEP START' "$log") -eq 1 ]] && grep -q 'STEP START test-i686-lib' "$log" || { echo "FAIL[red]: GATE_STEPS=test-i686-lib did not run exactly that step: $(grep 'STEP START' "$log" | tr '\n' ' ')"; fail=1; }
grep -q 'RESULT: FAIL' "$log" || { echo "FAIL[red]: footer not FAIL"; fail=1; }
grep -q 'RESULT: FAIL .*test-i686-lib' <<< "$out" || { echo "FAIL[red]: RESULT line does not name test-i686-lib: $(grep RESULT <<< "$out")"; fail=1; }
[[ $rc -ne 0 ]] || { echo "FAIL[red]: red i686 step did not fail the gate (exit 0): not blocking"; fail=1; }
# 3. GATE_STEPS=test-i686-lib, green: PASS, exit 0
gate green GATE_STEPS=test-i686-lib
grep -q 'STEP END test-i686-lib verdict=PASS' "$log" && grep -q 'RESULT: PASS' "$log" && [[ $rc -eq 0 ]] || { echo "FAIL[green]: green i686 run is not PASS (rc=$rc)"; fail=1; }
# 4. exact names: a prefix selects nothing and an unknown name fails the gate
gate prefix GATE_STEPS=test-i686
! grep -q 'STEP START' "$log" && [[ $rc -ne 0 ]] || { echo "FAIL[prefix]: GATE_STEPS=test-i686 ran a step or exited 0 (rc=$rc)"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: full_gate.sh test-i686-lib (default list, blocking, exact GATE_STEPS)"
exit "$fail"
