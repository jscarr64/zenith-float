#!/usr/bin/env bash
# Red test: tools/full_gate.sh runs the i686 tests as a blocking step `test-i686` (cargo test of the whole
# workspace minus the oracle tests, --target i686-unknown-linux-musl, own target subdir), in the default
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
if [[ "$1" == metadata ]]; then
  # workspace fixture: a lib, one plain integration test, one behind a feature, and the five oracle tests
  t='{"name":"zf","kind":["lib"],"required-features":[]}'
  for n in integ_a oracle specials_oracle methods_oracle accuracy_table elliptic_phi; do t+=",{\"name\":\"$n\",\"kind\":[\"test\"],\"required-features\":[]}"; done
  t+=',{"name":"needs_random","kind":["test"],"required-features":["random"]}'
  printf '{"packages":[{"name":"zf","targets":[%s]}]}\n' "$t"; exit 0
fi
if [[ "$1" == test && "$*" == *i686* && "$*" != *--message-format=json* && "${STUB_I686_INTEG_RED:-0}" == 1 ]]; then
  d="$CARGO_TARGET_DIR/debug/deps"; mkdir -p "$d"; : > "$d/zf-0"; : > "$d/integ_a-0"
  echo "     Running unittests src/lib.rs ($d/zf-0)"; echo "test result: ok. 5 passed; 0 failed; 1 ignored; 0 measured"
  echo "     Running tests/integ_a.rs ($d/integ_a-0)"; echo "test result: FAILED. 3 passed; 1 failed; 0 ignored; 0 measured"; exit 101
fi
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
grep -A3 'STEP END build-release' "$log" | grep -q 'STEP START test-i686' || { echo "FAIL[default]: no test-i686 step right after build-release"; fail=1; }
grep -q 'STEP END test-i686 verdict=FAIL' "$log" || { echo "FAIL[default]: red i686 run is not a FAIL step"; fail=1; }
grep -q 'test --workspace --lib --bins --test integ_a --target i686-unknown-linux-musl --no-run --message-format=json' "$tmp/cargo.args" || { echo "FAIL[default]: no cargo test --workspace (minus oracles) --target i686-unknown-linux-musl via fresh_test_run.sh"; fail=1; }
[[ -f "$tmp/target/i686/.worktree-owner" ]] || { echo "FAIL[default]: test-i686 did not use its own target subdir \$CARGO_TARGET_DIR/i686"; fail=1; }
[[ $rc -ne 0 ]] || { echo "FAIL[default]: gate exited 0"; fail=1; }
# 2. GATE_STEPS=test-i686, red: only that step runs, the gate FAILs at it (blocking)
gate red STUB_I686_RED=1 GATE_STEPS=test-i686
[[ $(grep -c 'STEP START' "$log") -eq 1 ]] && grep -q 'STEP START test-i686' "$log" || { echo "FAIL[red]: GATE_STEPS=test-i686 did not run exactly that step: $(grep 'STEP START' "$log" | tr '\n' ' ')"; fail=1; }
grep -q 'RESULT: FAIL' "$log" || { echo "FAIL[red]: footer not FAIL"; fail=1; }
grep -q 'RESULT: PARTIAL FAIL .*test-i686' <<< "$out" || { echo "FAIL[red]: RESULT line does not name test-i686: $(grep RESULT <<< "$out")"; fail=1; }
[[ $rc -ne 0 ]] || { echo "FAIL[red]: red i686 step did not fail the gate (exit 0): not blocking"; fail=1; }
# 3. GATE_STEPS=test-i686, green: PASS, exit 0
gate green GATE_STEPS=test-i686
grep -q 'STEP END test-i686 verdict=PASS' "$log" && grep -q 'RESULT: PASS' "$log" && [[ $rc -eq 0 ]] || { echo "FAIL[green]: green i686 run is not PASS (rc=$rc)"; fail=1; }
# 4. exact names: a prefix selects nothing and an unknown name fails the gate
gate prefix GATE_STEPS=test-i68
! grep -q 'STEP START' "$log" && [[ $rc -ne 0 ]] || { echo "FAIL[prefix]: GATE_STEPS=test-i68 ran a step or exited 0 (rc=$rc)"; fail=1; }
# 5. a selective run is PARTIAL in stdout, the .rc and the poplog footer; a default run is plain
gate partial GATE_STEPS=test-i686
grep -q '^RESULT: PARTIAL PASS' <<< "$out" && grep -q '^RESULT: PARTIAL PASS' "$tmp/gate-partial/i686-partial.rc" \
  && sed -n '/POPLOG FOOTER/,$p' "$log" | grep -q 'RESULT: PARTIAL PASS' || { echo "FAIL[partial]: green GATE_STEPS run is not PARTIAL PASS in stdout/.rc/footer: $(grep RESULT <<< "$out")"; fail=1; }
gate partialred STUB_I686_RED=1 GATE_STEPS=test-i686
grep -q '^RESULT: PARTIAL FAIL .*test-i686' <<< "$out" && sed -n '/POPLOG FOOTER/,$p' "$log" | grep -q 'RESULT: PARTIAL FAIL' || { echo "FAIL[partialred]: red GATE_STEPS run is not PARTIAL FAIL naming the step"; fail=1; }
gate plain STUB_I686_RED=1
grep -q '^RESULT: FAIL ' <<< "$out" && ! grep -q PARTIAL <<< "$out" || { echo "FAIL[plain]: default run is not a plain RESULT: FAIL"; fail=1; }
# 6. GATE_STEPS set but matching no step ("," or " "): FAIL with a clear message, never PASS
for sel in "," " "; do
  gate empty GATE_STEPS="$sel"
  [[ $rc -ne 0 ]] && grep -q '^RESULT: PARTIAL FAIL .*no step matched' <<< "$out" || { echo "FAIL[empty]: GATE_STEPS='$sel' ran nothing but rc=$rc: $(grep RESULT <<< "$out")"; fail=1; }
  rm -rf "$tmp/pop-empty" "$tmp/gate-empty"
done
# 7. the whole workspace minus the oracle tests: an integration-test (non-lib) failure on i686 FAILs the gate,
#    and the gate log has per-binary counts for the step
gate integ STUB_I686_INTEG_RED=1 GATE_STEPS=test-i686
grep -q 'STEP END test-i686 verdict=FAIL' "$log" && [[ $rc -ne 0 ]] || { echo "FAIL[integ]: an i686 integration-test failure did not FAIL the gate (rc=$rc)"; fail=1; }
grep -q 'i686 binary integ_a: passed=3 failed=1 ignored=0' "$log" && grep -q 'i686 binary zf: passed=5 failed=0 ignored=1' "$log" || { echo "FAIL[integ]: no per-binary i686 counts in the gate log"; fail=1; }
a=$(grep -m1 'i686-unknown-linux-musl --no-run' "$tmp/cargo.args" || true)
[[ " $a " == *" --test integ_a "* && " $a " != *" --test oracle "* && " $a " != *" --test needs_random "* ]] || { echo "FAIL[integ]: i686 args do not include integ_a and exclude the oracle tests by name: $a"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: full_gate.sh test-i686 (default list, blocking, exact GATE_STEPS)"
exit "$fail"
