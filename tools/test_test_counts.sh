#!/usr/bin/env bash
# Red test (item 3, host): every host test step of tools/full_gate.sh (build-debug, build-release,
# features-random-serde, oracle-fast, oracle-deep, mpfr) FAILs on a binary that ran 0 tests, was filtered to 0,
# had all tests ignored, or lost its summary; features-random-serde also FAILs when a test target with
# required-features did not run (feature-skipped). Stub cargo on PATH, throwaway POPLOG_ROOT.
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY CARGO_BUILD_JOBS RUST_TEST_THREADS \
  GATE_STEPS
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$("$root/tools/selftest_scratch.sh" test_test_counts) || exit 2; trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin"
cat > "$tmp/bin/cargo" <<'SH'
#!/usr/bin/env bash
echo "$*" >> "${0%/bin/cargo}/cargo.args"
if [[ "$1" == metadata ]]; then
  t='{"name":"zf","kind":["lib"],"required-features":[]},{"name":"integ_a","kind":["test"],"required-features":[]}'
  t+=',{"name":"needs_random","kind":["test"],"required-features":["random"]}'
  printf '{"packages":[{"name":"zf","features":{"random":[],"default":[]},"targets":[%s]}]}\n' "$t"; exit 0
fi
if [[ "$1" == test && "$*" == *--message-format=json* ]]; then
  b="$CARGO_TARGET_DIR/debug/deps/stub-0"; mkdir -p "${b%/*}"; : > "$b"
  printf '{"reason":"compiler-artifact","executable":"%s","profile":{"test":true}}\n' "$b"; exit 0
fi
if [[ "$1" == test ]]; then
  d="$CARGO_TARGET_DIR/debug/deps"; mkdir -p "$d"; : > "$d/integ_a-0"; : > "$d/needs_random-0"
  echo "     Running tests/integ_a.rs ($d/integ_a-0)"
  case "${STUB_HOST_LOSS:-}" in
    zero) echo "test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out";;
    filtered) echo "test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 7 filtered out";;
    ignored) echo "test result: ok. 0 passed; 0 failed; 4 ignored; 0 measured";;
    nosummary) :;;
    *) echo "test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured";;
  esac
  if [[ "$*" == *random* && "${STUB_HOST_LOSS:-}" != featskip ]]; then
    echo "     Running tests/needs_random.rs ($d/needs_random-0)"; echo "test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured"
  fi
fi
exit 0
SH
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"
chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
fail=0
gate() {
  local c=$1; shift; rc=0; : > "$tmp/cargo.args"
  out=$(env "$@" HOME="$tmp" PATH="$tmp/bin:$PATH" POPLOG_ROOT="$tmp/pop-$c" T="$tmp/gate-$c" \
    bash "$root/tools/full_gate.sh" "host-$c" 2>&1) || rc=$?
  log=$(find "$tmp/pop-$c" -name '*.log' 2>/dev/null | sort | sed -n 1p || true)
  [[ -n "$log" ]] || { echo "FAIL[$c]: no poplog"; fail=1; log=/dev/null; }
}
HOST="build-debug build-release features-random-serde oracle-fast oracle-deep mpfr"
for s in $HOST; do
  gate "green-$s" GATE_STEPS=$s
  grep -q "STEP END $s verdict=PASS" "$log" && [[ $rc -eq 0 ]] || { echo "FAIL[green-$s]: green host step did not PASS (rc=$rc)"; fail=1; }
  for loss in zero filtered ignored nosummary; do
    gate "$loss-$s" STUB_HOST_LOSS=$loss GATE_STEPS=$s
    grep -q "STEP END $s verdict=FAIL" "$log" && [[ $rc -ne 0 ]] || { echo "FAIL[$loss-$s]: coverage loss did not FAIL $s (rc=$rc)"; fail=1; }
    grep -q "$s coverage: FAIL.*integ_a" "$tmp/gate-$loss-$s/host-$loss-$s-$s.log" || { echo "FAIL[$loss-$s]: step log does not name integ_a"; fail=1; }
  done
done
gate featskip STUB_HOST_LOSS=featskip GATE_STEPS=features-random-serde
grep -q 'STEP END features-random-serde verdict=FAIL' "$log" && [[ $rc -ne 0 ]] || { echo "FAIL[featskip]: a skipped required-features test did not FAIL (rc=$rc)"; fail=1; }
grep -q 'coverage: FAIL.*needs_random(feature-skipped)' "$tmp/gate-featskip/host-featskip-features-random-serde.log" || { echo "FAIL[featskip]: step log does not name needs_random(feature-skipped)"; fail=1; }
# script level: an unverifiable zero-test entry is exit 2, not a pass
printf '     Running tests/a.rs (/x/a-0)\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured\n' > "$tmp/one.out"
mkdir -p "$tmp/emptyroot"
rc=0; TEST_COUNTS_ROOT="$tmp/emptyroot" bash "$root/tools/host_test_counts.sh" x "$tmp/one.out" >/dev/null 2>&1 || rc=$?
[[ $rc -eq 2 ]] || { echo "FAIL[zerolist]: unverifiable zero-test list gave rc=$rc, want 2"; fail=1; }
bash "$root/tools/host_test_counts.sh" x "$tmp/one.out" >/dev/null 2>&1 || { echo "FAIL[zerolist]: one green binary rejected"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: host test steps FAIL on 0-run, filtered, all-ignored, no-summary and feature-skipped binaries"
exit $fail
