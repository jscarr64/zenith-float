#!/usr/bin/env bash
# Red test for Jeff's logging rule (2026-10-23): tools/full_gate.sh must write a poplog with a header,
# per-step counts, a footer, and must FAIL a test step whose cargo test reports zero tests.
# Uses stub cargo/rustup on PATH and a throwaway POPLOG_ROOT; never touches ~/pop-logs.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin"
cat > "$tmp/bin/cargo" <<'SH'
#!/usr/bin/env bash
echo "$*" >> "${0%/bin/cargo}/cargo.args"
case "$*" in audit\ --version) echo "cargo-audit 0.0"; exit 0;; esac
# fresh_test_run.sh builds with --no-run --message-format=json first: report one test binary inside CARGO_TARGET_DIR
if [[ "$1" == test && "$*" == *--message-format=json* ]]; then
  b="$CARGO_TARGET_DIR/debug/deps/stub-0"; mkdir -p "${b%/*}"; : > "$b"
  printf '{"reason":"compiler-artifact","executable":"%s","profile":{"test":true}}\n' "$b"; exit 0
fi
if [[ "$1" == test ]]; then echo "test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured"; fi
exit 0
SH
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"
chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" POPLOG_ROOT="$tmp/pop" T="$tmp/gate" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
# HOME points at the tmp dir so the vendored scripts/poplog.sh fallback (used by CI) is exercised.
rc=0; HOME="$tmp" PATH="$tmp/bin:$PATH" bash "$root/tools/full_gate.sh" redtest > /dev/null 2>&1 || rc=$?
log=$(find "$tmp/pop" -name '*.log' 2>/dev/null | head -1 || true)
fail=0
[[ -n "$log" ]] || { echo "FAIL: full_gate.sh wrote no poplog"; exit 1; }
grep -q '^commit:' "$log" || { echo "FAIL: no header"; fail=1; }
grep -q 'STEP END build-debug verdict=FAIL .*run=0 passed=0' "$log" || { echo "FAIL: zero-test build-debug not FAIL"; fail=1; }
# slot geometry: POPSLOT_SLOTS x POPSLOT_JOBS, default 8 x 4 (Jeff 2026-10-23: 32 cores); read from the env, not assumed
ns=${POPSLOT_SLOTS:-8} nj=${POPSLOT_JOBS:-4}
grep -q "SLOT TAKEN slot=1/$ns .*CARGO_BUILD_JOBS=$nj RUST_TEST_THREADS=$nj" "$log" || { echo "FAIL: gate did not take a popslot slot 1/$ns at $nj jobs: $(grep -m1 'SLOT TAKEN' "$log")"; fail=1; }
# early mpfr-tests compile (Paul 2026-10-23, 6ab75ce prepublish): right after build-debug, cargo test --no-run with the feature
grep -A3 'STEP END build-debug' "$log" | grep -q 'STEP START mpfr-tests-build' || { echo "FAIL: no mpfr-tests-build step right after build-debug"; fail=1; }
# through tools/fresh_test_run.sh (Paul, one target dir per worktree): the --no-run build with JSON output, and no test run after it
grep -qx 'test -p zenith-float-num --release --features mpfr-tests --no-run --message-format=json-render-diagnostics' "$tmp/cargo.args" || { echo "FAIL: mpfr-tests-build did not run cargo test --no-run --features mpfr-tests via fresh_test_run.sh"; fail=1; }
# target-dir (its selftest needs a real cargo, so it fails under the stub) runs before build-debug
grep -m1 -o 'STEP START [a-z0-9-]*' "$log" | grep -qx 'STEP START target-dir' || { echo "FAIL: target-dir is not the first step"; fail=1; }
grep -q 'RESULT: FAIL' "$log" || { echo "FAIL: footer not FAIL"; fail=1; }
[[ $rc -ne 0 ]] || { echo "FAIL: gate exited 0 with zero tests"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: full_gate.sh poplog (zero tests -> FAIL, log $log)"
exit "$fail"
