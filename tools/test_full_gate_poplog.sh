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
case "$*" in audit\ --version) echo "cargo-audit 0.0"; exit 0;; esac
if [[ "$1" == test ]]; then echo "test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured"; fi
exit 0
SH
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"
chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export POPSLOT_PREFIX="$tmp/slot" POPLOG_ROOT="$tmp/pop" T="$tmp/gate" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
# HOME points at the tmp dir so the vendored scripts/poplog.sh fallback (used by CI) is exercised.
rc=0; HOME="$tmp" PATH="$tmp/bin:$PATH" bash "$root/tools/full_gate.sh" redtest > /dev/null 2>&1 || rc=$?
log=$(find "$tmp/pop" -name '*.log' 2>/dev/null | head -1 || true)
fail=0
[[ -n "$log" ]] || { echo "FAIL: full_gate.sh wrote no poplog"; exit 1; }
grep -q '^commit:' "$log" || { echo "FAIL: no header"; fail=1; }
grep -q 'STEP END build-debug verdict=FAIL .*run=0 passed=0' "$log" || { echo "FAIL: zero-test build-debug not FAIL"; fail=1; }
grep -q 'SLOT TAKEN slot=1/4 .*CARGO_BUILD_JOBS=6 RUST_TEST_THREADS=6' "$log" || { echo "FAIL: gate did not take a popslot slot at 6 jobs"; fail=1; }
grep -q 'RESULT: FAIL' "$log" || { echo "FAIL: footer not FAIL"; fail=1; }
[[ $rc -ne 0 ]] || { echo "FAIL: gate exited 0 with zero tests"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: full_gate.sh poplog (zero tests -> FAIL, log $log)"
exit "$fail"
