#!/usr/bin/env bash
# Red test (Jeff 2026-11-23 23:04 CT: "debug still runs every test for correctness"): no test may branch on
# debug_assertions to skip rows or tests, and no per-row wait may fail a debug row for being slow (raw recv_timeout /
# wait_timeout; Spence R1 on 82e4afb) (tools/debug_skips_scan.py). Allowed: the shared timing switch in
# zenith-float-testkit/src/test_time.rs and the hang-watchdog budgets. A fixture check proves the scan bites.
# usage: tools/test_no_debug_skips.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
scan="$root/tools/debug_skips_scan.py"
fail=0
bad() { echo "FAIL: $*"; fail=1; }

zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/test-no-debug-skips"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp (Spence n3)
fx=$(mktemp -d "$zf_tmp/no-debug-skips.XXXXXX") || exit 2
trap 'rm -rf "$fx"' EXIT
git -C "$fx" init -q
mkdir -p "$fx/tests"
printf '%s\n' '#[test]' 'fn watchdog() {' '    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {' \
  '        std::time::Duration::from_secs(300)' '    } else {' '        std::time::Duration::from_secs(20)' '    };' \
  '    let b = BUDGET.saturating_mul(if cfg!(debug_assertions) { 10 } else { 1 });' \
  '    // cfg!(debug_assertions) in a comment is fine' '    let r = test_time::wait_row(&rx, BUDGET);' '}' > "$fx/tests/ok.rs"
python3 "$scan" "$fx" > /dev/null || bad "scan flags the allowed watchdog fixture"
for skip in 'if cfg!(debug_assertions) { return; }' '#[cfg(not(debug_assertions))]' \
  "#[cfg_attr(debug_assertions, ignore)]" 'let p_max = if cfg!(debug_assertions) { 5 } else { 8 };' \
  'let outcome = match rx.recv_timeout(ROW_BUDGET) {' 'let got = rx.recv_timeout(std::time::Duration::from_secs(60));'; do
  printf '%s\n' "$skip" > "$fx/tests/skip.rs"
  python3 "$scan" "$fx" > /dev/null && bad "scan misses: $skip"
done

if ! out=$(python3 "$scan" "$root"); then
  bad "tests run less in debug, or fail a debug row for being slow:"
  printf '%s\n' "$out"
fi
[[ $fail -eq 0 ]] && { echo "test_no_debug_skips: ok"; exit 0; }
exit 1
