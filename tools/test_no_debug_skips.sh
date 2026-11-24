#!/usr/bin/env bash
# Red test (Jeff 2026-11-23 23:04 CT: "debug still runs every test for correctness"): no test may branch on
# debug_assertions to skip rows or tests (tools/debug_skips_scan.py). Allowed: the shared timing switch in
# zenith-float-num/src/test_time.rs and the hang-watchdog budgets. A fixture check proves the scan bites.
# usage: tools/test_no_debug_skips.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
scan="$root/tools/debug_skips_scan.py"
fail=0
bad() { echo "FAIL: $*"; fail=1; }

fx=$(mktemp -d "${TMPDIR:-/tmp}/no-debug-skips.XXXXXX") || exit 2
trap 'rm -rf "$fx"' EXIT
git -C "$fx" init -q
mkdir -p "$fx/tests"
printf '%s\n' '#[test]' 'fn watchdog() {' '    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {' \
  '        std::time::Duration::from_secs(300)' '    } else {' '        std::time::Duration::from_secs(20)' '    };' \
  '    let b = BUDGET.saturating_mul(if cfg!(debug_assertions) { 10 } else { 1 });' \
  '    // cfg!(debug_assertions) in a comment is fine' '}' > "$fx/tests/ok.rs"
python3 "$scan" "$fx" > /dev/null || bad "scan flags the allowed watchdog fixture"
for skip in 'if cfg!(debug_assertions) { return; }' '#[cfg(not(debug_assertions))]' \
  "#[cfg_attr(debug_assertions, ignore)]" 'let p_max = if cfg!(debug_assertions) { 5 } else { 8 };'; do
  printf '%s\n' "$skip" > "$fx/tests/skip.rs"
  python3 "$scan" "$fx" > /dev/null && bad "scan misses: $skip"
done

if ! out=$(python3 "$scan" "$root"); then
  bad "tests branch on debug_assertions to run less in debug:"
  printf '%s\n' "$out"
fi
[[ $fail -eq 0 ]] && { echo "test_no_debug_skips: ok"; exit 0; }
exit 1
