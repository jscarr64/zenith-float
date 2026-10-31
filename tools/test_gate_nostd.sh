#!/usr/bin/env bash
# Red test (ZF-T03): the full gate runs scripts/ci_nostd.sh (nostd-host build + run) as step ci-nostd,
# and a failing ci_nostd.sh FAILs the gate (stub cargo that exits 1 => step verdict FAIL, gate rc != 0).
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY CARGO_BUILD_JOBS RUST_TEST_THREADS \
  GATE_STEPS
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$("$root/tools/selftest_scratch.sh" test_gate_nostd) || exit 2; trap 'rm -rf "$tmp"' EXIT
fail=0
grep -Fqx 'run ci-nostd "bash scripts/ci_nostd.sh"' "$root/tools/full_gate.sh" || { echo "FAIL: no ci-nostd step running scripts/ci_nostd.sh"; fail=1; }
mkdir -p "$tmp/bin"
printf '#!/usr/bin/env bash\nexit 1\n' > "$tmp/bin/cargo"
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"; chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
rc=0
env GATE_STEPS=ci-nostd HOME="$tmp" PATH="$tmp/bin:$PATH" POPLOG_ROOT="$tmp/pop" T="$tmp/gate" \
  bash "$root/tools/full_gate.sh" t03-ci-nostd >/dev/null 2>&1 || rc=$?
log=$(find "$tmp/pop" -name '*.log' 2>/dev/null | sort | sed -n 1p || true)
if [[ -n "$log" ]] && grep -q "STEP END ci-nostd verdict=FAIL" "$log" && [[ $rc -ne 0 ]]; then :; else
  echo "FAIL: broken ci_nostd.sh did not FAIL the gate as step ci-nostd (rc=$rc)"; fail=1; fi
[[ $fail -eq 0 ]] || exit 1
echo "PASS test_gate_nostd"
