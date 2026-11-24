#!/usr/bin/env bash
# Red test (official f17d6ee gate, 2026-11-17): debug workspace suites in
# tools/full_gate.sh must default to 40m when STEP_TIMEOUT is unset (build-debug: 60m since
# Jeff 2026-11-23 23:04 CT; features-random-serde stays 40m).
# features-random-serde runs `cargo test --workspace --features random,serde`
# in debug (no --release), the same heavy suite as build-debug. Official
# full_gate at f17d6ee was fail-fast killed at features-random-serde=124
# (20m cap) with tests still passing, mid s3_euler_gamma.
# usage: tools/test_full_gate_debug_suite_timeouts.sh
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
fail=0
bad() { echo "FAIL: $*"; fail=1; }

# Evaluate the timeout-selection body from run() in full_gate.sh for a step
# name, so this tracks the gate itself rather than a stale copy.
timeout_for() {
  local name=$1
  local to body
  body=$(awk '
    /^run\(\) \{/ { in_run=1 }
    in_run && /to=\$\{STEP_TIMEOUT/ { grab=1 }
    grab {
      if ($0 ~ /timeout "\$to"/) exit
      print
    }
  ' "$root/tools/full_gate.sh")
  if [[ -z "$body" ]]; then
    echo "FAIL: could not extract timeout selection from run() in tools/full_gate.sh"
    exit 1
  fi
  eval "$body"
  printf '%s\n' "$to"
}

# features-random-serde must stay a debug workspace suite: that is why it
# shares build-debug's 40m default. build-release on this SHA finished in
# 4m 40s under the 20m cap.
feat="$(grep -E '^run features-random-serde ' "$root/tools/full_gate.sh" || true)"
[[ -n "$feat" ]] || bad "no run features-random-serde line in tools/full_gate.sh"
if [[ -n "$feat" ]]; then
  [[ "$feat" == *'cargo test --workspace'* ]] || bad "features-random-serde is not a workspace cargo test: $feat"
  [[ "$feat" == *--release* ]] && bad "features-random-serde is a --release suite; the 40m default is for the debug workspace tests: $feat"
fi

# Defaults when STEP_TIMEOUT / MPFR_STEP_TIMEOUT are unset.
unset STEP_TIMEOUT || true
unset MPFR_STEP_TIMEOUT || true
# Jeff 2026-11-23 23:04 CT: "timing limits are enforced in the release step only; debug still runs every
# test for correctness, and the debug step limit goes to 60 min". build-debug is the only limit that moved.
[[ "$(timeout_for build-debug)" == 60m ]] || bad "build-debug default timeout is $(timeout_for build-debug), expected 60m (Jeff 2026-11-23 23:04 CT)"
[[ "$(timeout_for features-random-serde)" == 40m ]] || bad "features-random-serde default timeout is $(timeout_for features-random-serde), expected 40m"
[[ "$(timeout_for build-release)" == 20m ]] || bad "build-release default timeout is $(timeout_for build-release), expected 20m"
[[ "$(timeout_for oracle-fast)" == 20m ]] || bad "oracle-fast default timeout is $(timeout_for oracle-fast), expected 20m"
[[ "$(timeout_for mpfr)" == 40m ]] || bad "mpfr default timeout is $(timeout_for mpfr), expected 40m"

[[ "$(timeout_for oracle-deep)" == 20m ]] || bad "oracle-deep default timeout is $(timeout_for oracle-deep), expected 20m"
[[ "$(timeout_for prepublish)" == 20m ]] || bad "prepublish default timeout is $(timeout_for prepublish), expected 20m"
[[ "$(timeout_for clippy-x86_64)" == 20m ]] || bad "clippy-x86_64 default timeout is $(timeout_for clippy-x86_64), expected 20m"

# The gate sets no case wall: ZFS_WALL (120 s in the T2 runs) belongs to the T2 harness, not to full_gate.sh.
grep -q 'ZFS_WALL' "$root/tools/full_gate.sh" && bad "tools/full_gate.sh sets or mentions ZFS_WALL; the case wall is the T2 harness's"

# Every wall-clock budget in the tests is pinned (file, native value, count) in tools/time_limits.tsv; the
# release-only rule changes where limits apply, never their values. tools/time_limits_scan.sh lists them.
if ! diff <(bash "$root/tools/time_limits_scan.sh") "$root/tools/time_limits.tsv" > /dev/null 2>&1; then
  bad "wall-clock budgets differ from tools/time_limits.tsv: $(diff <(bash "$root/tools/time_limits_scan.sh") "$root/tools/time_limits.tsv" 2>&1 | head -5 | tr '\n' ' ')"
fi

# Jeff 2026-11-23 23:04 CT: debug still runs every test; no test may branch on debug_assertions to run less.
bash "$root/tools/test_no_debug_skips.sh" > /dev/null 2>&1 || bad "tools/test_no_debug_skips.sh: $(bash "$root/tools/test_no_debug_skips.sh" 2>&1 | head -8 | tr '\n' ' ')"

# Explicit STEP_TIMEOUT must still win, including for the debug suites
# (selftests set STEP_TIMEOUT=60s).
STEP_TIMEOUT=60s
[[ "$(timeout_for build-debug)" == 60s ]] || bad "explicit STEP_TIMEOUT=60s not honored for build-debug: $(timeout_for build-debug)"
[[ "$(timeout_for features-random-serde)" == 60s ]] || bad "explicit STEP_TIMEOUT=60s not honored for features-random-serde: $(timeout_for features-random-serde)"
[[ "$(timeout_for build-release)" == 60s ]] || bad "explicit STEP_TIMEOUT=60s not honored for build-release: $(timeout_for build-release)"
unset STEP_TIMEOUT || true

if [[ $fail -eq 0 ]]; then
  echo "test_full_gate_debug_suite_timeouts: ok"
  exit 0
fi
exit 1
