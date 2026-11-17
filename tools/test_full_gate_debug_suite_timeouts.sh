#!/usr/bin/env bash
# Red test (official f17d6ee gate, 2026-11-17): debug workspace suites in
# tools/full_gate.sh must default to 40m when STEP_TIMEOUT is unset.
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
[[ "$(timeout_for build-debug)" == 40m ]] || bad "build-debug default timeout is $(timeout_for build-debug), expected 40m"
[[ "$(timeout_for features-random-serde)" == 40m ]] || bad "features-random-serde default timeout is $(timeout_for features-random-serde), expected 40m"
[[ "$(timeout_for build-release)" == 20m ]] || bad "build-release default timeout is $(timeout_for build-release), expected 20m"
[[ "$(timeout_for oracle-fast)" == 20m ]] || bad "oracle-fast default timeout is $(timeout_for oracle-fast), expected 20m"
[[ "$(timeout_for mpfr)" == 40m ]] || bad "mpfr default timeout is $(timeout_for mpfr), expected 40m"

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
