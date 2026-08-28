#!/usr/bin/env bash
# Default public test slew. MPFR oracle tests require: cargo test -p zenith-float-num --features mpfr-tests
#
# Wall-time budgets (override with CI_DEBUG_SECS / CI_MPFR_SECS). Skip checks with CI_SKIP_TIME_BUDGET=1.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

DEBUG_BUDGET="${CI_DEBUG_SECS:-600}"
MPFR_BUDGET="${CI_MPFR_SECS:-1800}"

run_budget() {
  local name="$1"
  local budget="$2"
  shift 2
  local start=$SECONDS
  "$@"
  local elapsed=$((SECONDS - start))
  echo "ci: ${name} ${elapsed}s (budget ${budget}s)"
  if [[ "${CI_SKIP_TIME_BUDGET:-}" == 1 ]]; then
    return 0
  fi
  if (( elapsed > budget )); then
    echo "error: ${name} exceeded wall-time budget (${elapsed}s > ${budget}s)" >&2
    exit 1
  fi
}

debug_slew() {
  if rg 'f32|f64' --glob '*.rs' --glob '*.md' --glob 'CHANGELOG*' . --glob '!doc/BUILD_CHECKLIST.md'; then
    echo "error: f32/f64 identifiers are forbidden in zenith-float"
    exit 1
  fi

  cargo test --workspace
  cargo test -p zenith-float-num --features random --test radix_roundtrip
  cargo test -p zenith-float-num --lib --release
  cargo test -p zenith-float-num --lib --no-default-features --features std
  cargo test --workspace --features random,serde
}

mpfr_slew() {
  cargo test -p zenith-float-num --features mpfr-tests --release -- --test-threads=1
}

run_budget debug "${DEBUG_BUDGET}" debug_slew

# MPFR bit-oracle gate (Linux x86_64 only; release, single-threaded).
if [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]]; then
  run_budget mpfr "${MPFR_BUDGET}" mpfr_slew
fi

# Optional nightly: Criterion vs doc/bench-baselines.tsv (slow).
if [[ "${CI_BENCH:-}" == 1 ]]; then
  bash "$root/scripts/bench-compare.sh"
fi

echo "zenith-float ci ok"
