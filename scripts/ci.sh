#!/usr/bin/env bash
# Default public test slew. MPFR oracle tests require: cargo test -p zenith-float-num --features mpfr-tests
#
# Wall-time budgets (override with CI_DEBUG_SECS / CI_MPFR_SECS). Skip checks with CI_SKIP_TIME_BUDGET=1.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"

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
  # Hardware float tokens, arch intrinsics, inline asm, and kernel unsafe.
  bash "$root/scripts/ci_no_intrinsics.sh"
  # CHANGELOG.md states the 2.0.0 MSRV in the words the release note requires.
  if rg -n -i 'accumath' --glob '!scripts/ci.sh' --glob '!CHANGELOG.md' .; then
    echo "error: this public crate must not mention that name"
    exit 1
  fi

  cargo test --workspace
  cargo test -p zenith-float-num --features random --test radix_roundtrip
  cargo test -p zenith-float-num --lib --release
  cargo test -p zenith-float-num --lib --no-default-features --features std
  cargo test --workspace --features random,serde
  cargo build -p zenith-float-num --no-default-features
  rustup target add --toolchain 1.99 thumbv7em-none-eabihf >/dev/null
  cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf
  cargo build --no-default-features --target thumbv7em-none-eabihf
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
