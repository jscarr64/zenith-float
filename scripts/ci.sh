#!/usr/bin/env bash
# Default public test slew. MPFR oracle tests require: cargo test -p zenith-float-num --features mpfr-tests
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if rg 'f32|f64' --glob '*.rs' --glob '*.md' --glob 'CHANGELOG*' . --glob '!doc/BUILD_CHECKLIST.md'; then
  echo "error: f32/f64 identifiers are forbidden in zenith-float"
  exit 1
fi

cargo test --workspace
cargo test -p zenith-float-num --features random --test radix_roundtrip
cargo test -p zenith-float-num --lib --release
cargo test -p zenith-float-num --lib --no-default-features --features std
cargo test --workspace --features random,serde

# MPFR bit-oracle gate (Linux x86_64 only; release, single-threaded).
if [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]]; then
  cargo test -p zenith-float-num --features mpfr-tests --release -- --test-threads=1
fi

# Optional nightly: Criterion vs doc/bench-baselines.tsv (slow).
if [[ "${CI_BENCH:-}" == 1 ]]; then
  bash "$root/scripts/bench-compare.sh"
fi

echo "zenith-float ci ok"
