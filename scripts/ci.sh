#!/usr/bin/env bash
# Default public test slew. MPFR oracle tests require: cargo test -p zenith-float-num --features mpfr-tests
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if rg 'f32|f64' --glob '*.rs' --glob '*.md' --glob 'CHANGELOG*' .; then
  echo "error: f32/f64 identifiers are forbidden in zenith-float"
  exit 1
fi

cargo test --workspace
cargo test -p zenith-float-num --lib --release
cargo test -p zenith-float-num --lib --no-default-features --features std
cargo test --workspace --features random,serde
echo "zenith-float ci ok"
