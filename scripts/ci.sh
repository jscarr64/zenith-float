#!/usr/bin/env bash
# Default public test slew. Does not enable mpfr-tests (needs rug / MPFR).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
cargo test --workspace
cargo test -p zenith-float-num --lib --no-default-features --features std
cargo test --workspace --features random,serde
echo "zenith-float ci ok"
