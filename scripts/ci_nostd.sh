#!/usr/bin/env bash
# no_std gate for zenith-float-num.
#
# 1. The library builds with --no-default-features (no std; core + alloc) on the host.
# 2. nostd-tests (a #![no_std] crate: core + alloc only) runs the mpmath reference cases on the host
#    as a #![no_std]/#![no_main] binary linking only libc. It exits non-zero if any case misses its bit target.
# The bare-metal thumbv7em (32-bit) QEMU runner was retired with P6 (64-bit targets only, 2026-10-30).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"

echo "== cargo build -p zenith-float-num --no-default-features --lib"
cargo build -p zenith-float-num --no-default-features --lib
echo "== cargo build -p zenith-float --no-default-features --lib"
cargo build -p zenith-float --no-default-features --lib

cd nostd-tests
echo "== nostd-host (x86_64-unknown-linux-gnu, #![no_std] + libc)"
cargo build --release --features host --bin nostd-host
"${CARGO_TARGET_DIR:-target}/release/nostd-host"
