#!/usr/bin/env bash
# no_std gate for zenith-float-num.
#
# 1. The library builds with --no-default-features for a bare-metal target (no std at all).
# 2. nostd-tests (a #![no_std] crate: core + alloc only) runs the mpmath reference cases
#    - on the host as a #![no_std]/#![no_main] binary linking only libc, and
#    - on thumbv7em-none-eabihf (Cortex-M4F, 32-bit limbs) under qemu-system-arm mps2-an386
#      with semihosting.
# Each runner exits non-zero if any case misses its bit target.
#
# Needs: rustup target thumbv7em-none-eabihf; qemu-system-arm on PATH (or QEMU=/path/to/it).
set -euo pipefail
cd "$(dirname "$0")/.."
QEMU="${QEMU:-qemu-system-arm}"

echo "== cargo build -p zenith-float-num --no-default-features --lib --target thumbv7em-none-eabihf"
cargo build -p zenith-float-num --no-default-features --lib --target thumbv7em-none-eabihf
echo "== cargo build -p zenith-float --no-default-features --lib --target thumbv7em-none-eabihf"
cargo build -p zenith-float --no-default-features --lib --target thumbv7em-none-eabihf

cd nostd-tests
echo "== nostd-host (x86_64-unknown-linux-gnu, #![no_std] + libc)"
cargo build --release --features host --bin nostd-host
"${CARGO_TARGET_DIR:-target}/release/nostd-host"

echo "== nostd-qemu (thumbv7em-none-eabihf, Cortex-M4F under QEMU)"
cargo build --release --features qemu --bin nostd-qemu --target thumbv7em-none-eabihf
"$QEMU" -machine mps2-an386 -cpu cortex-m4 -nographic -semihosting-config enable=on,target=native \
  -kernel "${CARGO_TARGET_DIR:-target}/thumbv7em-none-eabihf/release/nostd-qemu"
