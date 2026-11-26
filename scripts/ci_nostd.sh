#!/usr/bin/env bash
# no_std gate for zenith-float-num (64-bit only: P6/P8; Jeff 2026-11-10 "Yes 64 bit only please").
#
# 1. The library builds with --no-default-features for bare-metal 64-bit targets (no std at all):
#    x86_64-unknown-none and aarch64-unknown-none.
# 2. nostd-tests (a #![no_std] crate: core + alloc only) runs the mpmath reference cases on the host as a
#    #![no_std]/#![no_main] binary linking only libc. It exits non-zero if any case misses its bit target.
# The former thumbv7em-none-eabihf (32-bit) builds and its QEMU runner were removed with the 32-bit targets.
#
# Needs: rustup targets x86_64-unknown-none and aarch64-unknown-none.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"
for t in x86_64-unknown-none aarch64-unknown-none; do
  rustup target add --toolchain 1.99 "$t" >/dev/null
  echo "== cargo build -p zenith-float-num --no-default-features --lib --target $t"
  cargo build -p zenith-float-num --no-default-features --lib --target "$t"
  echo "== cargo build -p zenith-float --no-default-features --lib --target $t"
  cargo build -p zenith-float --no-default-features --lib --target "$t"
done

cd nostd-tests
echo "== nostd-host (x86_64-unknown-linux-gnu, #![no_std] + libc)"
cargo build --release --features host --bin nostd-host
"${CARGO_TARGET_DIR:-target}/release/nostd-host"
