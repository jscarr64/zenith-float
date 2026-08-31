#!/usr/bin/env bash
# Hex limb golds on aarch64 (static musl + qemu-user). Plan §20.2.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
target=aarch64-unknown-linux-musl
rustup target add "$target" >/dev/null
lld="$(rustc --print sysroot)/lib/rustlib/x86_64-unknown-linux-gnu/bin/rust-lld"
export CARGO_TARGET_AARCH64_UNKNOWN_LINUX_MUSL_LINKER="$lld"
export RUSTFLAGS="${RUSTFLAGS:-} -C target-feature=+crt-static"
cargo build -p zenith-float --release --bin hex_golds --target "$target"
bin="$(find "${CARGO_TARGET_DIR:-$root/target}/${target}/release" -name hex_golds -type f | head -1)"
[[ -n "$bin" ]] || { echo "error: aarch64 hex_golds not found" >&2; exit 1; }

qemu=""
if command -v qemu-aarch64 >/dev/null; then
  qemu="$(command -v qemu-aarch64)"
elif command -v qemu-aarch64-static >/dev/null; then
  qemu="$(command -v qemu-aarch64-static)"
elif [[ -x "$root/.tools/qemu-aarch64-static" ]]; then
  qemu="$root/.tools/qemu-aarch64-static"
fi
if [[ -z "$qemu" ]]; then
  echo "error: qemu-aarch64 or qemu-aarch64-static is required to run $target hex golds" >&2
  exit 1
fi
tmp="$(mktemp)"
trap 'rm -f "$tmp" "${tmp}.ref"' EXIT
"$qemu" "$bin" | sed '/^#/d' | grep -E '^[a-z]' >"$tmp"
sed '/^#/d' "$root/golds/hex/reference.txt" | grep -E '^[a-z]' >"${tmp}.ref"
if ! diff -u "${tmp}.ref" "$tmp"; then
  echo "error: aarch64 hex limb gold mismatch" >&2
  exit 1
fi
echo "ci_hex_arm: ok ($(wc -l < "${tmp}.ref") rows)"
