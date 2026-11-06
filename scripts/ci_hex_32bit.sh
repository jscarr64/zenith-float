#!/usr/bin/env bash
# Hex limb golds on 32-bit (`WORD_BIT_SIZE = 32`). Plan §20.2.
# Static musl so the host does not need gcc-multilib. Interchange bytes match x86_64
# (`to_bytes` uses big-endian u32 limbs).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"
target=i686-unknown-linux-musl
rustup target add --toolchain 1.99 "$target" >/dev/null
export RUSTFLAGS="${RUSTFLAGS:-} -C target-feature=+crt-static"
bash "$root/scripts/ci_hex_common.sh" --target "$target"
# Clippy on the 32-bit word size: lints that only fire when WORD_BIT_SIZE = 32.
cargo clippy -p zenith-float-num --all-targets --target "$target" -- -D warnings
echo "ci_hex_32bit: ok"
