#!/usr/bin/env bash
# Hex limb golds on 32-bit (`WORD_BIT_SIZE = 32`). Plan §20.2.
# Static musl so the host does not need gcc-multilib. Interchange bytes match x86_64
# (`to_bytes` uses big-endian u32 limbs).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
target=i686-unknown-linux-musl
rustup target add "$target" >/dev/null
export RUSTFLAGS="${RUSTFLAGS:-} -C target-feature=+crt-static"
bash "$root/scripts/ci_hex_common.sh" --target "$target"
echo "ci_hex_32bit: ok"
