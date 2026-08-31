#!/usr/bin/env bash
# Hex limb golds on wasm32-wasip1 via wasmtime. Plan §20.2.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
target=wasm32-wasip1
rustup target add "$target" >/dev/null

ensure_wasmtime() {
  if command -v wasmtime >/dev/null; then
    command -v wasmtime
    return
  fi
  local ver="v24.0.2"
  local dir="$root/.tools/wasmtime-${ver}"
  local bin="$dir/wasmtime"
  if [[ ! -x "$bin" ]]; then
    mkdir -p "$dir"
    local url="https://github.com/bytecodealliance/wasmtime/releases/download/${ver}/wasmtime-${ver}-x86_64-linux.tar.xz"
    curl -fsSL "$url" | tar -xJ -C "$dir" --strip-components=1
  fi
  echo "$bin"
}

wt="$(ensure_wasmtime)"
cargo build -p zenith-float --release --bin hex_golds --target "$target"
wasm="$(find "${CARGO_TARGET_DIR:-$root/target}/${target}/release" -name 'hex_golds.wasm' | head -1)"
[[ -n "$wasm" ]] || { echo "error: hex_golds.wasm not found" >&2; exit 1; }
tmp="$(mktemp)"
trap 'rm -f "$tmp" "${tmp}.ref"' EXIT
"$wt" run "$wasm" | sed '/^#/d' | grep -E '^[a-z]' >"$tmp"
sed '/^#/d' "$root/golds/hex/reference.txt" | grep -E '^[a-z]' >"${tmp}.ref"
if ! diff -u "${tmp}.ref" "$tmp"; then
  echo "error: wasm32 hex limb gold mismatch" >&2
  exit 1
fi
echo "ci_hex_wasm: ok ($(wc -l < "${tmp}.ref") rows)"
