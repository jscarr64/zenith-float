#!/usr/bin/env bash
# Diff example stdout (TSV rows) against golds/hex/reference.txt.
# Extra args are forwarded to `cargo run` (e.g. --target aarch64-unknown-linux-musl).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"
ref="$root/golds/hex/reference.txt"
[[ -f "$ref" ]] || { echo "error: missing $ref" >&2; exit 1; }
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/ci-hex-common"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp (Spence n3)
tmp="$(mktemp "$zf_tmp/f.XXXXXX")"
trap 'rm -f "$tmp" "${tmp}.ref"' EXIT
cargo run -q -p zenith-float --release --bin hex_golds "$@" | sed '/^#/d' | grep -E '^[a-z]' >"$tmp"
sed '/^#/d' "$ref" | grep -E '^[a-z]' >"${tmp}.ref"
if ! diff -u "${tmp}.ref" "$tmp"; then
  echo "error: hex limb gold mismatch" >&2
  exit 1
fi
echo "hex limb golds match $(basename "$ref") ($(wc -l < "${tmp}.ref") rows)"
