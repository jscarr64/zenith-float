#!/usr/bin/env bash
# Diff example stdout (TSV rows) against golds/hex/reference.txt.
# Extra args are forwarded to `cargo run` (e.g. --target i686-unknown-linux-gnu).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
ref="$root/golds/hex/reference.txt"
[[ -f "$ref" ]] || { echo "error: missing $ref" >&2; exit 1; }
tmp="$(mktemp)"
trap 'rm -f "$tmp" "${tmp}.ref"' EXIT
cargo run -p zenith-float --release --bin hex_golds "$@" | sed '/^#/d' | grep -E '^[a-z]' >"$tmp"
sed '/^#/d' "$ref" | grep -E '^[a-z]' >"${tmp}.ref"
if ! diff -u "${tmp}.ref" "$tmp"; then
  echo "error: hex limb gold mismatch" >&2
  exit 1
fi
echo "hex limb golds match $(basename "$ref") ($(wc -l < "${tmp}.ref") rows)"
