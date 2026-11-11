#!/usr/bin/env bash
# Red test (Jeff 2026-11-10 19:15 CT, "Yes 64 bit only please"; P6): the release gate and its helpers build, lint and test
# only on the 64-bit host. tools/full_gate.sh has no step for a 32-bit target, scripts/zenith_prepublish.sh and scripts/ci.sh
# build no 32-bit target, and the step list is exactly the 64-bit list below (a step added or dropped must update it).
# Each check also runs on a planted copy (a 32-bit step put back) and must catch it.
# usage: tools/test_full_gate_64bit_only.sh    (no cargo; text checks only)
set -u
root=$(cd "$(dirname "$0")/.." && pwd)
PAT='thumbv[0-9]|i686|i586|armv7|arm-unknown|wasm32|hex-32bit|ci_hex_32bit|--target'
WANT="target-dir golden-sums arm-tests-release qemu-time-factor build-debug mpfr-tests-build zero-test build-release build-nostd features-random-serde oracle-fast oracle-deep mpfr rules-scan no-clone-progress fmt clippy-x86_64 audit no-intrinsics unsafe-scan prepublish"
fail=0; bad() { echo "FAIL: $*"; fail=1; }
check() {  # DIR -> 0 when 64-bit only
  local d=$1 r=0 hits steps
  hits=$(grep -nE "^run [a-z0-9_-]+ \".*($PAT)" "$d/tools/full_gate.sh")
  [ -z "$hits" ] || { echo "  32-bit gate step: $hits"; r=1; }
  steps=$(grep -oE '^run [a-z0-9_-]+' "$d/tools/full_gate.sh" | awk '{print $2}' | tr '\n' ' ' | sed 's/ $//')
  [ "$steps" = "$WANT" ] || { echo "  gate steps: '$steps'"; r=1; }
  for f in scripts/zenith_prepublish.sh scripts/ci.sh; do
    hits=$(grep -nE "^[[:space:]]*(cargo|rustup)[^#]*($PAT)" "$d/$f")
    [ -z "$hits" ] || { echo "  32-bit build in $f: $hits"; r=1; }
  done
  return $r
}
check "$root" || bad "the gate is not 64-bit only"
t=$(mktemp -d "${TMPDIR:-$HOME/.cache}/gate64.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/tools" "$t/scripts"
plant() {  # FILE LINE
  cp "$root/tools/full_gate.sh" "$t/tools/"; cp "$root/scripts/zenith_prepublish.sh" "$root/scripts/ci.sh" "$t/scripts/"
  printf '%s\n' "$2" >> "$t/$1"
  if check "$t" > /dev/null; then bad "planted '$2' in $1 was not caught"; fi
}
plant tools/full_gate.sh 'run clippy-i686 "cargo clippy -p zenith-float-num --target i686-unknown-linux-musl"'
plant scripts/zenith_prepublish.sh 'cargo build --no-default-features --target thumbv7em-none-eabihf'
plant scripts/ci.sh '  cargo build --no-default-features --target thumbv7em-none-eabihf'
[ $fail = 0 ] && { echo "test_full_gate_64bit_only: ok"; exit 0; }
exit 1
