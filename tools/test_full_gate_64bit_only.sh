#!/usr/bin/env bash
# Red test (Jeff 2026-11-10 19:15 CT, "Yes 64 bit only please"; P6): the release gate and its helpers build, lint and test
# only on the 64-bit host. tools/full_gate.sh has no step for a 32-bit target, scripts/zenith_prepublish.sh and scripts/ci.sh
# build no 32-bit target, and the step list is exactly the 64-bit list below (a step added or dropped must update it).
# Each check also runs on a planted copy (a 32-bit step put back) and must catch it.
# ACCEPTANCE-pr39 F1/F5/F6 (Spence 2026-11-10): every library crate root (facade, num, bound) refuses a non-64-bit
# target with compile_error!, no script under scripts/ builds or runs a 32-bit or wasm target (comments excepted),
# and nostd-tests/.cargo/ configures no 32-bit target.
# usage: tools/test_full_gate_64bit_only.sh    (no cargo; text checks only)
set -u
root=$(cd "$(dirname "$0")/.." && pwd)
PAT='thumbv[0-9]|i686|i586|armv7|arm-unknown|wasm32|hex-32bit|ci_hex_32bit|--target'
PAT32='thumbv[0-9]|i686|i586|armv7|arm-unknown|wasm|hex-32bit|ci_hex_32bit'  # targets, not the --target flag
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
  for f in src/lib.rs zenith-float-num/src/lib.rs zenith-float-bound/src/lib.rs; do
    grep -A1 -E '^#\[cfg\(not\(target_pointer_width = "64"\)\)\]$' "$d/$f" | grep -qE '^compile_error!\(".*64-bit' \
      || { echo "  $f: no 64-bit compile_error!"; r=1; }
  done
  hits=$(grep -nE "$PAT32" "$d"/scripts/*.sh 2>/dev/null | grep -vE '^[^:]+:[0-9]+:[[:space:]]*#')
  [ -z "$hits" ] || { echo "  32-bit/wasm target in scripts/: $hits"; r=1; }
  hits=$(cat "$d"/nostd-tests/.cargo/config* 2>/dev/null | grep -E "$PAT32")
  [ -z "$hits" ] || { echo "  32-bit target in nostd-tests/.cargo: $hits"; r=1; }
  return $r
}
check "$root" || bad "the gate is not 64-bit only"
t=$(mktemp -d "${TMPDIR:-$HOME/.cache}/gate64.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/tools" "$t/scripts" "$t/src" "$t/zenith-float-num/src" "$t/zenith-float-bound/src" "$t/nostd-tests/.cargo"
fresh() {  # copy every checked file into the planted tree
  rm -rf "$t/scripts" "$t/nostd-tests/.cargo"; mkdir -p "$t/scripts" "$t/nostd-tests/.cargo"
  cp "$root/tools/full_gate.sh" "$t/tools/"; cp "$root"/scripts/*.sh "$t/scripts/"
  for f in src/lib.rs zenith-float-num/src/lib.rs zenith-float-bound/src/lib.rs; do cp "$root/$f" "$t/$f"; done
  cp "$root"/nostd-tests/.cargo/config* "$t/nostd-tests/.cargo/" 2>/dev/null || true
}
plant() {  # FILE LINE
  fresh; printf '%s\n' "$2" >> "$t/$1"
  if check "$t" > /dev/null; then bad "planted '$2' in $1 was not caught"; fi
}
unplant() {  # FILE : drop the 64-bit compile_error! pair from a crate root
  fresh; grep -vE '^(#\[cfg\(not\(target_pointer_width = "64"\)\)\]|compile_error!\(".*64-bit.*\);)$' "$root/$1" > "$t/$1"
  if check "$t" > /dev/null; then bad "removed compile_error! in $1 was not caught"; fi
}
plant tools/full_gate.sh 'run clippy-i686 "cargo clippy -p zenith-float-num --target i686-unknown-linux-musl"'
plant scripts/zenith_prepublish.sh 'cargo build --no-default-features --target thumbv7em-none-eabihf'
plant scripts/ci.sh '  cargo build --no-default-features --target thumbv7em-none-eabihf'
plant scripts/ci_hex_32bit.sh 'target=i686-unknown-linux-musl'
plant scripts/ci_nostd.sh 'cargo build --release --features qemu --bin nostd-qemu --target thumbv7em-none-eabihf'
plant scripts/ci_hex_wasm.sh 'target=wasm32-wasip1'
plant nostd-tests/.cargo/config.toml '[target.thumbv7em-none-eabihf]'
for f in src/lib.rs zenith-float-num/src/lib.rs zenith-float-bound/src/lib.rs; do unplant "$f"; done
[ $fail = 0 ] && { echo "test_full_gate_64bit_only: ok"; exit 0; }
exit 1
