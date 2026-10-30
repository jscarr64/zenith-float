#!/usr/bin/env bash
# P6 (64-bit minimum): zenith-float builds only for target_pointer_width = "64".
# Asserts that `cargo check` for each 32-bit target FAILS with the exact compile_error! message, and that
# 64-bit targets still check. i686-unknown-linux-musl is required (missing = FAIL); other targets run
# only when installed (reported as SKIP). Usage: scripts/check_64bit_only.sh   (CARGO_TARGET_DIR honoured)
set -u
MSG='zenith-float 2.0 requires a 64-bit target (target_pointer_width = 64); 32-bit targets are not supported.'
cd "$(dirname "$0")/.." || exit 2
installed=$(rustup target list --installed)
has() { grep -qx "$1" <<<"$installed"; }
pass=0 fail=0 skip=0
ok() { echo "PASS $*"; pass=$((pass + 1)); }
bad() { echo "FAIL $*"; fail=$((fail + 1)); }
refuse() { # target pkg [extra cargo args]
  local t=$1 p=$2; shift 2
  local out rc; out=$(cargo check -p "$p" --target "$t" "$@" 2>&1); rc=$?
  if [ $rc -eq 0 ]; then bad "$t $p: cargo check succeeded (must refuse 32-bit)"
  elif grep -qF "$MSG" <<<"$out"; then ok "$t $p: refused with the message (exit $rc)"
  else bad "$t $p: exit $rc without the message"; tail -n 8 <<<"$out"; fi
}
accept() { # target [extra cargo args]
  local t=$1; shift
  local out rc; out=$(cargo check --target "$t" "$@" 2>&1); rc=$?
  if [ $rc -eq 0 ]; then ok "$t: cargo check $* succeeded"; else bad "$t: cargo check $* exit $rc"; tail -n 8 <<<"$out"; fi
}
# 32-bit: must refuse.
if has i686-unknown-linux-musl; then
  for p in zenith-float-num zenith-float zenith-float-compare; do refuse i686-unknown-linux-musl "$p"; done
else bad "i686-unknown-linux-musl not installed (required for this check)"; fi
if has thumbv7em-none-eabihf; then refuse thumbv7em-none-eabihf zenith-float-num --no-default-features --lib
else echo "SKIP thumbv7em-none-eabihf: target not installed"; skip=$((skip + 1)); fi
# 64-bit: must build.
accept x86_64-unknown-linux-gnu --workspace
if has aarch64-unknown-linux-gnu; then accept aarch64-unknown-linux-gnu --workspace
else echo "SKIP aarch64-unknown-linux-gnu: target not installed"; skip=$((skip + 1)); fi
echo "check_64bit_only: $pass passed, $fail failed, $skip skipped"
[ $fail -eq 0 ]
