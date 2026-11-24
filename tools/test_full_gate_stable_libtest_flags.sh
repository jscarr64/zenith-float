#!/usr/bin/env bash
# Red test (official 545cb3c gate, 2026-11-17): every libtest flag in
# tools/full_gate.sh TT must be accepted by this toolchain. --fail-fast is
# nightly-only; build-debug exited 101 immediately with:
#   error: The "fail-fast" flag is only accepted on the nightly compiler
#   with -Z unstable-options
# usage: tools/test_full_gate_stable_libtest_flags.sh
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
line="$(grep -E '^TT=' "$root/tools/full_gate.sh")" || {
  echo "FAIL: no TT= line in tools/full_gate.sh"
  exit 1
}
# Expand $RUST_TEST_THREADS the same way full_gate.sh does.
RUST_TEST_THREADS="${RUST_TEST_THREADS:-2}"
eval "$line"
# TT is of the form: -- --test-threads=N [more libtest flags]
if [[ "$TT" != --\ --* ]]; then
  echo "FAIL: TT= is not in the expected '-- --flags' form: $TT"
  exit 1
fi
# shellcheck disable=SC2086
set -- $TT
# $1 is the cargo/libtest separator "--"
shift
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/test-full-gate-stable-libtest-flags"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp (Spence n3)
t="$(mktemp -d "$zf_tmp/run.XXXXXX")"
trap 'rm -rf "$t"' EXIT
printf '#[test] fn t() {}\n' > "$t/t.rs"
rustc --test "$t/t.rs" -o "$t/tbin"
if ! out="$("$t/tbin" "$@" 2>&1)"; then
  echo "FAIL: full_gate.sh TT libtest flags are rejected by $(rustc --version | awk '{print $1, $2}'): $*"
  printf '%s\n' "$out"
  exit 1
fi
echo "test_full_gate_stable_libtest_flags: ok"
