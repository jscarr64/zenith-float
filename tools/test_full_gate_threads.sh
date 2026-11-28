#!/usr/bin/env bash
# Red test (P14, Jeff 2026-11-10 21:15 CT "The computer is all your tonight"; Paul: no test step single-threaded
# unless required): tools/full_gate.sh, scripts/zenith_prepublish.sh and scripts/ci.sh never pin a test step to one
# thread, and full_gate.sh's own fallback for CARGO_BUILD_JOBS / RUST_TEST_THREADS (no popslot slot value) is the
# host's core count, not a fixed number. Measured: the mpfr-tests suite passes at 32 threads, 649 passed, 3:08 wall
# (~/pop-logs/rusty/2026-11-10/fold-pr1-69-mpfr-32threads.log), against 19:06 at one thread (gate3 step mpfr).
# Each check also runs on a planted copy and must catch it.
# usage: tools/test_full_gate_threads.sh    (no cargo; text checks only)
set -u
root=$(cd "$(dirname "$0")/.." && pwd)
fail=0; bad() { echo "FAIL: $*"; fail=1; }
check() {  # DIR -> 0 when no step is pinned to one thread and the fallback is the core count
  local d=$1 r=0 hits f
  for f in tools/full_gate.sh scripts/zenith_prepublish.sh scripts/ci.sh; do
    hits=$(grep -nE -- '--test-threads[= ]+1([^0-9]|$)|RUST_TEST_THREADS=1([^0-9]|$)' "$d/$f")
    [ -z "$hits" ] || { echo "  single-threaded test step in $f: $hits"; r=1; }
  done
  hits=$(grep -nE '(CARGO_BUILD_JOBS|RUST_TEST_THREADS)=\$\{(CARGO_BUILD_JOBS|RUST_TEST_THREADS):-\$\{POPSLOT_JOBS:-[0-9]+\}\}' "$d/tools/full_gate.sh")
  [ -z "$hits" ] || { echo "  fixed job fallback in tools/full_gate.sh: $hits"; r=1; }
  grep -qE 'RUST_TEST_THREADS=\$\{RUST_TEST_THREADS:-\$\{POPSLOT_JOBS:-\$\(nproc\)\}\}' "$d/tools/full_gate.sh" \
    || { echo "  tools/full_gate.sh: RUST_TEST_THREADS does not fall back to \$(nproc)"; r=1; }
  return $r
}
check "$root" || bad "a gate test step is pinned to one thread, or the job fallback is not the core count"
t=$(mktemp -d "${TMPDIR:-$HOME/.cache}/gatethreads.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/tools" "$t/scripts"
plant() {  # FILE LINE
  cp "$root/tools/full_gate.sh" "$t/tools/"; cp "$root/scripts/zenith_prepublish.sh" "$root/scripts/ci.sh" "$t/scripts/"
  printf '%s\n' "$2" >> "$t/$1"
  if check "$t" > /dev/null; then bad "planted '$2' in $1 was not caught"; fi
}
plant tools/full_gate.sh 'run mpfr "tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests -- --test-threads=1"'
plant tools/full_gate.sh 'export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-${POPSLOT_JOBS:-4}}'
plant scripts/zenith_prepublish.sh '"$root/tools/fresh_test_run.sh" cargo test -p zenith-float-num --test proptest_props -- --test-threads 1'
plant scripts/ci.sh '  RUST_TEST_THREADS=1 cargo test -p zenith-float-num --features mpfr-tests --release'
[ $fail = 0 ] && { echo "test_full_gate_threads: ok"; exit 0; }
exit 1
