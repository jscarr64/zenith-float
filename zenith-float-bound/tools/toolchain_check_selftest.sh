#!/usr/bin/env bash
# Self-test of tools/toolchain_check.sh with fake rustc/cargo on PATH and fixture rust-toolchain.toml files: the pinned channel passes,
# any other version (older, a longer prefix such as 1.990), a missing or unpinned toolchain file, or a missing rustc fails.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap toolchain_check_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
here=$(cd "$(dirname "$0")" && pwd)
fail() { echo "toolchain_check selftest FAIL: $*"; exit 1; }
t=$(mktemp -d) || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
mkfake() { # mkfake <dir> <rustc version line> <cargo version line>
  mkdir -p "$1"
  # shellcheck disable=SC2016 # shim body expands when the shim runs, not here
  printf '#!/bin/sh\necho "%s"\n' "$2" > "$1/rustc"; printf '#!/bin/sh\necho "%s"\n' "$3" > "$1/cargo"; chmod +x "$1/rustc" "$1/cargo"
}
printf '[toolchain]\nchannel = "1.99"\ncomponents = ["clippy"]\n' > "$t/pinned.toml"
printf '[toolchain]\nchannel = "stable"\n' > "$t/stable.toml"
mkfake "$t/ok" "rustc 1.99.1 (58cf8f9c3 2026-10-13)" "cargo 1.99.1 (5f94df478 2026-08-27)"
mkfake "$t/old" "rustc 1.95.0 (59807616e 2026-04-14)" "cargo 1.99.1 (5f94df478 2026-08-27)"
mkfake "$t/oldcargo" "rustc 1.99.1 (58cf8f9c3 2026-10-13)" "cargo 1.95.0 (aaaaaaaaa 2026-04-14)"
mkfake "$t/long" "rustc 1.990.0 (0 2027-01-01)" "cargo 1.990.0 (0 2027-01-01)"
sys="$t/sys"; mkdir -p "$sys"; for c in bash sed; do ln -s "$(command -v "$c")" "$sys/$c"; done
# Fixture runs of the real toolchain_check.sh log through poplog into a throwaway POPLOG_ROOT, never the live ~/pop-logs
# (RUNNING.tsv rows, per-run logs). Only the fixture runs get it; this selftest's own run log stays live.
plr=$(mktemp -d) || fail "mktemp"; runlog_on_exit "rm -rf \"$plr\""
chk() { POPLOG_ROOT="$plr" PATH="$1:$sys" "$here/toolchain_check.sh" "$2" > "$t/out" 2>&1; }
ck; chk "$t/ok" "$t/pinned.toml" || fail "pinned 1.99.1 rejected: $(cat "$t/out")"
ck; grep -q 'rustc 1.99.1' "$t/out" || fail "version not printed"
ck; chk "$t/old" "$t/pinned.toml" && fail "rustc 1.95.0 accepted for channel 1.99"
ck; grep -q 'TOOLCHAIN MISMATCH' "$t/out" || fail "mismatch not named"
ck; chk "$t/oldcargo" "$t/pinned.toml" && fail "cargo 1.95.0 accepted for channel 1.99"
ck; chk "$t/long" "$t/pinned.toml" && fail "1.990.0 accepted for channel 1.99"
ck; chk "$t/ok" "$t/missing.toml" && fail "missing toolchain file accepted"
ck; chk "$t/ok" "$t/stable.toml" && fail "unpinned channel accepted"
ck; mkdir -p "$t/none"; chk "$t/none" "$t/pinned.toml" && fail "missing rustc accepted"
# (Under RUNLOG_PLAIN=1, e.g. full_gate_selftest's fixture gate, nothing logs anywhere, so there is nothing to check.)
ck; [ "${RUNLOG_PLAIN:-0}" = 1 ] || [ -n "$(find "$plr" -name '*-toolchain_check*.log' -print -quit)" ] || fail "fixture runs did not log under POPLOG_ROOT=$plr (isolation not exercised)"
# Both mktemp guards: rerun this selftest (plain, no logging, no recursion) with a mktemp that fails on call N.
if [ -z "${MKTEMP_GUARD_INNER:-}" ]; then
  real=$(command -v mktemp); mkdir -p "$t/mkshim"
  # shellcheck disable=SC2016 # shim body expands when the shim runs, not here
  printf '#!/bin/sh\nn=$(( $(cat "$MKSHIM_COUNT") + 1 )); echo "$n" > "$MKSHIM_COUNT"\n[ "$n" -lt "$MKSHIM_FAIL" ] && exec "%s" "$@"\nexit 1\n' "$real" > "$t/mkshim/mktemp"
  chmod +x "$t/mkshim/mktemp"
  for n in 1 2; do
    echo 0 > "$t/mkcount"
    MKTEMP_GUARD_INNER=1 RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1 MKSHIM_COUNT="$t/mkcount" MKSHIM_FAIL=$n PATH="$t/mkshim:$PATH" \
      bash "$here/toolchain_check_selftest.sh" > "$t/mkout" 2>&1 && fail "mktemp failure on call $n not caught"
    ck; grep -q 'FAIL: mktemp$' "$t/mkout" || fail "mktemp failure on call $n: wrong failure: $(cat "$t/mkout")"
    ck; [ "$(cat "$t/mkcount")" = "$n" ] || fail "mktemp failure on call $n: ran on after it ($(cat "$t/mkcount") calls)"
  done
fi
st_done toolchain_check 14   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "toolchain_check selftest: ok"
