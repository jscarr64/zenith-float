#!/usr/bin/env bash
# Red test for Paul's rule (2026-10-23, one target dir per worktree) and tools/fresh_test_run.sh:
# an unset CARGO_TARGET_DIR, a target dir owned by another worktree, a foreign (outside the target dir) test
# binary, a missing binary and zero binaries each FAIL; a real cargo fixture builds, passes the binary check
# before any test runs, then runs its test; every cargo test step in tools/full_gate.sh and every cargo test run in
# scripts/zenith_prepublish.sh goes through the script.
# usage: tools/fresh_test_run_selftest.sh   (one tiny throwaway crate under mktemp; never the repo's target dir)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
ftr="$root/tools/fresh_test_run.sh"
t=$(mktemp -d); trap 'rm -rf "$t"' EXIT
git init -q "$t/wa"; git init -q "$t/wb"
wa=$(cd "$t/wa" && pwd -P)
bad=0
fail() { echo "FAIL: $*" >&2; bad=1; }
# want pass|fail NAME DIR CMD...: run CMD in DIR (a worktree) and check its exit status
want() {
  local w=$1 name=$2 dir=$3 rc=0; shift 3
  (cd "$dir" && "$@") > "$t/out" 2>&1 || rc=$?
  if { [[ $w == pass ]] && [[ $rc -ne 0 ]]; } || { [[ $w == fail ]] && [[ $rc -eq 0 ]]; }; then
    fail "$name wanted $w, rc=$rc"; sed 's/^/  | /' "$t/out" >&2
  else echo "fresh_test_run_selftest: $name $w ok"; fi
}
has() { grep -c -- "$1" "$t/out" > /dev/null || fail "$2: output lacks '$1'"; }
td="$t/target"
want fail unset-dir "$t/wa" env -u CARGO_TARGET_DIR "$ftr" --claim
has 'CARGO_TARGET_DIR is unset' unset-dir
want fail not-a-worktree "$t" env CARGO_TARGET_DIR="$td" "$ftr" --claim
want pass claim "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" --claim
[[ "$(cat "$td/.worktree-owner" 2>/dev/null)" == "$wa" ]] || fail "no owner marker naming $wa"
want pass claim-again "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" --claim
want fail shared-dir "$t/wb" env CARGO_TARGET_DIR="$td" "$ftr" --claim
has 'SHARED' shared-dir; has "$wa" shared-dir
# Spence's hole: a non-empty dir without the marker (stale binaries, or a deleted marker) must not be claimed
mkdir -p "$t/stale/debug/deps"; printf '#!/bin/sh\n' > "$t/stale/debug/deps/harness-deadbeef"
want fail unmarked-stale "$t/wa" env CARGO_TARGET_DIR="$t/stale" "$ftr" --claim
has 'UNMARKED' unmarked-stale
[[ ! -e "$t/stale/.worktree-owner" ]] || fail "unmarked-stale: marker written"
want pass claim-empty "$t/wa" env CARGO_TARGET_DIR="$t/empty-target" "$ftr" --claim
(cd "$t/wa" && CARGO_TARGET_DIR="$t/tdm" "$ftr" --claim > /dev/null); mkdir -p "$t/tdm/debug"
rm "$t/tdm/.worktree-owner"
want fail marker-deleted "$t/wb" env CARGO_TARGET_DIR="$t/tdm" "$ftr" --claim
has 'UNMARKED' marker-deleted
# a caller-set override target dir takes the same check (fails before any build)
want fail override-flag-shared "$t/wb" env CARGO_TARGET_DIR="$t/wb-own" "$ftr" cargo test --target-dir "$td"
has 'SHARED' override-flag-shared
want fail override-eq-unmarked "$t/wb" env CARGO_TARGET_DIR="$t/wb-own" "$ftr" cargo test --target-dir="$t/stale"
has 'UNMARKED' override-eq-unmarked
want fail override-env-shared "$t/wb" env CARGO_TARGET_DIR="$t/wb-own" CARGO_BUILD_TARGET_DIR="$td" "$ftr" --claim
has 'SHARED' override-env-shared
mkdir -p "$td/debug/deps" "$t/elsewhere"
for b in "$td/debug/deps/t-1" "$t/elsewhere/t-2"; do printf '#!/bin/sh\n' > "$b"; chmod +x "$b"; done
printf '%s\n' "$td/debug/deps/t-1" > "$t/good.txt"
printf '%s\n%s\n' "$td/debug/deps/t-1" "$t/elsewhere/t-2" > "$t/foreign.txt"
printf '%s\n' "$td/debug/deps/gone" > "$t/missing.txt"
: > "$t/empty.txt"
want pass bins-inside "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" --bins "$t/good.txt"
want fail foreign-binary "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" --bins "$t/foreign.txt"
has 'FOREIGN' foreign-binary
want fail missing-binary "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" --bins "$t/missing.txt"
want fail zero-binaries "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" --bins "$t/empty.txt"
want fail shared-bins "$t/wb" env CARGO_TARGET_DIR="$td" "$ftr" --bins "$t/good.txt"
want fail no-args "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr"
# real cargo: build, check the binaries, then run; --no-run stops after the check; a failing test fails
mkdir -p "$t/wa/src"
printf '[package]\nname = "ftr_fixture"\nversion = "0.0.0"\nedition = "2021"\n[workspace]\n' > "$t/wa/Cargo.toml"
printf '#[test]\nfn passes() {}\n#[test]\nfn fails_when_asked() { assert!(std::env::var_os("FTR_FAIL").is_none()); }\n' > "$t/wa/src/lib.rs"
want pass cargo-run "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" cargo test --lib -- --test-threads=1
has 'test binaries, all inside' cargo-run; has 'test result: ok. 2 passed' cargo-run
want pass cargo-no-run "$t/wa" env CARGO_TARGET_DIR="$td" "$ftr" cargo test --lib --no-run
grep -c 'test result:' "$t/out" > /dev/null && fail "cargo-no-run ran tests"
want fail cargo-red-test "$t/wa" env FTR_FAIL=1 CARGO_TARGET_DIR="$td" "$ftr" cargo test --lib -- --test-threads=1
want fail cargo-shared "$t/wb" env CARGO_TARGET_DIR="$td" "$ftr" cargo test --manifest-path "$t/wa/Cargo.toml" --lib
if grep -c 'test binaries' "$t/out" > /dev/null; then fail "cargo-shared built before the claim failed"; fi
# wiring: target-dir runs --claim first; every cargo test step goes through tools/fresh_test_run.sh
g="$root/tools/full_gate.sh"
first=$(grep -m1 -E '^run [a-z0-9-]+ "' "$g")
[[ $first == 'run target-dir '*'tools/fresh_test_run.sh --claim'* ]] || fail "first gate step is not target-dir --claim: $first"
while IFS= read -r line; do
  [[ $line == *'tools/fresh_test_run.sh cargo test'* ]] || fail "cargo test step not through fresh_test_run: $line"
done < <(grep -E '^run [a-z0-9-]+ ".*cargo test' "$g")
# the prepublish step (scripts/zenith_prepublish.sh) runs cargo test itself: each one goes through the script too
pp="$root/scripts/zenith_prepublish.sh"; n=0
while IFS= read -r line; do
  n=$((n + 1))
  [[ $line == *'tools/fresh_test_run.sh" cargo test'* ]] || fail "prepublish cargo test not through fresh_test_run: $line"
done < <(grep -E '^[^#]*cargo test( |$)' "$pp" | grep -v '^ *step ')
[[ $n -ge 3 ]] || fail "prepublish: expected >= 3 cargo test runs, found $n"
[[ $bad == 0 ]] && echo "fresh_test_run_selftest: ok"
exit "$bad"
