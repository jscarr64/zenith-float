#!/usr/bin/env bash
# One target dir per worktree (Paul, 2026-10-23; ported from zenith-bound 9791ccc). A gate run uses its own
# CARGO_TARGET_DIR (/tmp/target-<agent>-<worktree>) and runs only test binaries built inside it.
# usage: tools/fresh_test_run.sh --claim
#          CARGO_TARGET_DIR must be set and the cwd inside a git worktree; the first worktree to use the dir is
#          written to $CARGO_TARGET_DIR/.worktree-owner and any other worktree on that dir FAILs (SHARED);
#          a non-empty dir without that marker FAILs (UNMARKED). Override dirs (CARGO_BUILD_TARGET_DIR,
#          --target-dir) take the same check.
#        tools/fresh_test_run.sh cargo test ARGS... [-- TEST-ARGS...]
#          --claim, then build with `cargo test ARGS --no-run --message-format=json`, then the binary check:
#          at least one test binary, every one exists and lies inside CARGO_TARGET_DIR (else FOREIGN), all
#          BEFORE any test runs; then the original command runs the tests (skipped when ARGS has --no-run).
#        tools/fresh_test_run.sh --bins FILE
#          --claim, then the binary check on FILE (one test binary path per line); used by the selftest.
# Self-test: tools/fresh_test_run_selftest.sh. No new tools: bash, git, cargo, python3 (all in full_gate's preflight).
set -euo pipefail
die() { echo "fresh_test_run: FAIL $*" >&2; exit 1; }
usage="usage: tools/fresh_test_run.sh --claim | --bins FILE | cargo test ARGS... [-- TEST-ARGS...]"
[[ $# -gt 0 ]] || die "$usage"
case $1 in
  --claim) [[ $# -eq 1 ]] || die "$usage" ;;
  --bins) [[ $# -eq 2 ]] || die "$usage" ;;
  cargo) [[ ${2:-} == test ]] || die "$usage" ;;
  *) die "$usage" ;;
esac
[[ -n ${CARGO_TARGET_DIR:-} ]] || die "CARGO_TARGET_DIR is unset (one target dir per worktree: /tmp/target-<agent>-<worktree>)"
wt=$(git rev-parse --show-toplevel 2> /dev/null) || die "$(pwd) is not inside a git worktree"
wt=$(cd "$wt" && pwd -P)
# claim_dir DIR: print DIR's real path if this worktree owns it. A dir with no marker is claimed only when
# absent or empty; a non-empty unmarked dir (stale binaries, or a deleted marker) FAILs (UNMARKED).
claim_dir() {
  local d=$1 r m owner
  if [[ -d $d && ! -f $d/.worktree-owner && -n $(ls -A "$d") ]]; then
    die "UNMARKED target dir $d is not empty and has no .worktree-owner (stale or re-claimed; use a fresh dir)"
  fi
  mkdir -p "$d" || die "cannot create target dir $d"
  r=$(cd "$d" && pwd -P)
  m=$r/.worktree-owner
  if [[ -f $m ]]; then
    owner=$(cat "$m")
    [[ $owner == "$wt" ]] || die "SHARED target dir $r is owned by worktree $owner, not $wt (one target dir per worktree)"
  else
    printf '%s\n' "$wt" > "$m" || die "cannot write $m"
  fi
  printf '%s\n' "$r"
}
td=$(claim_dir "$CARGO_TARGET_DIR") || exit 1
# a caller-set override target dir (CARGO_BUILD_TARGET_DIR, --target-dir) takes the same ownership check;
# --target-dir wins over the env in cargo, so the binary check then uses it
if [[ -n ${CARGO_BUILD_TARGET_DIR:-} ]]; then claim_dir "$CARGO_BUILD_TARGET_DIR" > /dev/null || exit 1; fi
prev=
for a in "${@:2}"; do
  [[ $a == -- ]] && break
  o=
  if [[ $prev == --target-dir ]]; then o=$a; elif [[ $a == --target-dir=* ]]; then o=${a#--target-dir=}; fi
  if [[ -n $o ]]; then td=$(claim_dir "$o") || exit 1; fi
  prev=$a
done
if [[ $1 == --claim ]]; then echo "fresh_test_run: target dir $td belongs to worktree $wt"; exit 0; fi
# check_bins FILE: >= 1 binary, each present and inside $td
check_bins() {
  local p real n=0
  while IFS= read -r p; do
    [[ -n $p ]] || continue
    n=$((n + 1))
    [[ -f $p ]] || die "test binary $p missing"
    real=$(cd "$(dirname "$p")" && pwd -P)/$(basename "$p")
    case $real in "$td"/*) ;; *) die "FOREIGN test binary $real is not inside $td" ;; esac
  done < "$1"
  [[ $n -gt 0 ]] || die "no test binaries (zero binaries is a failure)"
  echo "fresh_test_run: $n test binaries, all inside $td"
}
if [[ $1 == --bins ]]; then [[ -r $2 ]] || die "cannot read $2"; check_bins "$2"; exit 0; fi
# cargo test: split at the first `--` so --no-run --message-format=json go to cargo, not the test harness
pre=() no_run=0 seen=0
for a in "${@:3}"; do
  if [[ $seen == 0 && $a == -- ]]; then seen=1; fi
  if [[ $seen == 0 ]]; then
    if [[ $a == --no-run ]]; then no_run=1; else pre+=("$a"); fi
  fi
done
w=$(mktemp -d); trap 'rm -rf "$w"' EXIT
echo "fresh_test_run: build: cargo test ${pre[*]} --no-run --message-format=json"
rc=0; cargo test "${pre[@]}" --no-run --message-format=json-render-diagnostics > "$w/build.json" || rc=$?
[[ $rc == 0 ]] || die "build failed (exit $rc)"
python3 -c '
import json, sys
for line in open(sys.argv[1]):
    try:
        m = json.loads(line)
    except ValueError:
        continue
    if m.get("reason") == "compiler-artifact" and m.get("executable") and m.get("profile", {}).get("test"):
        print(m["executable"])
' "$w/build.json" | sort -u > "$w/bins.txt"
check_bins "$w/bins.txt"
[[ $no_run == 1 ]] && exit 0
exec "${@:1}"
