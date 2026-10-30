#!/usr/bin/env bash
# usage: tools/fresh_test_run.sh CMD...   (full_gate's build-debug/build-release: CMD is the cargo test run)
#        tools/fresh_test_run.sh --claim  (full_gate's target-dir step: only the target-dir checks below)
# Runs CMD with its output passed through, then fails unless CMD passed, it ran at least one test binary ("Running ... (PATH)" lines), and
# every binary it ran lies inside CARGO_TARGET_DIR (required; must belong to this worktree, see below). Catches a test step that ran a binary from elsewhere (another
# worktree's target). No age test (Spence's ruling, 2026-10-23): the gate may reuse a target dir; per-job subdirs are optional.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap fresh_test_run "$@"   # run log (tools/runlog.sh)
set -u
[ $# -gt 0 ] || { echo "usage: tools/fresh_test_run.sh CMD... | --claim"; exit 2; }
# One target dir per worktree (Paul, 2026-10-23): CARGO_TARGET_DIR must be set; the first worktree to use it writes its path to
# $CARGO_TARGET_DIR/.worktree-owner, and any other worktree on that dir fails (SHARED). Only an empty or missing dir is claimed;
# a non-empty dir with no marker fails (UNMARKED). The worktree is the git top level of the cwd.
[ -n "${CARGO_TARGET_DIR:-}" ] || { echo "fresh_test_run: FAIL CARGO_TARGET_DIR is unset (one target dir per worktree: $HOME/<agent>-work/target-<agent>-<worktree>)"; exit 1; }
wt=$(git rev-parse --show-toplevel 2> /dev/null) && wt=$(cd "$wt" && pwd -P) || { echo "fresh_test_run: FAIL $(pwd) is not inside a git worktree"; exit 1; }
mkdir -p "$CARGO_TARGET_DIR" && td=$(cd "$CARGO_TARGET_DIR" && pwd -P) || { echo "fresh_test_run: FAIL cannot create target dir $CARGO_TARGET_DIR"; exit 1; }
mk=$td/.worktree-owner
if [ -f "$mk" ]; then
  owner=$(cat "$mk")
  [ "$owner" = "$wt" ] || { echo "fresh_test_run: FAIL SHARED target dir $td is owned by worktree $owner, not $wt (one target dir per worktree)"; exit 1; }
elif [ -n "$(ls -A "$td")" ]; then
  # Spence (2026-10-23): a non-empty dir with no marker may hold another worktree's (or a pre-rule) build; never adopt it.
  echo "fresh_test_run: FAIL UNMARKED target dir $td is not empty and has no .worktree-owner marker (use a fresh per-worktree dir)"; exit 1
else
  printf '%s\n' "$wt" > "$mk" || { echo "fresh_test_run: FAIL cannot write $mk"; exit 1; }
fi
[ "$1" = --claim ] && { echo "fresh_test_run: target dir $td belongs to worktree $wt"; exit 0; }
out=$(mktemp "$(zb_scratch)"/zb-fresh.XXXXXX) || exit 2
"$@" 2>&1 | tee "$out"; rc=${PIPESTATUS[0]}
bad=0 n=0
while IFS= read -r p; do
  n=$((n + 1))
  real=$(cd "$(dirname "$p")" 2> /dev/null && pwd -P)/$(basename "$p")
  case "$real" in "$td"/*) ;; *) echo "fresh_test_run: FOREIGN binary $p is not inside $td"; bad=1; continue ;; esac
  [ -f "$real" ] || { echo "fresh_test_run: binary $p missing"; bad=1; continue; }
done < <(sed -n 's/^ *Running .*(\(\/[^)]*\))$/\1/p' "$out")
rm -f "$out"
[ "$n" -gt 0 ] || { echo "fresh_test_run: no test binary ran"; bad=1; }
[ "$rc" = 0 ] || { echo "fresh_test_run: test command failed (exit $rc)"; exit "$rc"; }
[ "$bad" = 0 ] || exit 1
echo "fresh_test_run: $n test binaries, all inside $td"
