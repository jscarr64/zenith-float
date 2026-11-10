#!/usr/bin/env bash
. "${BASH_SOURCE[0]%/*}/scratch.sh"   # zb_scratch: temp files and dirs under $HOME, never the system temp dir
# zb's thin layer over the shared Pop logger ~/pop-logs/bin/poplog.sh (Rusty's; see ~/pop-logs/README.md; Jeff's rules of 9:42 and 9:59 AM CT
# 2026-10-23, and 10:33 AM CT: rely on poplog, no workarounds). poplog writes the log (header, steps, counts, heartbeat, footer), owns the
# canonical RUNNING.tsv row and its lock, tees step output to stdout, joins nested runs to the parent's log (SUB-RUN sections) and fails the
# run (footer FAIL, non-zero exit) when any step or sub-run failed. This file only adds what zb needs on top:
#   - a wrapped script registers its own exit cleanup with runlog_on_exit CMD, never a raw `trap ... EXIT`, and never execs a command
#   - runlog_wrap opens the run (or the sub-run, when started inside poprun or another wrapped script) and, at exit, finishes it and exits
#     non-zero if any step failed, even when the script itself ends with exit 0;
#   - runlog_step always passes --tests when asked, so a test step that ran zero tests is FAIL;
#   - runlog_step --tests also enforces zb's own 0-tests=FAIL rule (sum of "test result:" counts), with or without poplog;
#   - a missing or unloadable poplog.sh is a HARD failure (exit 2): no run without a run log. RUNLOG_PLAIN=1 (self-test fixtures only) runs
#     the steps plain on purpose, and only with the selftest-only marker RUNLOG_SELFTEST=1 too (exported only by self-test fixtures);
#     RUNLOG_PLAIN=1 alone is refused (exit 2), so a real gate cannot be run unlogged by a stray RUNLOG_PLAIN.
# Use, first thing in a script:   . "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap JOB "$@"
#   runlog_step [--tests] NAME -- CMD...     runlog_warn MESSAGE
# One-shot (verify and test runs): tools/runlog.sh JOB [--tests] -- CMD...
# Env: RUNLOG_AGENT (else poprun's POPLOG_AGENT, else popslot's POPSLOT_AGENT; none: refuse), RUNLOG_POPLOG (default ~/pop-logs/bin/poplog.sh); poplog's POPLOG_ROOT, POPLOG_HEARTBEAT, POPLOG_ECHO apply.
RUNLOG_AGENT=${RUNLOG_AGENT:-${POPLOG_AGENT:-${POPSLOT_AGENT:-}}}   # from the caller, never a hardcoded agent (Spence 2026-10-24)
RUNLOG_POPLOG=${RUNLOG_POPLOG:-$HOME/pop-logs/bin/poplog.sh}
RUNLOG_ON=""
RUNLOG_EXIT_CMDS=()

_rl_exit() { # EXIT trap: finish the (sub-)run; a failed step makes the script exit non-zero
  local rc=$? c
  for c in "${RUNLOG_EXIT_CMDS[@]}"; do eval "$c"; done
  [ -n "$RUNLOG_ON" ] || exit "$rc"
  poplog_finish "$rc"
  if [ "$rc" = 0 ] && [ "${POPLOG_OVERALL:-0}" != 0 ]; then exit 1; fi
  exit "$rc"
}

runlog_wrap() { # runlog_wrap JOB "$@"
  local job=$1; shift
  trap _rl_exit EXIT
  set -o pipefail   # Spence review 5: a failing stage before a pipe must not be masked by the last stage
  if [ "${RUNLOG_PLAIN:-0}" = 1 ]; then   # fixture-only plain mode: refused unless the selftest-only marker is also set
    if [ "${RUNLOG_SELFTEST:-0}" != 1 ]; then
      echo "runlog: FAIL RUNLOG_PLAIN=1 without RUNLOG_SELFTEST=1; plain mode (no run log) is for self-test fixtures only, refusing to run $job" >&2; exit 2
    fi
    return 0
  fi
  [ -n "$RUNLOG_AGENT" ] || { echo "runlog: FAIL no agent for $job: set RUNLOG_AGENT or run under poprun/popslot --agent" >&2; exit 2; }
  if [ ! -f "$RUNLOG_POPLOG" ]; then echo "runlog: FAIL $RUNLOG_POPLOG missing; refusing to run $job without a run log" >&2; exit 2; fi
  # shellcheck source=/dev/null
  . "$RUNLOG_POPLOG" || { echo "runlog: FAIL cannot source $RUNLOG_POPLOG; refusing to run $job without a run log" >&2; exit 2; }
  poplog_init "$RUNLOG_AGENT" "$job" "$(git rev-parse --show-toplevel 2> /dev/null || pwd)" "$0 $*"
  trap _rl_exit EXIT
  RUNLOG_ON=1
  [ "${POPLOG_NESTED:-0}" = 1 ] || echo "runlog: $POPLOG_FILE" >&2
}

runlog_on_exit() { # runlog_on_exit CMD: run CMD at exit, before the run log is finished (a raw `trap ... EXIT` would replace the finish)
  RUNLOG_EXIT_CMDS+=("$1")
}

runlog_detach_children() { # for self-tests: fixture runs (some fail on purpose) open their own logs instead of joining this run as sub-runs
  export -n POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE 2> /dev/null || true
}

runlog_warn() {
  if [ -n "$RUNLOG_ON" ]; then POPLOG_WARNINGS+=("$*"); _pl_log "WARNING $*"; fi
  echo "warning: $*" >&2
}

_rl_tee() { local f=$1; shift; "$@" 2>&1 | tee "$f"; return "${PIPESTATUS[0]}"; }
_rl_ran() { awk '/test result: /{for(i=1;i<NF;i++) if($(i+1) ~ /^(passed|failed);?$/) s+=$i} END{print s+0}' "$1"; }

runlog_step() { # runlog_step [--tests] NAME -- CMD...
  local tflag=() name f rc n
  [ "${1:-}" = --tests ] && { tflag=(--tests); shift; }
  name=$1; shift; [ "${1:-}" = -- ] && shift
  if [ "${#tflag[@]}" = 0 ]; then
    if [ -z "$RUNLOG_ON" ]; then "$@"; return $?; fi
    poplog_run_step "$name" -- "$@"; return $?
  fi
  # zb layer: a test step must run at least one test, independently of poplog's own check
  f=$(mktemp "$(zb_scratch)"/zb-runlog-tests.XXXXXX) || return 2
  if [ -z "$RUNLOG_ON" ]; then _rl_tee "$f" "$@"; else poplog_run_step "$name" --tests -- _rl_tee "$f" "$@"; fi; rc=$?
  n=$(_rl_ran "$f"); rm -f "$f"
  if [ "$rc" = 0 ] && [ "$n" = 0 ]; then echo "runlog: step $name ran 0 tests: FAIL" >&2; return 1; fi
  return "$rc"
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  job=${1:?usage: tools/runlog.sh JOB [--tests] -- CMD...}; shift
  tf=(); [ "${1:-}" = --tests ] && { tf=(--tests); shift; }; [ "${1:-}" = -- ] && shift
  runlog_wrap "$job" "$@"
  runlog_step "${tf[@]}" "$job" -- "$@"
  exit $?
fi
