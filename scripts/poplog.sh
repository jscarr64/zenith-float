#!/usr/bin/env bash
# poplog.sh - shared Pop!_OS run-logging library (Jeff's rule, 2026-10-23).
# Source it, then:
#   poplog_init AGENT JOB [REPO] [CMD-DESCRIPTION]   header + RUNNING.tsv row + heartbeat + traps
#   log_step_start NAME                              step start line
#   poplog_run_step NAME [--tests] -- cmd...         run cmd as a step (live output, counts, pause detection)
#   log_counts FILE                                  parse cargo test / nextest summaries -> COUNT_* vars
#   log_step_end NAME RC [FILE] [--tests]            end line, duration, rc, counts; 0 tests on a test step = FAIL
#   poplog_progress "record 12/500"                  current record/progress shown by the heartbeat
#   poplog_finish [RC]                               footer PASS/FAIL, removes RUNNING row
# Env: POPLOG_ROOT (~/pop-logs), POPLOG_HEARTBEAT (300 s), POPLOG_PROGRESS_FILE (a step may write it directly).
# This file must not use `set -e` effects on callers; it works under `set -euo pipefail`.

POPLOG_ROOT=${POPLOG_ROOT:-$HOME/pop-logs}
POPLOG_HEARTBEAT=${POPLOG_HEARTBEAT:-300}
POPLOG_RUNNING="$POPLOG_ROOT/RUNNING.tsv"
POPLOG_OVERALL=0
POPLOG_WARNINGS=()
POPLOG_DONE=0
POPLOG_HB_PID=""
POPLOG_CHILD_PID=""

POPLOG_RUNNING_HEADER=$'#agent\tjob\tpid\tlog\tstart\tlast_heartbeat\tstate\tslot'
_pl_iso() { date '+%Y-%m-%dT%H:%M:%S%:z'; }
_pl_now() { date '+%Y-%m-%d %H:%M:%S %Z'; }
_pl_log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >> "$POPLOG_FILE"; }
_pl_fmt() { local s=$1; printf '%dh%02dm%02ds' $((s / 3600)) $((s % 3600 / 60)) $((s % 60)); }

_pl_running_add() {
  mkdir -p "$POPLOG_ROOT"
  (
    flock -x 9
    _pl_running_migrate
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$POPLOG_AGENT" "$POPLOG_JOB" "$POPLOG_PID" "$POPLOG_FILE" \
      "$(_pl_iso)" "$(_pl_iso)" RUNNING - >> "$POPLOG_RUNNING"
  ) 9> "$POPLOG_RUNNING.lock"
}
# Canonical RUNNING.tsv (see README): header line, 8 tab-separated columns, ISO-8601 local times with offset.
# _pl_running_migrate (caller holds the lock): ensure the canonical header; pad old 7-column rows with slot "-".
_pl_running_migrate() {
  if [[ ! -s "$POPLOG_RUNNING" ]]; then printf '%s\n' "$POPLOG_RUNNING_HEADER" > "$POPLOG_RUNNING"; return 0; fi
  [[ "$(head -1 "$POPLOG_RUNNING")" == "$POPLOG_RUNNING_HEADER" ]] && return 0
  awk -F'\t' -v OFS='\t' -v h="$POPLOG_RUNNING_HEADER" 'NR==1 {print h; if ($1 !~ /^#?agent$/) {if (NF<8) $8="-"; print}; next} {if (NF<8) $8="-"; print}' \
    "$POPLOG_RUNNING" > "$POPLOG_RUNNING.tmp.$$" && mv "$POPLOG_RUNNING.tmp.$$" "$POPLOG_RUNNING"
}
# _pl_running_update STATE [SLOT]  (refresh last_heartbeat + state (+slot) for our pid)
_pl_running_update() {
  (
    flock -x 9
    [[ -f "$POPLOG_RUNNING" ]] || exit 0
    awk -F'\t' -v OFS='\t' -v p="$POPLOG_PID" -v t="$(_pl_iso)" -v s="$1" -v sl="${2:-}" \
      'NR>1 && $3==p {$6=t; $7=s; if (sl != "") $8=sl} {print}' "$POPLOG_RUNNING" > "$POPLOG_RUNNING.tmp.$$" \
      && mv "$POPLOG_RUNNING.tmp.$$" "$POPLOG_RUNNING"
  ) 9> "$POPLOG_RUNNING.lock"
}
_pl_running_remove() {
  (
    flock -x 9
    [[ -f "$POPLOG_RUNNING" ]] || exit 0
    awk -F'\t' -v p="$POPLOG_PID" 'NR==1 || $3!=p' "$POPLOG_RUNNING" > "$POPLOG_RUNNING.tmp.$$" \
      && mv "$POPLOG_RUNNING.tmp.$$" "$POPLOG_RUNNING"
  ) 9> "$POPLOG_RUNNING.lock"
  rm -f "$POPLOG_ROOT/STALE/$POPLOG_PID.flag"
}

# poplog_progress TEXT - record current progress for the heartbeat.
poplog_progress() { printf '%s\n' "$*" > "$POPLOG_PROGRESS_FILE"; }

# Descendant pids of $1 (inclusive).
_pl_tree() { local p=$1 c; echo "$p"; for c in $(pgrep -P "$p" 2>/dev/null); do _pl_tree "$c"; done; }
_pl_paused() {
  local p st
  [[ -n "${1:-}" ]] || return 1
  for p in $(_pl_tree "$1"); do
    st=$(awk '{print $3}' "/proc/$p/stat" 2>/dev/null || true)
    [[ "$st" == T || "$st" == t ]] && return 0
  done
  return 1
}

_pl_heartbeat_loop() {
  local step cur prog child state
  while sleep "$POPLOG_HEARTBEAT"; do
    kill -0 "$POPLOG_PID" 2>/dev/null || exit 0
    step=$(cat "$POPLOG_STATE_DIR/step" 2>/dev/null || echo '-')
    child=$(cat "$POPLOG_STATE_DIR/child" 2>/dev/null || true)
    cur=$(grep -a -E '^test .* \.\.\. ' "$POPLOG_STATE_DIR/out" 2>/dev/null | tail -1 || true)
    [[ -n "$cur" ]] || cur=$(grep -a -E '^ *(PASS|FAIL|START) \[' "$POPLOG_STATE_DIR/out" 2>/dev/null | tail -1 || true)
    prog=$(tail -1 "$POPLOG_PROGRESS_FILE" 2>/dev/null || true)
    state=RUNNING
    [[ -f "$POPLOG_STATE_DIR/waiting" ]] && state=WAITING
    if _pl_paused "$child" || _pl_paused "$POPLOG_PID"; then state=PAUSED; fi
    _pl_log "HEARTBEAT state=$state step=$step elapsed=$(_pl_fmt $(( $(date +%s) - POPLOG_T0 ))) current_test=${cur:--} progress=${prog:--}"
    [[ "$state" == PAUSED ]] && _pl_log "PAUSED (SIGSTOP detected; process state T)"
    _pl_running_update "$state"
  done
}

_pl_on_signal() {
  local sig=$1 code=$2
  [[ "$POPLOG_DONE" == 1 ]] && exit "$code"
  if [[ -n "$POPLOG_CHILD_PID" ]]; then kill -TERM "$POPLOG_CHILD_PID" 2>/dev/null || true; kill -CONT "$POPLOG_CHILD_PID" 2>/dev/null || true; fi
  _pl_log "KILLED by SIG$sig (exit $code)"
  poplog_finish "$code"
  exit "$code"
}
_pl_on_exit() {
  local rc=$?
  [[ "$POPLOG_DONE" == 1 ]] && return
  _pl_log "EXIT without poplog_finish (rc=$rc)"
  [[ $rc -eq 0 ]] && rc=1
  poplog_finish "$rc"
}

poplog_init() {
  POPLOG_AGENT=$1; POPLOG_JOB=$2; POPLOG_REPO=${3:-$PWD}; local desc=${4:-${0##*/}}
  POPLOG_PID=$$
  POPLOG_T0=$(date +%s)
  local day; day=$(date +%F)
  mkdir -p "$POPLOG_ROOT/$POPLOG_AGENT/$day"
  # nested runs (poprun -> full_gate.sh) get their own log; never reuse an inherited POPLOG_FILE
  POPLOG_FILE=$POPLOG_ROOT/$POPLOG_AGENT/$day/$(date +%H%M)-$POPLOG_JOB.log
  [[ -e "$POPLOG_FILE" ]] && POPLOG_FILE="${POPLOG_FILE%.log}-$$.log"
  POPLOG_STATE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/poplog.XXXXXX")
  POPLOG_PROGRESS_FILE=$POPLOG_STATE_DIR/progress
  export POPLOG_FILE POPLOG_PROGRESS_FILE POPLOG_PID POPLOG_STATE_DIR POPLOG_ROOT POPLOG_RUNNING
  : > "$POPLOG_STATE_DIR/out"; echo '-' > "$POPLOG_STATE_DIR/step"
  local sha='(not a git repo)'
  sha=$(git -C "$POPLOG_REPO" rev-parse HEAD 2>/dev/null || echo '(not a git repo)')
  local dirty=''; [[ -n "$(git -C "$POPLOG_REPO" status --porcelain 2>/dev/null | head -1)" ]] && dirty=' (dirty)'
  {
    echo "===== POPLOG HEADER ====="
    echo "agent:            $POPLOG_AGENT"
    echo "job:              $POPLOG_JOB"
    echo "command:          $desc"
    echo "repo:             $POPLOG_REPO"
    echo "commit:           $sha$dirty"
    echo "rustc:            $(rustc --version 2>/dev/null || echo 'not found')"
    echo "CARGO_TARGET_DIR: ${CARGO_TARGET_DIR:-(unset)}"
    echo "host/pid:         $(hostname)/$POPLOG_PID"
    echo "start:            $(_pl_now)"
    echo "heartbeat:        every ${POPLOG_HEARTBEAT}s; progress file $POPLOG_PROGRESS_FILE"
    echo "========================="
  } > "$POPLOG_FILE"
  _pl_running_add
  _pl_heartbeat_loop & POPLOG_HB_PID=$!
  trap '_pl_on_signal INT 130' INT
  trap '_pl_on_signal TERM 143' TERM
  trap '_pl_on_signal HUP 129' HUP
  trap '_pl_on_exit' EXIT
}

log_step_start() {
  export POPLOG_STEP=$1; POPLOG_STEP_T0=$(date +%s)
  echo "$1" > "$POPLOG_STATE_DIR/step"; : > "$POPLOG_STATE_DIR/out"; : > "$POPLOG_PROGRESS_FILE"
  _pl_log "STEP START $1"
}

# log_counts FILE - sums every cargo-test and nextest summary in FILE.
log_counts() {
  local f=$1 line
  COUNT_PASSED=0; COUNT_FAILED=0; COUNT_IGNORED=0; COUNT_SUMMARIES=0
  while IFS= read -r line; do
    if [[ "$line" =~ test\ result:\ [A-Za-z]+\.\ ([0-9]+)\ passed\;\ ([0-9]+)\ failed\;\ ([0-9]+)\ ignored ]]; then
      COUNT_PASSED=$((COUNT_PASSED + BASH_REMATCH[1])); COUNT_FAILED=$((COUNT_FAILED + BASH_REMATCH[2]))
      COUNT_IGNORED=$((COUNT_IGNORED + BASH_REMATCH[3])); COUNT_SUMMARIES=$((COUNT_SUMMARIES + 1))
    elif [[ "$line" =~ Summary\ \[.*\]\ +([0-9]+)\ tests?\ run: ]]; then
      COUNT_SUMMARIES=$((COUNT_SUMMARIES + 1))
      [[ "$line" =~ ([0-9]+)\ passed ]] && COUNT_PASSED=$((COUNT_PASSED + BASH_REMATCH[1]))
      [[ "$line" =~ ([0-9]+)\ failed ]] && COUNT_FAILED=$((COUNT_FAILED + BASH_REMATCH[1]))
      [[ "$line" =~ ([0-9]+)\ skipped ]] && COUNT_IGNORED=$((COUNT_IGNORED + BASH_REMATCH[1]))
    fi
  done < <(sed 's/\x1b\[[0-9;]*m//g' "$f")
  COUNT_RUN=$((COUNT_PASSED + COUNT_FAILED))
  return 0
}

# log_step_end NAME RC [OUTFILE] [--tests]
log_step_end() {
  local name=$1 rc=$2 f=${3:-$POPLOG_STATE_DIR/out} tests=${4:-} verdict=PASS dur
  dur=$(( $(date +%s) - ${POPLOG_STEP_T0:-$POPLOG_T0} ))
  log_counts "$f"
  [[ $rc -ne 0 ]] && verdict=FAIL
  if [[ "$tests" == --tests || $COUNT_SUMMARIES -gt 0 ]] && [[ $((COUNT_RUN + COUNT_IGNORED)) -eq 0 ]]; then
    verdict=FAIL; POPLOG_WARNINGS+=("step $name: ZERO tests reported"); _pl_log "WARNING step $name reported zero tests -> FAIL"
  fi
  [[ $COUNT_FAILED -gt 0 ]] && verdict=FAIL
  _pl_log "STEP END $name verdict=$verdict exit=$rc duration=$(_pl_fmt "$dur") tests: run=$COUNT_RUN passed=$COUNT_PASSED failed=$COUNT_FAILED ignored=$COUNT_IGNORED (summaries=$COUNT_SUMMARIES)"
  if [[ $verdict == FAIL ]]; then POPLOG_OVERALL=1; POPLOG_WARNINGS+=("step $name FAIL (exit $rc)"); return 1; fi
  return 0
}

# poplog_run_step NAME [--tests] -- cmd...   (tests auto-detected for cargo test / nextest)
poplog_run_step() {
  local name=$1 tests=''; shift
  [[ "${1:-}" == --tests ]] && { tests=--tests; shift; }
  [[ "${1:-}" == -- ]] && shift
  [[ "$*" =~ (cargo[[:space:]]+(test|nextest)|nextest[[:space:]]+run) ]] && tests=--tests
  log_step_start "$name"
  _pl_log "CMD $*"
  local rc=0
  "$@" > >(tee -a "$POPLOG_STATE_DIR/out" >> "$POPLOG_FILE") 2>&1 &
  POPLOG_CHILD_PID=$!; echo "$POPLOG_CHILD_PID" > "$POPLOG_STATE_DIR/child"
  wait "$POPLOG_CHILD_PID" || rc=$?
  POPLOG_CHILD_PID=""; : > "$POPLOG_STATE_DIR/child"
  sleep 0.2  # let tee drain
  log_step_end "$name" "$rc" "$POPLOG_STATE_DIR/out" "$tests" || return 1
}

# pop_slot_acquire [N] - take the first free of N (default 3) counting-lock slots $POPSLOT_PREFIX1..N
# (default /tmp/pop-gate.slot). Holds it on fd $POPSLOT_FD until this process exits. If none is free the run is
# WAITING (RUNNING.tsv + heartbeat); waiters queue on $POPSLOT_PREFIX.queue (one poller at a time, FIFO-ish) and
# poll every $POPSLOT_POLL s (default 3). Sets POPSLOT_SLOT and POPSLOT_WAITED. Works in or out of a poplog run.
_ps_note() { if [[ -n "${POPLOG_FILE:-}" && -f "${POPLOG_FILE:-}" ]]; then _pl_log "$*"; fi; echo "popslot: $*" >&2; }
_ps_try() {
  local i fd
  for ((i = 1; i <= POPSLOT_N; i++)); do
    exec {fd}>> "${POPSLOT_PREFIX}$i"
    if flock -n "$fd"; then POPSLOT_FD=$fd; POPSLOT_SLOT=$i; return 0; fi
    exec {fd}>&-
  done
  return 1
}
pop_slot_acquire() {
  POPSLOT_N=${1:-${POPSLOT_SLOTS:-3}}; POPSLOT_PREFIX=${POPSLOT_PREFIX:-/tmp/pop-gate.slot}
  local t0 qfd ctx=''; t0=$(date +%s)
  [[ -n "${POPLOG_PID:-}" && -n "${POPLOG_STATE_DIR:-}" && -d "${POPLOG_STATE_DIR:-/nonexistent}" ]] && ctx=1
  if ! _ps_try; then
    _ps_note "SLOT WAITING: all $POPSLOT_N slots (${POPSLOT_PREFIX}1..$POPSLOT_N) busy"
    if [[ -n "$ctx" ]]; then : > "$POPLOG_STATE_DIR/waiting"; _pl_running_update WAITING "wait/$POPSLOT_N"; fi
    exec {qfd}>> "${POPSLOT_PREFIX}.queue"
    flock "$qfd"
    until _ps_try; do sleep "${POPSLOT_POLL:-3}"; done
    exec {qfd}>&-
  fi
  POPSLOT_WAITED=$(( $(date +%s) - t0 ))
  export POPSLOT_SLOT POPSLOT_FD POPSLOT_HELD=$POPSLOT_SLOT
  # Per-gate budget (Paul/Jeff 2026-10-23): 3 slots x 8 jobs = 24 of 32 cores, unless the caller overrides.
  export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-8} RUST_TEST_THREADS=${RUST_TEST_THREADS:-8}
  if [[ -n "$ctx" ]]; then rm -f "$POPLOG_STATE_DIR/waiting"; _pl_running_update RUNNING "$POPSLOT_SLOT/$POPSLOT_N"; fi
  _ps_note "SLOT TAKEN slot=$POPSLOT_SLOT/$POPSLOT_N (${POPSLOT_PREFIX}$POPSLOT_SLOT) waited=$(_pl_fmt "$POPSLOT_WAITED") CARGO_BUILD_JOBS=$CARGO_BUILD_JOBS RUST_TEST_THREADS=$RUST_TEST_THREADS"
}

poplog_finish() {
  local rc=${1:-$POPLOG_OVERALL} w
  [[ "$POPLOG_DONE" == 1 ]] && return 0
  POPLOG_DONE=1
  if [[ -n "$POPLOG_HB_PID" ]]; then kill "$POPLOG_HB_PID" 2>/dev/null || true; fi
  [[ $rc -ne 0 ]] && POPLOG_OVERALL=1
  {
    echo "===== POPLOG FOOTER ====="
    if [[ $POPLOG_OVERALL -eq 0 ]]; then echo "RESULT: PASS"; else echo "RESULT: FAIL (exit $rc)"; fi
    echo "end:      $(_pl_now)"
    echo "duration: $(_pl_fmt $(( $(date +%s) - POPLOG_T0 )))"
    echo "warnings: ${#POPLOG_WARNINGS[@]}"
    for w in "${POPLOG_WARNINGS[@]}"; do echo "  - $w"; done
    echo "========================="
  } >> "$POPLOG_FILE"
  _pl_running_remove
  rm -rf "$POPLOG_STATE_DIR"
  return 0
}
