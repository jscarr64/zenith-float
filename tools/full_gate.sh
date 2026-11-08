#!/usr/bin/env bash
# Canonical full gate. usage: tools/full_gate.sh TAG
# Every step runs under `timeout` (STEP_TIMEOUT, default 20m; the mpfr step uses MPFR_STEP_TIMEOUT, default 40m); the rc of each step is written to $T/TAG.rc and its
# output to $T/TAG-<step>.log (T, else ZF_GATE_LOG_DIR, else $HOME/rusty-work/logs/zf-gate/<sha>-<time>; never /tmp). REPO selects the checkout (default: this one).
# CARGO_TARGET_DIR is REQUIRED and belongs to one worktree (Paul, 2026-10-23: /tmp/target-<agent>-<worktree>). The first
# step, target-dir, runs the tools/fresh_test_run.sh selftest and --claim: it fails if CARGO_TARGET_DIR is unset or owned
# by another worktree (marker $CARGO_TARGET_DIR/.worktree-owner). Every cargo test step (mpfr-tests-build included) runs
# through tools/fresh_test_run.sh: build with --no-run, then fail on zero test binaries or one outside CARGO_TARGET_DIR,
# before any test runs; only then the tests run.
#
# Clone ratchet "in progress" mode: `rules-scan` scans the added lines since the ratchet started (RATCHET_BASE,
# default c74e834^ = the parent of the first commit with rules_scan and the baseline) and enforces the baseline
# count; `no-clone-progress` prints the remaining no_clone.sh violations (a tracked count, it never fails the gate).
# The zero-clone goal is unchanged: the baseline in tools/clone_count_baseline.txt only goes down.
# `rules-scan` also enforces (via rules_scan.py count checks) the helper-copy baseline (tools/helper_copy_count_baseline.txt),
# the exception allowlist (stale entries fail, tools/exception_count_baseline.txt) and the library-scoped ref ratchets;
# and the neg_via_mul use-count baseline (tools/neg_via_mul_count_baseline.txt, all scanned code incl. tests; it only goes down).
# the clone count includes test code. `no-clone-progress` also prints the helper-copy count. All baselines only go down.
set -u
gate_dir=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
T=${T:-${ZF_GATE_LOG_DIR:-$HOME/rusty-work/logs/zf-gate/$(git -C "$gate_dir" rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)}}; mkdir -p "$T"
cd "${REPO:-$(cd "$(dirname "$0")/.." && pwd)}" || exit 2
tag=${1:?usage: tools/full_gate.sh TAG}; : > "$T/$tag.rc"
# ---- parallel mode (Jeff 2026-11-07 via Paul: up to 4 top-level popslot leaves) ----------------------------------
# ZF_GATE_PARALLEL=N (1..4) makes this script a lightweight COORDINATOR: it holds no slot, runs no step itself and
# launches each leaf of ZF_GATE_LEAVES as its own top-level `poprun --agent A --job JOB-LEAF -- popslot --agent A --
# tools/full_gate.sh TAG-LEAF` with ZF_GATE_ONLY=<the leaf's steps> (run() skips every other step), a private poplog
# (inherited POPLOG_FILE/STATE_DIR/OUTER_* logger state is not passed on, so no leaf shares the coordinator's out/step/
# child files), its own CARGO_TARGET_DIR=$ZF_GATE_TARGET_ROOT/A-JOB/LEAF and T=$T/LEAF. Slot-holder variables are
# never cleared: inside a live slot holder (POPSLOT_HELD / live POPSLOT_HOLDER_PID) it refuses (exit 2); run it outside
# a slot. Unset/empty ZF_GATE_PARALLEL defaults to 4 on Pop (poprun/popslot present) and to the serial gate off Pop; 0 is an explicit serial opt-out; a leaf (ZF_GATE_ONLY set) never takes the default.
# Fails closed: the leaves must cover every `run` step of this file exactly once (target-dir: in every leaf, it claims
# that leaf's target); every leaf's rc file must list exactly its steps, each =0, plus `done`; every leaf's exit code
# is collected with `if wait PID` (a failed leaf never hides behind another). At most N leaves run at once.
# ZF_GATE_LEAF_TIMEOUT (s, default 14400) bounds a leaf; TERM/INT/HUP stop launches and TERM (then after
# ZF_GATE_KILL_GRACE s, KILL) only the coordinator's own leaf process groups, and reap them.
# Target cleanup: the run root A-JOB must not exist (mkdir, never -p; an existing dir is refused and left alone). It is
# removed only when this run created it (owner marker = this pid + start time), every leaf is reaped, no process is
# left in a leaf's process group, no popslot scope named in a leaf log is still active and no process has its cwd or
# an fd under it; otherwise it is kept and reported. ZF_GATE_KEEP_TARGETS=1 always keeps it.
# Leaves share one worktree, so a step that writes a tracked file gets a per-leaf path instead: COMPARE_OUT (prepublish's
# compare-bench timings, normally doc/compare-results.tsv) goes to $T/LEAF. A leaf's prepublish starts from an empty target and
# shares the CPUs, so its timeout is PREPUBLISH_STEP_TIMEOUT (default 40m in leaves, like mpfr; the serial gate keeps STEP_TIMEOUT).
ZF_GATE_LEAVES_DEFAULT='debug:target-dir golden-sums arm-tests-release qemu-time-factor build-debug zero-test build-nostd build-thumbv7em features-random-serde rules-scan no-clone-progress fmt clippy-x86_64 clippy-i686 audit no-intrinsics hex-32bit unsafe-scan
release:target-dir build-release
mpfr:target-dir mpfr-tests-build oracle-fast oracle-deep mpfr
prepublish:target-dir prepublish'
zfp_say() { printf 'zf-parallel: %s\n' "$*"; [[ -n ${zfp_rc_file:-} ]] && printf 'parallel: %s\n' "$*" >> "$zfp_rc_file"; return 0; }
zfp_start() { local s; s=$(cat "/proc/$1/stat" 2>/dev/null) || return 1; s=${s##*) }; set -- $s; echo "${20}"; }
zfp_zombie_or_gone() { local s; s=$(cat "/proc/$1/stat" 2>/dev/null) || return 0; s=${s##*) }; [[ ${s%% *} == Z ]]; }
zfp_group_alive() { pgrep -g "$1" >/dev/null 2>&1; }
zfp_kill_leaf() { # PGID: TERM the leaf's process group, then KILL after the grace period
  local g=$1 i=0 grace=${ZF_GATE_KILL_GRACE:-30}
  kill -TERM -- "-$g" 2>/dev/null; kill -CONT -- "-$g" 2>/dev/null
  while zfp_group_alive "$g" && (( i < grace * 10 )); do sleep 0.1; i=$((i + 1)); done
  if zfp_group_alive "$g"; then zfp_say "leaf group $g still alive after ${grace}s: KILL"; kill -KILL -- "-$g" 2>/dev/null; fi
  return 0
}
zfp_users_of() { # DIR: pids with a cwd or fd under DIR
  local l t
  for l in /proc/[0-9]*/cwd /proc/[0-9]*/fd/*; do
    t=$(readlink "$l" 2>/dev/null) || continue
    [[ $t == "$1" || $t == "$1"/* ]] && { l=${l#/proc/}; echo "${l%%/*}"; }
  done
  return 0
}
zfp_main() {
  local tag=$1 n=$ZF_GATE_PARALLEL bin=${ZF_GATE_POPBIN:-$HOME/pop-logs/bin} self=${ZF_GATE_SELF:-$gate_dir/full_gate.sh}
  local agent=${ZF_GATE_AGENT:-${POPLOG_AGENT:-}} troot=${ZF_GATE_TARGET_ROOT:-/mnt/hub/cargo-targets}
  local leaves=${ZF_GATE_LEAVES:-$ZF_GATE_LEAVES_DEFAULT} limit=${ZF_GATE_LEAF_TIMEOUT:-14400}
  zfp_rc_file="$T/$tag.rc"
  [[ $n =~ ^[1-4]$ ]] || { zfp_say "REFUSING: ZF_GATE_PARALLEL=$n (1..4: the per-agent popslot cap)"; return 2; }
  [[ $limit =~ ^[0-9]+$ ]] || { zfp_say "REFUSING: ZF_GATE_LEAF_TIMEOUT=$limit"; return 2; }
  [[ $agent =~ ^[a-z][a-z0-9_-]*$ ]] || { zfp_say "REFUSING: agent '$agent' (set ZF_GATE_AGENT or POPLOG_AGENT)"; return 2; }
  if [[ -n ${POPSLOT_HELD:-} ]] || { [[ -n ${POPSLOT_HOLDER_PID:-} ]] && kill -0 "$POPSLOT_HOLDER_PID" 2>/dev/null; }; then
    zfp_say "REFUSING: called inside a held popslot (POPSLOT_HELD=${POPSLOT_HELD:-} POPSLOT_HOLDER_PID=${POPSLOT_HOLDER_PID:-}); leaves would reuse that one slot. Run the coordinator outside popslot: poprun --agent $agent --job JOB -- env ZF_GATE_PARALLEL=$n tools/full_gate.sh TAG"
    return 2
  fi
  [[ -x $bin/poprun && -x $bin/popslot ]] || { zfp_say "REFUSING: no poprun/popslot in $bin (parallel mode is Pop-only; set ZF_GATE_PARALLEL=0 for the serial gate)"; return 2; }
  [[ -d $troot ]] || { zfp_say "REFUSING: target root $troot is not a directory"; return 2; }
  local job=${ZF_GATE_JOB:-zfgate-$tag-$(date +%Y%m%d-%H%M%S)-$$}
  [[ $job =~ ^[A-Za-z0-9._-]+$ && $tag =~ ^[A-Za-z0-9._-]+$ ]] || { zfp_say "REFUSING: job '$job' / tag '$tag' must be [A-Za-z0-9._-]"; return 2; }
  # coverage: every `run` step of the gate exactly once (target-dir in every leaf), no unknown or repeated step
  local -a canon=() names=(); local -A seen=() want=() lsteps_of=(); local line lname s bad=0
  mapfile -t canon < <(grep -oE '^run [A-Za-z0-9_-]+' "$self" | cut -d' ' -f2)
  (( ${#canon[@]} > 0 )) || { zfp_say "REFUSING: no run steps found in $self"; return 2; }
  for s in "${canon[@]}"; do want[$s]=1; done
  while IFS= read -r line; do
    [[ -z $line ]] && continue
    lname=${line%%:*}; [[ $line == *:* && $lname =~ ^[a-z0-9-]+$ && -z ${lsteps_of[$lname]:-} ]] || { zfp_say "REFUSING: bad or repeated leaf line '$line'"; return 2; }
    names+=("$lname"); lsteps_of[$lname]=${line#*:}
    [[ " ${line#*:} " == *" target-dir "* ]] || { zfp_say "REFUSING: leaf $lname lacks target-dir (it claims the leaf's target)"; bad=1; }
    for s in ${line#*:}; do
      [[ -n ${want[$s]:-} ]] || { zfp_say "REFUSING: leaf $lname: unknown step $s"; bad=1; continue; }
      [[ $s == target-dir ]] && continue
      [[ -n ${seen[$s]:-} ]] && { zfp_say "REFUSING: step $s in leaves ${seen[$s]} and $lname"; bad=1; }
      seen[$s]=$lname
    done
  done <<< "$leaves"
  for s in "${canon[@]}"; do [[ $s == target-dir || -n ${seen[$s]:-} ]] || { zfp_say "REFUSING: step $s is in no leaf"; bad=1; }; done
  (( bad == 0 && ${#names[@]} > 0 )) || return 2
  local root=$troot/$agent-$job
  mkdir -- "$root" 2>/dev/null || { zfp_say "REFUSING: target root $root exists or cannot be created (never reused, never deleted)"; return 2; }
  local me; me="$$ $(zfp_start $$) $job"; printf '%s\n' "$me" > "$root/.zf-gate-owner"
  zfp_say "coordinator pid $$ job $job leaves ${names[*]} parallel $n root $root"
  local -A pid_of=() t0_of=() rc_of=() why_of=() log_of=(); local -a queue=("${names[@]}") live=()
  zfp_cancel=0
  trap 'zfp_cancel=129' HUP; trap 'zfp_cancel=130' INT; trap 'zfp_cancel=143' TERM
  local p l now
  while (( ${#queue[@]} > 0 || ${#live[@]} > 0 )); do
    if (( zfp_cancel )) && (( ${#queue[@]} > 0 )); then for l in "${queue[@]}"; do rc_of[$l]=$zfp_cancel; why_of[$l]="not started (cancelled)"; done; queue=(); fi
    while (( ! zfp_cancel && ${#queue[@]} > 0 && ${#live[@]} < n )); do
      l=${queue[0]}; queue=("${queue[@]:1}")
      mkdir -p "$T/$l"
      setsid env -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE -u POPLOG_STEP \
        -u POPLOG_OUTER_LOG -u POPLOG_OUTER_PID -u POPLOG_OUTER_DEPTH -u POPLOG_OUTER_START -u POPLOG_OUTER_AGENT \
        -u ZF_GATE_PARALLEL ZF_GATE_ONLY="${lsteps_of[$l]}" CARGO_TARGET_DIR="$root/$l" T="$T/$l" COMPARE_OUT="$T/$l/compare-results.tsv" PREPUBLISH_STEP_TIMEOUT="${PREPUBLISH_STEP_TIMEOUT:-40m}" REPO="$PWD" POPLOG_AGENT="$agent" \
        "$bin/poprun" --agent "$agent" --job "$job-$l" --repo "$PWD" -- "$bin/popslot" --agent "$agent" -- bash "$self" "$tag-$l" \
        > "$T/$l.leaf.out" 2>&1 < /dev/null &
      p=$!; pid_of[$l]=$p; t0_of[$l]=$(date +%s); live+=("$l"); zfp_say "leaf $l started pid/pgid $p steps: ${lsteps_of[$l]}"
    done
    sleep 1
    local -a still=()
    for l in "${live[@]}"; do
      p=${pid_of[$l]}
      if zfp_zombie_or_gone "$p"; then
        if wait "$p"; then rc_of[$l]=0; else rc_of[$l]=$?; fi
        zfp_group_alive "$p" && { zfp_say "leaf $l exited with processes left in group $p"; zfp_kill_leaf "$p"; }
        zfp_say "leaf $l exit ${rc_of[$l]} after $(( $(date +%s) - t0_of[$l] ))s"
        continue
      fi
      now=$(date +%s)
      if (( zfp_cancel )); then why_of[$l]="cancelled (SIG $((zfp_cancel - 128)))"; zfp_kill_leaf "$p"
      elif (( now - t0_of[$l] > limit )); then why_of[$l]="timeout ${limit}s"; zfp_say "leaf $l over ${limit}s: stopping"; zfp_kill_leaf "$p"; fi
      still+=("$l")
    done
    live=("${still[@]}")
  done
  trap - HUP INT TERM
  # verdict: every leaf exit 0 and its rc file lists exactly its steps, each =0, and done
  local fail=0 f got
  for l in "${names[@]}"; do
    f=$T/$l/$tag-$l.rc
    log_of[$l]=$(grep -oE 'log [^ ]+\.log' "$T/$l.leaf.out" 2>/dev/null | tail -1 | cut -d' ' -f2)
    [[ -n ${log_of[$l]} ]] && zfp_say "leaf $l poplog ${log_of[$l]}"
    if [[ ${rc_of[$l]:-x} != 0 ]]; then zfp_say "LEAF FAIL $l exit=${rc_of[$l]:-none} ${why_of[$l]:-}"; fail=1; continue; fi
    got=$(grep -E '^[A-Za-z0-9_-]+=[0-9]+$' "$f" 2>/dev/null | tr '\n' ' ')
    local exp=''; for s in ${lsteps_of[$l]}; do exp+="$s=0 "; done
    # rc-file order is the gate's run order; compare as sets of exact lines
    if [[ "$(tr ' ' '\n' <<<"$got" | sed '/^$/d' | sort)" != "$(tr ' ' '\n' <<<"$exp" | sed '/^$/d' | sort)" ]] || ! grep -qx done "$f" 2>/dev/null; then
      zfp_say "LEAF FAIL $l: rc file $f has [${got}] (want [${exp}] and done)"; fail=1
    else zfp_say "leaf $l PASS: ${got}"; fi
  done
  # cleanup only what this run made, only when nothing can still write into it
  local keep=''
  [[ ${ZF_GATE_KEEP_TARGETS:-0} == 1 ]] && keep='ZF_GATE_KEEP_TARGETS=1'
  [[ "$(cat "$root/.zf-gate-owner" 2>/dev/null)" == "$me" ]] || keep+=' owner marker changed'
  for l in "${names[@]}"; do
    [[ -n ${pid_of[$l]:-} ]] && zfp_group_alive "${pid_of[$l]}" && keep+=" group ${pid_of[$l]} alive"
    for s in $(grep -ohE 'popslot-[0-9-]+\.scope' "$T/$l.leaf.out" ${log_of[$l]:+"${log_of[$l]}"} 2>/dev/null | sort -u); do
      systemctl --user is-active --quiet "$s" 2>/dev/null && keep+=" scope $s active"
    done
  done
  got=$(zfp_users_of "$root" | tr '\n' ' '); [[ -n $got ]] && keep+=" in use by pids $got"
  if [[ -n $keep ]]; then zfp_say "KEEPING $root ($keep)"; else rm -rf -- "$root" && zfp_say "removed own target root $root"; fi
  if (( zfp_cancel )); then zfp_say "RESULT: CANCELLED"; return "$zfp_cancel"; fi
  if (( fail )); then zfp_say "RESULT: FAIL"; return 1; fi
  zfp_say "RESULT: PASS (${#names[@]} leaves)"; return 0
}
# Parallel is the default on Pop (Jeff 2026-11-07): unset/empty ZF_GATE_PARALLEL means 4 when
# poprun/popslot exist, serial otherwise; 0 is an explicit serial opt-out. A leaf (ZF_GATE_ONLY
# set) never takes the default, or it would re-enter the coordinator (POPSLOT_HELD refuses there,
# but the leaf must stay serial on its own terms).
if [[ -z ${ZF_GATE_PARALLEL:-} && -z ${ZF_GATE_ONLY:-} && -x ${ZF_GATE_POPBIN:-$HOME/pop-logs/bin}/poprun && -x ${ZF_GATE_POPBIN:-$HOME/pop-logs/bin}/popslot ]]; then ZF_GATE_PARALLEL=4; fi
if [[ -n ${ZF_GATE_PARALLEL:-} && ${ZF_GATE_PARALLEL:-} != 0 ]]; then zfp_main "$tag"; exit $?; fi
# shellcheck disable=SC1091
if [[ -f "$HOME/pop-logs/bin/poplog.sh" ]]; then source "$HOME/pop-logs/bin/poplog.sh"; else source "$gate_dir/../scripts/poplog.sh"; fi
poplog_init "${POPLOG_AGENT:-rusty}" "full_gate-$tag" "$PWD" "tools/full_gate.sh $tag"
echo "poplog: $POPLOG_FILE" | tee -a "$T/$tag.rc" >&2
# Pop scheduling (Paul, 2026-10-23): heavy gates share 3 counting-lock slots (/tmp/pop-gate.slot1..3, replaces the
# single flock /tmp/pop-gate.lock), 8 jobs each (24 of 32 cores). The slot is held on an fd until this script exits;
# waiting shows as WAITING in ~/pop-logs/RUNNING.tsv. Skipped when a caller already holds a slot (POPSLOT_HELD).
# CARGO_BUILD_JOBS / RUST_TEST_THREADS default to the slot value POPSLOT_JOBS (4) (pop_slot_acquire) unless the caller sets them.
if [[ -z "${POPSLOT_HELD:-}" ]]; then pop_slot_acquire "${POPSLOT_SLOTS:-8}"; fi
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-${POPSLOT_JOBS:-4}} RUST_TEST_THREADS=${RUST_TEST_THREADS:-${POPSLOT_JOBS:-4}}
TT="-- --test-threads=$RUST_TEST_THREADS"
# Preflight: every external tool a step needs must exist; a missing tool fails the gate loudly instead of
# letting a step's check read "command not found" as "no matches".
missing=()
for tool in bash timeout cargo rustup git grep python3 tail; do command -v "$tool" >/dev/null 2>&1 || missing+=("$tool"); done
cargo audit --version >/dev/null 2>&1 || missing+=("cargo-audit")
if [[ "${#missing[@]}" -gt 0 ]]; then
  echo "error: missing tools: ${missing[*]}" | tee -a "$T/$tag.rc" >&2
  echo "preflight=2" >> "$T/$tag.rc"; POPLOG_WARNINGS+=("preflight: missing ${missing[*]}"); poplog_finish 2; exit 2
fi
# pipefail: a failing command inside a step's pipeline fails the step
# Jeff's logging rule (2026-10-23): every step is logged via poplog (~/pop-logs/bin/poplog.sh, else the vendored
# scripts/poplog.sh): header (cmd, repo, SHA, rustc, CARGO_TARGET_DIR), per-step start/end/duration/rc/test counts,
# 5-min heartbeats, RUNNING.tsv row, PASS/FAIL footer. A test step reporting zero tests FAILs. Exit is nonzero on FAIL.
TEST_STEPS=" build-debug build-release features-random-serde oracle-fast oracle-deep mpfr "
run() {
  # parallel leaf (ZF_GATE_ONLY from the coordinator, see zfp_main): only this leaf's steps run, in the usual order
  if [[ -n ${ZF_GATE_ONLY:-} && " $ZF_GATE_ONLY " != *" $1 "* ]]; then return 0; fi
  name=$1; shift; echo "start $name $(date +%T)" >> "$T/$tag.rc"
  rm -f "$(_pl_out)"  # the last step left it a symlink to its log; log_step_start truncates it
  log_step_start "$name"; ln -sf "$T/$tag-$name.log" "$(_pl_out)"  # nested: own out file, never the parent's
  to=${STEP_TIMEOUT:-20m}; [[ $name == mpfr ]] && to=${MPFR_STEP_TIMEOUT:-40m}; [[ $name == prepublish ]] && to=${PREPUBLISH_STEP_TIMEOUT:-$to}  # mpfr: integration_tests alone is ~1090 s
  timeout "$to" bash -o pipefail -c "$*" > "$T/$tag-$name.log" 2>&1 &
  POPLOG_CHILD_PID=$!; echo "$POPLOG_CHILD_PID" > "$POPLOG_STATE_DIR/child"
  wait "$POPLOG_CHILD_PID"; rc=$?; POPLOG_CHILD_PID=""
  echo "$name=$rc" >> "$T/$tag.rc"
  tests=""; [[ "$TEST_STEPS" == *" $name "* ]] && tests=--tests
  log_step_end "$name" "$rc" "$T/$tag-$name.log" "$tests" || true
  echo "  output: $T/$tag-$name.log" >> "$POPLOG_FILE"
  # forward the step's test-binary and summary lines to stdout (full output stays in $T) so a wrapping
  # fresh_test_run.sh (CMD... mode) sees which binaries ran and checks they lie inside CARGO_TARGET_DIR
  grep -E '^ *Running .*\(/[^)]*\)$|^test result: ' "$T/$tag-$name.log" || true
}
ORACLES="--test oracle --test specials_oracle --test methods_oracle --test accuracy_table --test elliptic_phi"
run target-dir "tools/fresh_test_run_selftest.sh && tools/fresh_test_run.sh --claim"
# goldens pinned by SHA256SUMS next to them (fails closed: missing/empty sums, unlisted golden, drift)
run golden-sums "tools/golden_sums_selftest.sh && tools/golden_sums_check.sh zenith-float-num/tests/oracle/goldens borrow_ops.jsonl"
# arm-tests.yml must run the workspace tests with --release, like build-release below (fails closed: missing workflow, step or flag)
run arm-tests-release "tools/test_arm_tests_release.sh"
# the qemu-aarch64 ZF_TIME_FACTOR pin must stay at or above the measured s9_erf ratio + margin (fails closed)
run qemu-time-factor "tools/test_qemu_aarch64_time_factor.sh"
run build-debug "tools/fresh_test_run.sh cargo test --workspace $TT"
# early compile of the mpfr-tests integration tests (6ab75ce prepublish broke here only at the last step)
run mpfr-tests-build "tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests --no-run"
# Zero-test check. Spence review of PR #18: the check was only in
# scripts/ci.sh, which set -e stops at the ci_no_name failure before reaching
# it, and full_gate never called it at all. Placed here, not between
# build-debug and mpfr-tests-build: tools/test_full_gate_poplog.sh asserts
# mpfr-tests-build follows build-debug directly. The debug --list still
# reuses the build-debug binaries (mpfr-tests-build is a release --no-run),
# so this does not force an extra debug rebuild; do not move it after
# build-release for the same reason.
run zero-test "bash scripts/check_zero_test_targets.sh"
run build-release "tools/fresh_test_run.sh cargo test --workspace --release $TT"
run build-nostd "cargo build --no-default-features && cargo build -p zenith-float-num --no-default-features"
run build-thumbv7em "cargo build --no-default-features --target thumbv7em-none-eabihf && cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf"
run features-random-serde "tools/fresh_test_run.sh cargo test --workspace --features random,serde $TT"
run oracle-fast "tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests $ORACLES -- --test-threads=1"
run oracle-deep "ORACLE_TIER=deep tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests $ORACLES -- --test-threads=1"
# MPFR_TIER selects the MPFR tier (tests/mpfr/common.rs): the gate runs `full` (every run count, the
# 100032-bit tests) unless the caller sets MPFR_TIER=fast.
run mpfr "MPFR_TIER=${MPFR_TIER:-full} tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests -- --test-threads=1"
# ratchet mode: the clone ratchet started at c74e834 (first commit with rules_scan + baseline); older added lines are tracked by no-clone-progress, not flagged here
run rules-scan "tools/rules_scan_selftest.sh && tools/rules_scan.sh ${RATCHET_BASE:-c74e834^}..HEAD"
# informational: remaining violations are a tracked count (must only go down; hard fail is for the final step)
run no-clone-progress "tools/no_clone_selftest.sh && { tools/no_clone.sh | tail -1; tools/no_clone.sh --helpers | tail -1; true; }"
run fmt "cargo fmt --all -- --check"
run clippy-x86_64 "cargo clippy --workspace --all-targets --all-features -- -D warnings"
run clippy-i686 "cargo clippy -p zenith-float-num --all-targets --target i686-unknown-linux-musl -- -D warnings"
run audit "cargo audit --deny warnings"
run no-intrinsics "bash scripts/ci_no_intrinsics.sh"
run hex-32bit "bash scripts/ci_hex_32bit.sh"
run unsafe-scan "bash tools/gate_unsafe_scan.sh"
run prepublish "bash scripts/zenith_prepublish.sh"
echo "done" >> "$T/$tag.rc"
poplog_finish
exit "$POPLOG_OVERALL"
