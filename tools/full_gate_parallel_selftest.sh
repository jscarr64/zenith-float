#!/usr/bin/env bash
# Selftest for tools/full_gate.sh parallel mode and the parallel default on Pop.
# Harper's coordinator harness (2026-11-07: 32 checks, 5 mutants killed) adapted into the
# repo by Vera (2026-11-08) as a self-contained file, plus the default-dispatch checks
# (unset on Pop, ZF_GATE_PARALLEL=0, unset off Pop, bad value, leaf guard).
# Default (Jeff 2026-11-07): unset/empty ZF_GATE_PARALLEL means 4 on Pop (poprun/popslot
# present) and the serial gate off Pop; 0 opts out to serial; a leaf (ZF_GATE_ONLY set)
# never takes the parallel default.
# usage: tools/full_gate_parallel_selftest.sh [GATE] [WORKDIR]
#   GATE defaults to the full_gate.sh beside this script; WORKDIR defaults to a mktemp dir
#   (removed on exit). Fake poprun/popslot and a fake leaf gate are generated under WORKDIR;
#   no real gate steps run in the coordinator checks; the serial checks and the leaf-order
#   check run the real gate's preflight and selected steps only, in a fixture, with contained
#   poplog state and no slot. On Pop, a fixture-repo poprun exercises the clean-tree check, and
#   a brief real popslot is taken only when this process does not already hold a slot.
#   Run on Pop (the serial checks need the gate's toolchain). Fast (minutes).
set -u
gate=${1:-$(cd "$(dirname "$0")" && pwd)/full_gate.sh}
if [[ -n ${2:-} ]]; then W=$2; mkdir "$W" || { echo "workdir exists"; exit 2; }
else W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; fi
mkdir -p "$W/bin" "$W/tools" "$W/troot"
cat > "$W/bin/poprun" <<'FIXTURE_POPRUN_EOF'
#!/usr/bin/env bash
job=?; while [[ $# -gt 0 && $1 != -- ]]; do [[ $1 == --job ]] && job=$2; shift; done; shift
echo "poprun: log $FAKE_STATE/fake-$job.log" >&2
echo "OOM guard: MemoryMax=6G MemorySwapMax=0 (own systemd --user scope popslot-1-2.scope)" >&2
: > "$FAKE_STATE/fake-$job.log"
exec "$@"
FIXTURE_POPRUN_EOF
cat > "$W/bin/popslot" <<'FIXTURE_POPSLOT_EOF'
#!/usr/bin/env bash
while [[ $# -gt 0 && $1 != -- ]]; do shift; done; shift
[[ -n ${POPSLOT_HOLDER_PID:-}${POPSLOT_SLOT:-}${POPSLOT_FD:-}${POPSLOT_HELD:-} ]] && echo "inherited holder=${POPSLOT_HOLDER_PID-} slot=${POPSLOT_SLOT-} fd=${POPSLOT_FD-} held=${POPSLOT_HELD-}" >> "$FAKE_STATE/holder-leak"
exec "$@"
FIXTURE_POPSLOT_EOF
chmod +x "$W/bin/"*
cat > "$W/fakegate-head" <<'FIXTURE_FAKEGATE_EOF'
#!/usr/bin/env bash
# fake leaf gate: behaviour per leaf from FAKE_<leaf>
leaf=${CARGO_TARGET_DIR##*/}; tag=$1; mkdir -p "$T" "$CARGO_TARGET_DIR"
echo "$CARGO_TARGET_DIR" >> "$FAKE_STATE/targets"
env | grep -E '^POPLOG_(FILE|STATE_DIR|OUTER_LOG)=' >> "$FAKE_STATE/logger-leak"
mode_var=FAKE_${leaf//-/_}; mode=${!mode_var:-pass}
( flock 9; n=$(ls "$FAKE_STATE/active" | wc -l); n=$((n+1)); : > "$FAKE_STATE/active/$$"; m=$(cat "$FAKE_STATE/max" 2>/dev/null || echo 0); (( n > m )) && echo "$n" > "$FAKE_STATE/max" ) 9> "$FAKE_STATE/lock"
sleep "${FAKE_SLEEP:-2}"
case $mode in
  sleep) sleep 300 ;;
  escape) ( cd "$CARGO_TARGET_DIR" && setsid sleep 20 < /dev/null > /dev/null 2>&1 & ) ;;
esac
rm -f "$FAKE_STATE/active/$$"
skip=''; [[ $mode == missing ]] && skip=${ZF_GATE_ONLY##* }
for s in $ZF_GATE_ONLY; do [[ $s == "$skip" ]] && continue; echo "$s=0" >> "$T/$tag.rc"; done
[[ $mode == badrc ]] && echo "${ZF_GATE_ONLY##* }=1" >> "$T/$tag.rc"
echo done >> "$T/$tag.rc"
[[ $mode == fail ]] && exit 3
exit 0
FIXTURE_FAKEGATE_EOF
cp "$gate" "$W/tools/full_gate.sh"
{ cat "$W/fakegate-head"; grep -E '^run ' "$gate"; } > "$W/fakegate"
chmod +x "$W/fakegate"
# Fallback poplog for off-Pop runs: the gate sources $HOME/pop-logs/bin/poplog.sh first,
# else $gate_dir/../scripts/poplog.sh; the fixture provides the latter.
if [[ -f "$(cd "$(dirname "$gate")" && pwd)/../scripts/poplog.sh" ]]; then
  mkdir -p "$W/scripts"; cp "$(cd "$(dirname "$gate")" && pwd)/../scripts/poplog.sh" "$W/scripts/poplog.sh"
fi
fails=0; n=0
ok() { n=$((n+1)); echo "ok   $n $*"; }
bad() { n=$((n+1)); fails=$((fails+1)); echo "FAIL $n $*"; }
grep -q zfp_main "$gate" || { echo "FAIL 0 gate has no parallel mode (ZF_GATE_PARALLEL coordinator)"; exit 1; }
SIX='a:target-dir golden-sums arm-tests-release qemu-time-factor parallel-selftest build-debug zero-test
b:target-dir build-nostd build-thumbv7em features-random-serde rules-scan no-clone-progress
c:target-dir fmt clippy-x86_64 clippy-i686 audit no-intrinsics hex-32bit unsafe-scan
release:target-dir build-release
mpfr:target-dir mpfr-tests-build oracle-fast oracle-deep mpfr
prepublish:target-dir prepublish'
# case: fresh state dir, run coordinator; extra env as args.
# Fixture poplog (A1): the coordinator now refuses without a live top-level poplog. Slot-holder
# vars from an outer poprun/popslot are unset so a gate-step run of this selftest does not
# trip the live-holder refusal or leak a slot into the fakes. Tests that need a holder set it
# again in "$@".
coord() { # NAME env...
  local c=$1; shift
  S=$W/state-$c; mkdir -p "$S/active" "$S/plogstate"; : > "$S/targets"; : > "$S/logger-leak"; : > "$S/plog.log"
  ( env -u POPSLOT_HELD -u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE \
      -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_NESTED -u POPLOG_OUTER_LOG -u POPLOG_OUTER_PID \
      FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" ZF_GATE_KILL_GRACE=2 FAKE_SLEEP=2 \
      POPLOG_FILE="$S/plog.log" POPLOG_STATE_DIR="$S/plogstate" POPLOG_NESTED=0 \
      "$@" bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
}
# leftover(): a fake leaf is still running. The `sleep 300` match is scoped to processes
# whose cwd is under this run's workdir: poplog's heartbeat also runs `sleep 300`
# (POPLOG_HEARTBEAT, default 300s) as a child of the poprun shell, with the invocation
# dir as its cwd, and an unscoped match false-positives on it under poprun. A fake leaf's
# `sleep 300` runs with cwd in the workdir (the leaf runs there), so it is still caught.
leftover() {
  pgrep -f "$W/fakegate" >/dev/null && return 0
  local p
  for p in $(pgrep -f 'sleep 300' -u "$(id -u)" 2>/dev/null); do
    [[ $(readlink "/proc/$p/cwd" 2>/dev/null) == "$W"* ]] && return 0
  done
  return 1
}
# 1 all pass, 4 default leaves
coord pass ZF_GATE_PARALLEL=4; rc=$?
[[ $rc == 0 ]] && ok "all leaves pass -> exit 0" || bad "pass exit $rc"
[[ ! -e $W/troot/harper-job-pass ]] && ok "own target root removed after success" || bad "own root not removed"
[[ $(sort "$S/targets" | wc -l) == 4 && $(sort -u "$S/targets" | wc -l) == 4 ]] && ok "4 leaves, 4 distinct targets" || bad "targets: $(tr '\n' ' ' < "$S/targets")"
grep -vq "^$W/troot/harper-job-pass/" "$S/targets" && bad "target outside run root" || ok "every target under the run root"
[[ ! -s $S/logger-leak ]] && ok "no inherited poplog state in leaves" || bad "logger leak: $(cat "$S/logger-leak")"
[[ ! -e $S/holder-leak ]] && ok "no slot holder inherited" || bad "holder leaked into leaf"
[[ $(cat "$S/max") -le 4 ]] && ok "max concurrent $(cat "$S/max") <= 4" || bad "max $(cat "$S/max")"
grep -q 'ROLE: parallel-coordinator' "$S/plog.log" && ok "coordinator wrote ROLE into its poplog" || bad "no coordinator ROLE in poplog"
grep -q 'RESULT: PASS' "$S/plog.log" && ok "coordinator poplog has RESULT: PASS" || bad "no coordinator RESULT in poplog"
# 2 cap 2 with six leaves
coord cap ZF_GATE_PARALLEL=2 "ZF_GATE_LEAVES=$SIX"; rc=$?
[[ $rc == 0 && $(cat "$S/max") == 2 ]] && ok "6 leaves, cap 2: max concurrent exactly 2" || bad "cap: rc $rc max $(cat "$S/max" 2>/dev/null)"
coord cap1 ZF_GATE_PARALLEL=1 "ZF_GATE_LEAVES=$SIX"; rc=$?
[[ $rc == 0 && $(cat "$S/max") == 1 ]] && ok "cap 1 is serial" || bad "cap1: rc $rc max $(cat "$S/max" 2>/dev/null)"
# 3 one leaf fails, others still complete and are reaped
coord fail ZF_GATE_PARALLEL=4 FAKE_release=fail; rc=$?
[[ $rc == 1 ]] && ok "failing leaf -> exit 1" || bad "fail exit $rc"
grep -q 'LEAF FAIL release exit=3' "$W/out-fail" && ok "failed leaf named with its exit 3" || bad "no LEAF FAIL release exit=3"
[[ $(grep -c 'PASS:' "$W/out-fail") == 3 ]] && ok "other 3 leaves still collected" || bad "others not collected"
coord faillast ZF_GATE_PARALLEL=4 FAKE_prepublish=fail FAKE_SLEEP=1; rc=$?
[[ $rc == 1 ]] && ok "failure in last leaf not masked" || bad "faillast exit $rc"
# 4 rc-file contract
coord missing ZF_GATE_PARALLEL=4 FAKE_mpfr=missing; rc=$?
[[ $rc == 1 ]] && grep -q 'LEAF FAIL mpfr: rc file' "$W/out-missing" && ok "leaf exit 0 but missing step -> FAIL" || bad "missing step rc $rc"
coord badrc ZF_GATE_PARALLEL=4 FAKE_debug=badrc; rc=$?
[[ $rc == 1 ]] && ok "step rc 1 inside an exit-0 leaf -> FAIL" || bad "badrc rc $rc"
# 5 timeout of a stuck leaf
coord timeout ZF_GATE_PARALLEL=4 FAKE_prepublish=sleep ZF_GATE_LEAF_TIMEOUT=5; rc=$?
[[ $rc == 1 ]] && grep -q 'timeout 5s' "$W/out-timeout" && ok "stuck leaf timed out -> FAIL" || bad "timeout rc $rc"
leftover && bad "processes left after timeout" || ok "no process left after timeout"
[[ ! -e $W/troot/harper-job-timeout ]] && ok "root removed after timeout" || bad "root kept after timeout"
# 6 cancellation
S=$W/state-cancel; mkdir -p "$S/active" "$S/plogstate"; : > "$S/targets"; : > "$S/logger-leak"; : > "$S/plog.log"
env -u POPSLOT_HELD -u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE \
  -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_NESTED \
  FAKE_STATE="$S" T="$W/T-cancel" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
  ZF_GATE_AGENT=harper ZF_GATE_JOB=job-cancel ZF_GATE_KILL_GRACE=2 ZF_GATE_PARALLEL=2 "ZF_GATE_LEAVES=$SIX" FAKE_a=sleep FAKE_b=sleep \
  POPLOG_FILE="$S/plog.log" POPLOG_STATE_DIR="$S/plogstate" POPLOG_NESTED=0 \
  bash "$W/tools/full_gate.sh" tcancel > "$W/out-cancel" 2>&1 &
cp=$!; sleep 4; kill -TERM "$cp"; wait "$cp"; rc=$?
[[ $rc == 143 ]] && ok "TERM -> exit 143" || bad "cancel exit $rc"
leftover && bad "processes left after cancel" || ok "no leaf process left after cancel"
grep -q 'not started (cancelled)' "$W/out-cancel" && ok "queued leaves not launched after cancel" || bad "queued leaves launched"
[[ $(wc -l < "$S/targets") == 2 ]] && ok "only the 2 running leaves had started" || bad "started $(wc -l < "$S/targets")"
[[ ! -e $W/troot/harper-job-cancel ]] && ok "own root removed after cancel (all reaped)" || bad "root kept after cancel"
# 7 existing root refused and untouched
mkdir "$W/troot/harper-job-exist"; echo keep > "$W/troot/harper-job-exist/sentinel"
coord exist ZF_GATE_PARALLEL=4; rc=$?
[[ $rc == 2 && $(cat "$W/troot/harper-job-exist/sentinel") == keep && ! -s $S/targets ]] && ok "existing root refused, untouched, nothing launched" || bad "exist rc $rc"
# 8 a process still using the root keeps it
coord escape ZF_GATE_PARALLEL=4 FAKE_debug=escape; rc=$?
[[ -d $W/troot/harper-job-escape ]] && grep -q KEEPING "$W/out-escape" && ok "root in use by a foreign process is kept" || bad "in-use root removed (rc $rc)"
for p in $(pgrep -f 'sleep 20' -u "$(id -u)"); do [[ $(readlink "/proc/$p/cwd") == "$W/troot/harper-job-escape"* ]] && kill "$p"; done
# 9 inside a held slot
coord holder ZF_GATE_PARALLEL=4 POPSLOT_HELD=3; rc=$?
[[ $rc == 2 && ! -e $W/troot/harper-job-holder && ! -s $S/targets ]] && grep -q 'outside popslot' "$W/out-holder" && ok "inside a held slot -> refuse, nothing launched" || bad "holder rc $rc"
coord holder2 ZF_GATE_PARALLEL=4 POPSLOT_HOLDER_PID=$$; rc=$?
[[ $rc == 2 ]] && ok "live POPSLOT_HOLDER_PID -> refuse" || bad "holder2 rc $rc"
# 10 coverage
coord nocover ZF_GATE_PARALLEL=4 "ZF_GATE_LEAVES=${SIX/ audit/}"; rc=$?
[[ $rc == 2 && ! -s $S/targets ]] && grep -q 'step audit is in no leaf' "$W/out-nocover" && ok "uncovered step -> refuse" || bad "nocover rc $rc"
coord dup ZF_GATE_PARALLEL=4 "ZF_GATE_LEAVES=$SIX
d:target-dir fmt"; rc=$?
[[ $rc == 2 ]] && grep -q 'step fmt in leaves' "$W/out-dup" && ok "duplicated step -> refuse" || bad "dup rc $rc"
coord unk ZF_GATE_PARALLEL=4 "ZF_GATE_LEAVES=$SIX
e:target-dir bogus"; rc=$?
[[ $rc == 2 ]] && ok "unknown step -> refuse" || bad "unk rc $rc"
coord badn ZF_GATE_PARALLEL=5; rc=$?
[[ $rc == 2 ]] && ok "ZF_GATE_PARALLEL=5 refused (cap 4)" || bad "badn rc $rc"
coord nobin ZF_GATE_PARALLEL=4 ZF_GATE_POPBIN=/nonexistent; rc=$?
[[ $rc == 2 ]] && ok "no poprun/popslot -> refuse" || bad "nobin rc $rc"
# 12 default dispatch (Jeff 2026-11-07: parallel is the default on Pop)
# run the gate with ZF_GATE_PARALLEL unset, like a bare `tools/full_gate.sh TAG` on Pop
coord_unset() { # NAME env...
  local c=$1; shift
  S=$W/state-$c; mkdir -p "$S/active" "$S/plogstate"; : > "$S/targets"; : > "$S/logger-leak"; : > "$S/plog.log"
  ( env -u ZF_GATE_PARALLEL -u POPSLOT_HELD -u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE \
      -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_NESTED \
      FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" ZF_GATE_KILL_GRACE=2 FAKE_SLEEP=2 \
      POPLOG_FILE="$S/plog.log" POPLOG_STATE_DIR="$S/plogstate" POPLOG_NESTED=0 \
      "$@" bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
}
# run the gate directly (no coordinator) with contained poplog state, for the dispatch
# checks that must NOT enter the coordinator. POPSLOT_HELD=1 skips slot acquisition;
# preset POPLOG_FILE/STATE_DIR keep poplog_init inside the fixture (no real log).
serial() { # NAME env...
  local c=$1; shift
  S=$W/state-$c; mkdir -p "$S/active" "$S/popstate"; : > "$S/targets"; : > "$S/logger-leak"; : > "$S/poplog.log"
  ( env -u ZF_GATE_PARALLEL FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" ZF_GATE_KILL_GRACE=2 FAKE_SLEEP=2 POPSLOT_HELD=1 POPLOG_FILE="$S/poplog.log" POPLOG_STATE_DIR="$S/popstate" POPLOG_OUTER_LOG= "$@" bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
}
# 33 unset ZF_GATE_PARALLEL on Pop -> coordinator, 4 leaves
coord_unset unset; rc=$?
[[ $rc == 0 ]] && ok "unset ZF_GATE_PARALLEL on Pop -> coordinator, exit 0" || bad "unset exit $rc"
[[ $(sort "$S/targets" | wc -l) == 4 && $(sort -u "$S/targets" | wc -l) == 4 ]] && ok "unset -> 4 leaves, 4 distinct targets" || bad "unset targets: $(tr '\n' ' ' < "$S/targets")"
[[ $(ls "$S"/fake-*.log 2>/dev/null | wc -l) == 4 ]] && ok "unset -> poprun invoked for 4 leaves" || bad "unset poprun markers: $(ls "$S"/fake-*.log 2>/dev/null | wc -l)"
# 34 ZF_GATE_PARALLEL=0 -> serial, never the coordinator
serial zero ZF_GATE_PARALLEL=0; rc=$?
[[ $rc != 0 ]] && ok "=0 -> serial, fails fast in the fixture (rc $rc)" || bad "=0 exit $rc"
[[ $(ls "$S"/fake-*.log 2>/dev/null | wc -l) == 0 ]] && ok "=0 -> poprun never invoked" || bad "=0 poprun invoked"
grep -q 'REFUSING' "$W/out-zero" && bad "=0 -> coordinator refused instead of serial" || ok "=0 -> no coordinator refusal"
grep -q '^start target-dir' "$W/T-zero/tzero.rc" && ok "=0 -> serial path ran target-dir" || bad "=0 -> no serial target-dir step"
# 35 unset off Pop (no poprun/popslot) -> serial, never a refusal
serial offpop ZF_GATE_POPBIN=/nonexistent; rc=$?
[[ $rc != 0 ]] && ok "unset off Pop -> serial, fails fast in the fixture (rc $rc)" || bad "offpop exit $rc"
[[ $(ls "$S"/fake-*.log 2>/dev/null | wc -l) == 0 ]] && ok "unset off Pop -> poprun never invoked" || bad "offpop poprun invoked"
grep -q 'REFUSING' "$W/out-offpop" && bad "unset off Pop -> refused instead of serial" || ok "unset off Pop -> no refusal"
grep -q '^start target-dir' "$W/T-offpop/toffpop.rc" && ok "unset off Pop -> serial path ran target-dir" || bad "offpop -> no serial target-dir step"
# 36 bad value refused
coord badval ZF_GATE_PARALLEL=abc; rc=$?
[[ $rc == 2 ]] && grep -q 'REFUSING: ZF_GATE_PARALLEL=abc' "$W/out-badval" && ok "ZF_GATE_PARALLEL=abc refused" || bad "badval rc $rc"
# 37 a leaf (ZF_GATE_ONLY set) never takes the parallel default: without this guard a
# leaf would re-enter the coordinator and the gate would break (POPSLOT_HELD refusal).
serial leafguard ZF_GATE_ONLY="target-dir"; rc=$?
[[ $rc != 0 ]] && ok "leaf (ZF_GATE_ONLY) -> serial, fails fast in the fixture (rc $rc)" || bad "leafguard exit $rc"
[[ $(ls "$S"/fake-*.log 2>/dev/null | wc -l) == 0 ]] && ok "leaf -> poprun never invoked" || bad "leafguard poprun invoked"
grep -q 'REFUSING' "$W/out-leafguard" && bad "leaf -> coordinator refused instead of serial" || ok "leaf -> no coordinator refusal"
grep -q '^start target-dir' "$W/T-leafguard/tleafguard.rc" && ok "leaf -> serial path ran target-dir" || bad "leafguard -> no serial target-dir step"
# 11 leaf filter keeps canonical order and only its steps (run() with ZF_GATE_ONLY)
# Use the real copied gate (not the fakegate) on the serial path with a fixture poplog.
order() {
  local c=order
  S=$W/state-$c; mkdir -p "$S/active" "$S/popstate"; : > "$S/targets"; : > "$S/logger-leak"; : > "$S/poplog.log"
  ( env -u ZF_GATE_PARALLEL FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" POPSLOT_HELD=1 POPLOG_FILE="$S/poplog.log" POPLOG_STATE_DIR="$S/popstate" POPLOG_OUTER_LOG= \
      ZF_GATE_ONLY="fmt audit target-dir" bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
}
order; rc=$?
rcf=$W/T-order/torder.rc
got=$(grep -E '^(start |[A-Za-z0-9_-]+=|done$)' "$rcf" 2>/dev/null | tr '\n' '|')
[[ $got == start\ target-dir\ *\|target-dir=*\|start\ fmt\ *\|fmt=*\|start\ audit\ *\|audit=*\|done\| ]] \
  && ok "leaf filter keeps canonical order (target-dir, fmt, audit) and only those steps" \
  || bad "leaf order rc=$rc got [$got]"
grep -q 'start build-debug' "$rcf" 2>/dev/null && bad "leaf order ran a step not in ZF_GATE_ONLY" || ok "leaf filter skipped steps not in ZF_GATE_ONLY"
grep -q 'ROLE: parallel-leaf' "$S/poplog.log" && grep -q 'NOT A GATE: this is a parallel leaf log' "$S/poplog.log" \
  && ok "leaf poplog is marked NOT A GATE / ROLE: parallel-leaf" \
  || bad "leaf poplog missing evidence marker"

# A1: a bare coordinator (no live poplog) is refused and launches nothing
coord_nolog() {
  local c=nolog
  S=$W/state-$c; mkdir -p "$S/active"; : > "$S/targets"; : > "$S/logger-leak"
  ( env -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_NESTED \
      -u POPSLOT_HELD -u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE \
      FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" ZF_GATE_PARALLEL=4 FAKE_SLEEP=1 \
      bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
}
coord_nolog; rc=$?
[[ $rc == 2 && ! -s $S/targets ]] && grep -q 'no live top-level poplog' "$W/out-nolog" \
  && ok "bare coordinator (no poplog) refused, nothing launched" || bad "nolog rc $rc targets=$(tr '\n' ' ' < "$S/targets")"

# A1: a nested sub-run is refused
coord nested POPLOG_NESTED=1; rc=$?
[[ $rc == 2 && ! -s $S/targets ]] && grep -q 'not a nested sub-run' "$W/out-nested" \
  && ok "coordinator inside a nested sub-run refused" || bad "nested rc $rc"

# A2: a stale dead POPSLOT_HOLDER_PID/SLOT/FD (HELD unset, so the coordinator does not refuse)
# must not be inherited by leaves
coord staleholder ZF_GATE_PARALLEL=4 POPSLOT_HOLDER_PID=999999 POPSLOT_SLOT=3 POPSLOT_FD=99; rc=$?
[[ $rc == 0 && ! -e $S/holder-leak ]] && ok "stale dead POPSLOT_HOLDER_PID/SLOT/FD not inherited by leaves" \
  || bad "staleholder rc $rc leak=$(tr '\n' ' ' < "$S/holder-leak" 2>/dev/null)"
grep -q 'REFUSING: called inside a held popslot' "$W/out-staleholder" && bad "stale dead holder refused the coordinator" || ok "stale dead holder did not refuse the coordinator"

# A3: systemctl error is not treated as 'scope inactive'
printf '#!/bin/sh\nexit 1\n' > "$W/bin/systemctl"; chmod +x "$W/bin/systemctl"
coord syserr ZF_GATE_PARALLEL=4 PATH="$W/bin:$PATH"; rc=$?
[[ $rc == 0 ]] && grep -q 'status-unknown (systemctl exit 1)' "$W/out-syserr" && [[ -d $W/troot/harper-job-syserr ]] \
  && ok "systemctl error keeps the target root (fail closed)" \
  || bad "syserr rc $rc out=$(grep -E 'KEEPING|removed' "$W/out-syserr")"
rm -f "$W/bin/systemctl"
# the kept root must not block later jobs; leave it (unique job name)

# A4: traps are still installed when cleanup runs (trap - comes after KEEPING/removed)
# and zfp_kill_leaf requires a recorded start time before signaling
awk_ord=$(awk '
  /zfp_kill_leaf\(\)/ {k=1}
  k && /zfp_start/ {ks=1}
  k && /kill -TERM/ {if(!ks) {print "kill-before-start-check"; exit}}
  /^  trap - HUP INT TERM$/ {if (c==0) early=1; t=NR}
  /KEEPING \$root/ {c=NR}
  END {if(t==0||c==0) print "missing"; else if(early) print "trap-cleared-before-cleanup"; else print "ok"}
' "$gate")
[[ $awk_ord == ok ]] && ok "traps stay until after cleanup; kill guards start time" || bad "script order: $awk_ord"

# A5: the selftest is a required run step of the gate, in the debug leaf, after qemu-time-factor
grep -q '^run parallel-selftest ' "$gate" && ok "full_gate.sh has run parallel-selftest" || bad "parallel-selftest is not a run step"
grep -A1 '^run qemu-time-factor ' "$gate" | grep -q '^run parallel-selftest ' \
  || grep -A2 '^run qemu-time-factor ' "$gate" | grep -q '^run parallel-selftest ' \
  && ok "parallel-selftest is after qemu-time-factor (not between build-debug and mpfr-tests-build)" \
  || bad "parallel-selftest placement"
grep -qE "^ZF_GATE_LEAVES_DEFAULT=.*parallel-selftest" "$gate" && ok "debug leaf lists parallel-selftest" || bad "leaves default missing parallel-selftest"

# Real poprun clean-tree check and real popslot, only when the real binaries exist.
# Clean-tree uses a fixture git repo (never the zenith-float worktree). poprun is invoked
# with POPLOG_* unset so it is a top-level log, and without popslot (no extra slot).
real_bin=${ZF_GATE_REAL_POPBIN:-$HOME/pop-logs/bin}
if [[ -x $real_bin/poprun && -x $real_bin/popslot ]]; then
  fix=$W/fixrepo
  mkdir -p "$fix"
  git -C "$fix" init -q
  git -C "$fix" config user.email harper@local
  git -C "$fix" config user.name Harper
  printf 'fixture\n' > "$fix/README"
  git -C "$fix" add README
  git -C "$fix" commit -q -m "fixture"
  fixsha=$(git -C "$fix" rev-parse HEAD)
  jobct=zfpr34-ct-$$
  ( env -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE -u POPLOG_STEP \
      -u POPLOG_OUTER_LOG -u POPLOG_OUTER_PID -u POPLOG_OUTER_DEPTH -u POPLOG_OUTER_START -u POPLOG_OUTER_AGENT \
      -u POPSLOT_HELD -u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE \
      "$real_bin/poprun" --agent harper --job "$jobct" --repo "$fix" -- true > "$W/out-ct-clean" 2>&1 )
  ctlog=$(grep -oE 'log [^ ]+\.log' "$W/out-ct-clean" | awk 'END{print $2}')
  if [[ -n $ctlog && -f $ctlog ]] && grep -q '^RESULT: PASS' "$ctlog" && grep -q "gate-commit: $fixsha" "$ctlog"; then
    ok "real poprun clean-tree on a clean fixture: RESULT PASS gate-commit=$fixsha"
  else
    bad "real poprun clean fixture: log=$ctlog out=$(tr '\n' ' ' < "$W/out-ct-clean")"
  fi
  printf 'dirty\n' >> "$fix/README"
  jobdt=zfpr34-dt-$$
  ( env -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE -u POPLOG_STEP \
      -u POPLOG_OUTER_LOG -u POPLOG_OUTER_PID -u POPLOG_OUTER_DEPTH -u POPLOG_OUTER_START -u POPLOG_OUTER_AGENT \
      -u POPSLOT_HELD -u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE \
      "$real_bin/poprun" --agent harper --job "$jobdt" --repo "$fix" -- true > "$W/out-ct-dirty" 2>&1 ) || true
  dtlog=$(grep -oE 'log [^ ]+\.log' "$W/out-ct-dirty" | awk 'END{print $2}')
  if [[ -n $dtlog && -f $dtlog ]] && grep -q 'NOT A GATE (tree dirty' "$dtlog"; then
    ok "real poprun clean-tree on a dirty fixture: NOT A GATE"
  else
    bad "real poprun dirty fixture: log=$dtlog"
  fi
  git -C "$fix" checkout -- README
  # Real popslot: a stale dead holder must not produce SLOT REUSED. This acquires one
  # slot briefly; skip when a slot is already held (this selftest as a gate step) so we
  # do not wait on a fifth Harper slot during a 4-leaf gate.
  if [[ -z ${POPSLOT_HELD:-} ]] && { [[ -z ${POPSLOT_HOLDER_PID:-} ]] || ! kill -0 "$POPSLOT_HOLDER_PID" 2>/dev/null; }; then
    ( env -u POPSLOT_HELD -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_SCOPE -u POPSLOT_HOLDER_PID \
        POPSLOT_HOLDER_PID=999999 POPSLOT_SLOT=3 \
        "$real_bin/popslot" --agent harper -- sh -c 'echo SLOT_NOW=$POPSLOT_SLOT HOLDER_NOW=$POPSLOT_HOLDER_PID' \
        > "$W/out-real-slot" 2>&1 )
    if grep -q 'SLOT REUSED' "$W/out-real-slot"; then
      bad "real popslot reused a stale dead holder: $(tr '\n' ' ' < "$W/out-real-slot")"
    elif grep -q 'SLOT TAKEN' "$W/out-real-slot"; then
      ok "real popslot with a stale dead holder took a fresh slot (not REUSED)"
    else
      bad "real popslot unexpected: $(tr '\n' ' ' < "$W/out-real-slot")"
    fi
  else
    # already inside a slot: confirm the launch line unsets holder vars (static, always)
    ok "skip real popslot acquire (a slot is already held); launch-line unset is checked below"
  fi
else
  ok "skip real poprun/popslot (not on Pop or ZF_GATE_REAL_POPBIN unset)"
  ok "skip real poprun dirty-tree (not on Pop)"
  ok "skip real popslot acquire (not on Pop)"
fi
# always: the leaf launch unsets leftover holder variables
grep -q -- '-u POPSLOT_HOLDER_PID -u POPSLOT_SLOT -u POPSLOT_FD -u POPSLOT_HELD -u POPSLOT_SCOPE' "$gate" \
  && ok "leaf launch unsets POPSLOT_HOLDER_PID/SLOT/FD/HELD/SCOPE" \
  || bad "leaf launch missing POPSLOT unsets"

echo "selftest: $n checks, $fails failed"
[[ $fails == 0 ]]
