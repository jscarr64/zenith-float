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
#   no real gate steps run in the coordinator checks; the serial checks run the real gate's
#   preflight and first step only, in a fixture repo, with contained poplog state and no slot.
#   Run on Pop (the serial checks need the gate's toolchain). Fast (minutes).
set -u
gate=${1:-$(cd "$(dirname "$0")" && pwd)/full_gate.sh}
if [[ -n ${2:-} ]]; then W=$2; mkdir "$W" || { echo "workdir exists"; exit 2; }
else W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; fi
mkdir -p "$W/bin" "$W/tools" "$W/troot"
cat > "$W/bin/poprun" <<'FIXTURE_POPRUN_EOF'
#!/usr/bin/env bash
job=?; while [[ $# -gt 0 && $1 != -- ]]; do [[ $1 == --job ]] && job=$2; shift; done; shift
echo "poprun: log $FAKE_STATE/fake-$job.log" >&2; : > "$FAKE_STATE/fake-$job.log"
exec "$@"
FIXTURE_POPRUN_EOF
cat > "$W/bin/popslot" <<'FIXTURE_POPSLOT_EOF'
#!/usr/bin/env bash
while [[ $# -gt 0 && $1 != -- ]]; do shift; done; shift
[[ -n ${POPSLOT_HOLDER_PID:-} ]] && echo "inherited holder" >> "$FAKE_STATE/holder-leak"
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
SIX='a:target-dir golden-sums arm-tests-release qemu-time-factor build-debug zero-test
b:target-dir build-nostd build-thumbv7em features-random-serde rules-scan no-clone-progress
c:target-dir fmt clippy-x86_64 clippy-i686 audit no-intrinsics hex-32bit unsafe-scan
release:target-dir build-release
mpfr:target-dir mpfr-tests-build oracle-fast oracle-deep mpfr
prepublish:target-dir prepublish'
# case: fresh state dir, run coordinator; extra env as args
coord() { # NAME env...
  local c=$1; shift
  S=$W/state-$c; mkdir -p "$S/active"; : > "$S/targets"; : > "$S/logger-leak"
  ( env FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" ZF_GATE_KILL_GRACE=2 FAKE_SLEEP=2 "$@" bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
}
leftover() { pgrep -f "$W/fakegate" >/dev/null || pgrep -f 'sleep 300' -u "$(id -u)" >/dev/null; }
# 1 all pass, 4 default leaves
coord pass ZF_GATE_PARALLEL=4; rc=$?
[[ $rc == 0 ]] && ok "all leaves pass -> exit 0" || bad "pass exit $rc"
[[ ! -e $W/troot/harper-job-pass ]] && ok "own target root removed after success" || bad "own root not removed"
[[ $(sort "$S/targets" | wc -l) == 4 && $(sort -u "$S/targets" | wc -l) == 4 ]] && ok "4 leaves, 4 distinct targets" || bad "targets: $(tr '\n' ' ' < "$S/targets")"
grep -vq "^$W/troot/harper-job-pass/" "$S/targets" && bad "target outside run root" || ok "every target under the run root"
[[ ! -s $S/logger-leak ]] && ok "no inherited poplog state in leaves" || bad "logger leak: $(cat "$S/logger-leak")"
[[ ! -e $S/holder-leak ]] && ok "no slot holder inherited" || bad "holder leaked into leaf"
[[ $(cat "$S/max") -le 4 ]] && ok "max concurrent $(cat "$S/max") <= 4" || bad "max $(cat "$S/max")"
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
S=$W/state-cancel; mkdir -p "$S/active"; : > "$S/targets"; : > "$S/logger-leak"
env FAKE_STATE="$S" T="$W/T-cancel" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
  ZF_GATE_AGENT=harper ZF_GATE_JOB=job-cancel ZF_GATE_KILL_GRACE=2 ZF_GATE_PARALLEL=2 "ZF_GATE_LEAVES=$SIX" FAKE_a=sleep FAKE_b=sleep \
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
  S=$W/state-$c; mkdir -p "$S/active"; : > "$S/targets"; : > "$S/logger-leak"
  ( env -u ZF_GATE_PARALLEL FAKE_STATE="$S" T="$W/T-$c" REPO="$W" ZF_GATE_POPBIN="$W/bin" ZF_GATE_SELF="$W/fakegate" ZF_GATE_TARGET_ROOT="$W/troot" \
      ZF_GATE_AGENT=harper ZF_GATE_JOB="job-$c" ZF_GATE_KILL_GRACE=2 FAKE_SLEEP=2 "$@" bash "$W/tools/full_gate.sh" "t$c" > "$W/out-$c" 2>&1 )
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
# 11 leaf filter keeps canonical order and only its steps (run() with ZF_GATE_ONLY): exercised by the real-gate run
echo "selftest: $n checks, $fails failed"
[[ $fails == 0 ]]
