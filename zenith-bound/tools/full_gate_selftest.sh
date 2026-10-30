#!/usr/bin/env bash
# Self-test of tools/full_gate.sh: no step may pass VACUOUSLY when the tool it calls is missing (a `! grep ...` step passes if grep is absent),
# and the scan steps must still catch a planted hit and pass a clean tree. Runs the gate on throwaway fixtures under $(zb_scratch) with GATE_STEPS.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap full_gate_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
export LC_ALL=C
root=$(cd "$(dirname "$0")/.." && pwd)
fail() { echo "full_gate selftest FAIL: $*"; exit 1; }
t=$(mktemp -d "$(zb_scratch)"/zb-gate-selftest.XXXXXX); runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
SCANS="unsafe-scan no-float-scan no-panic-scan no-clone-scan derive-scan manifest-scan"
# The fixture gates run outside this self-test's own run log (their planted failures must not fail it), with no poplog at all.
fe=(env -u GATE_STEPS -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE RUNLOG_POPLOG=/nonexistent RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1)
rcof() { sed -n "s/^$2=//p" "$t/logs/$1.rc"; }

# 1. stripped PATH: only what the gate driver itself needs; every tool a step calls (cargo, grep, python3, ...) is missing
mkdir -p "$t/bin" "$t/logs"
ck; for c in bash timeout mkdir date; do ln -s "$(command -v "$c")" "$t/bin/$c" || fail "cannot link $c"; done
ck; "${fe[@]}" T="$t/logs" PATH="$t/bin" "$t/bin/bash" "$root/tools/full_gate.sh" stripped > "$t/stripped.out" 2>&1 && fail "stripped run (every step failing) exited 0"
ck; grep -q '^failed$' "$t/logs/stripped.rc" || fail "stripped run did not finish with 'failed': $(cat "$t/stripped.out")"
ck; grep -q '^done$' "$t/logs/stripped.rc" && fail "stripped run wrote 'done' although steps failed"
steps=$(sed -n 's/^\([a-z0-9-]*\)=[0-9]*$/\1/p' "$t/logs/stripped.rc")
ck; [ "$(echo "$steps" | wc -w)" -eq 25 ] || fail "stripped run reported $(echo "$steps" | wc -w) steps, want exactly 25 (bump this with every tools/full_gate.sh step): $steps"
for need in gen-la-guard gen-la-heavy la4-red-markers; do ck; echo "$steps" | grep -qx "$need" || fail "gate step $need missing from tools/full_gate.sh"; done
for s in $steps; do
  ck; [ "$(rcof stripped "$s")" != 0 ] || fail "step $s PASSED with its tools missing (vacuous pass)"
done
ck; for s in $SCANS; do grep -q 'MISSING TOOL' "$t/logs/stripped-$s.log" || fail "step $s does not name the missing tool"; done

# 2. full PATH, fixture trees: a planted hit fails each scan, a clean tree passes, a missing scanned dir is not an error
mk() { mkdir -p "$t/$1/src" "$t/$1/tests"; printf 'pub fn f(x: u8) -> u8 {\n    x\n}\n' > "$t/$1/src/lib.rs"; }
mk clean
ck; "${fe[@]}" GATE_STEPS="$SCANS" T="$t/logs" REPO="$t/clean" "$root/tools/full_gate.sh" clean > /dev/null 2>&1 || fail "clean gate run exited non-zero"
ck; grep -q '^done$' "$t/logs/clean.rc" || fail "clean gate run did not write 'done'"
ck; for s in $SCANS; do [ "$(rcof clean "$s")" = 0 ] || fail "step $s fails on a clean tree: $(cat "$t/logs/clean-$s.log")"; done
plant() { # plant <step> <line>
  mk "p-$1"; printf '%s\n' "$2" >> "$t/p-$1/src/lib.rs"
  ck; "${fe[@]}" GATE_STEPS="$1" T="$t/logs" REPO="$t/p-$1" "$root/tools/full_gate.sh" "p-$1" > /dev/null 2>&1 && fail "gate with failed step $1 exited 0 (silent pass)"
  ck; grep -q '^done$' "$t/logs/p-$1.rc" && fail "gate with failed step $1 wrote 'done'"
  ck; [ "$(rcof "p-$1" "$1")" = 1 ] || fail "step $1 missed a planted hit ($(rcof "p-$1" "$1"))"
  ck; grep -q 'src/lib.rs' "$t/logs/p-$1-$1.log" || fail "step $1 does not show the hit"
}
plant unsafe-scan 'pub fn g() { unsafe { } }'
plant no-float-scan 'pub fn g(x: f64) -> f64 { x }'
plant no-panic-scan 'pub fn g(x: Option<u8>) -> u8 { x.unwrap() }'
plant no-clone-scan 'pub fn g(x: &[u8]) -> Vec<u8> { x.to_vec() }'
plant derive-scan '#[derive(Clone)]'
plant manifest-scan 'pub fn g() -> &'"'"'static str { env!("CARGO_MANIFEST_DIR") }'
# GATE_STEPS is honoured: only the named step ran
ck; [ "$(grep -c '^[a-z-]*=' "$t/logs/p-unsafe-scan.rc")" = 1 ] || fail "GATE_STEPS ran other steps"
# missing poplog: a hard failure (exit non-zero, nothing runs), not a warning; only RUNLOG_PLAIN=1 (self-test fixtures) runs without it
env -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE -u RUNLOG_PLAIN RUNLOG_POPLOG=/nonexistent GATE_STEPS=unsafe-scan T="$t/logs" REPO="$t/clean" \
  "$root/tools/full_gate.sh" nopoplog > "$t/nopoplog.out" 2>&1 && fail "gate ran with poplog.sh missing: $(cat "$t/nopoplog.out")"
ck; grep -q 'missing' "$t/nopoplog.out" || fail "missing poplog not named: $(cat "$t/nopoplog.out")"
ck; grep -q '^unsafe-scan=' "$t/logs/nopoplog.rc" 2>/dev/null && fail "a step ran with poplog.sh missing"
# RUNLOG_PLAIN without the selftest-only marker RUNLOG_SELFTEST=1 is refused (exit non-zero, nothing runs): plain mode is for fixtures only
env -u POPLOG_FILE -u POPLOG_STATE_DIR -u POPLOG_PID -u POPLOG_PROGRESS_FILE -u RUNLOG_SELFTEST RUNLOG_PLAIN=1 GATE_STEPS=unsafe-scan T="$t/logs" REPO="$t/clean" \
  "$root/tools/full_gate.sh" plainnomark > "$t/plainnomark.out" 2>&1 && fail "gate ran with RUNLOG_PLAIN=1 but no RUNLOG_SELFTEST=1: $(cat "$t/plainnomark.out")"
ck; grep -q 'RUNLOG_SELFTEST' "$t/plainnomark.out" || fail "refused plain run does not name RUNLOG_SELFTEST: $(cat "$t/plainnomark.out")"
ck; grep -q '^unsafe-scan=' "$t/logs/plainnomark.rc" 2>/dev/null && fail "a step ran with RUNLOG_PLAIN=1 but no marker"
# 3. toolchain step: a rustc that is not the pinned channel (Pop's /usr/bin/rustc 1.95 ahead of ~/.cargo/bin) fails the gate
mk tc; mkdir -p "$t/tc/tools"; cp "$root/tools/toolchain_check.sh" "$root/tools/toolchain_check_selftest.sh" "$t/tc/tools/"; printf '[toolchain]\nchannel = "1.99"\n' > "$t/tc/rust-toolchain.toml"
for v in 1.95.0 1.99.1; do
  mkdir -p "$t/fake-$v"; for c in rustc cargo; do printf '#!/bin/sh\necho "%s %s (x 2026-01-01)"\n' "$c" "$v" > "$t/fake-$v/$c"; chmod +x "$t/fake-$v/$c"; done
  "${fe[@]}" GATE_STEPS=toolchain T="$t/logs" REPO="$t/tc" PATH="$t/fake-$v:$PATH" "$root/tools/full_gate.sh" "tc-$v" > /dev/null 2>&1
done
ck; [ "$(rcof tc-1.95.0 toolchain)" = 1 ] || fail "toolchain step passed with rustc 1.95.0 for channel 1.99 ($(rcof tc-1.95.0 toolchain))"
ck; grep -q 'TOOLCHAIN MISMATCH' "$t/logs/tc-1.95.0-toolchain.log" || fail "toolchain step does not name the mismatch"
ck; [ "$(rcof tc-1.99.1 toolchain)" = 0 ] || fail "toolchain step fails on the pinned rustc: $(cat "$t/logs/tc-1.99.1-toolchain.log")"
# poplog treats any step whose command says `cargo test` as a test step, so a build-only `cargo test --no-run` step reads as zero tests -> FAIL
ck; grep -nE '^run [a-z-]+ .*cargo test[^"]*--no-run' "${BASH_SOURCE[0]%/*}/full_gate.sh" && fail "a gate step builds with cargo test --no-run; poplog counts it as zero tests (use cargo build --tests)"
# 4. mutants step: MUTATE_TARGET (caller-set or default) gets the same ownership claim as CARGO_TARGET_DIR (Spence, 2026-10-23); the fixture's
# mutate.sh/mutate_selftest.sh are stubs that pass, so only the claim can fail the step
mk mu; mkdir -p "$t/mu/tools" "$t/mu-bin"; cp "$root/tools/fresh_test_run.sh" "$root/tools/runlog.sh" "$root/tools/scratch.sh" "$t/mu/tools/"
printf '#!/bin/sh\necho stub mutate_selftest.sh\n' > "$t/mu/tools/mutate_selftest.sh"; chmod +x "$t/mu/tools/mutate_selftest.sh"
# fake tools/mutate.sh: a runner that ALWAYS exits 0 (the fail-open case the gate must not trust). Slices: central (001-a 002-b), s1 (201-x),
# s2 (empty). FAKE_MUT picks what the run prints: ok | timeout | missing | invalid | survived | floor (floor: central needs 3 kills) |
# slow (002-b takes 3s, then KILLED) | dup (002-b KILLED twice) | unknown (the 002-b run reports 999-z) | forged (002-b prints only an
# indented 'KILLED 002-b', as mutate.sh echoes suite output) | rc (exits 3) | oom (002-b SIGKILLs itself, exit 137) | floorbad | listbad |
# slicesbad (--slices exits 1 after two stderr lines).
cat > "$t/mu/tools/mutate.sh" <<'FAKE'
#!/usr/bin/env bash
# the real mutate.sh reads these (a stale one points it at another manifest/floor/repo); the gate must never pass them on
for v in MUTATE_REPO MUTATE_MANIFEST MUTATE_EQUIVALENT MUTATE_FLOOR MUTATE_TIMEOUT MUTATE_BUILD_TIMEOUT MUTATE_TMP MUTATE_PERSIST; do [ -z "${!v+x}" ] || { echo "fake mutate.sh: leaked $v=${!v}" >&2; exit 1; }; done
sc=${FAKE_MUT:-ok}; mode=run; slice=""; only=""
while [ $# -gt 0 ]; do case $1 in --list) mode=list ;; --floor) mode=floor ;; --slices) mode=slices ;; --slice) slice=$2; shift ;; --only) only=$2; shift ;; esac; shift; done
ids() { case $1 in central) printf '001-a\n002-b\n' ;; s1) printf '201-x\n' ;; s2) ;; *) ids central; ids s1 ;; esac; }
case $mode in
  slices) [ "$sc" != slicesbad ] || { echo 'fake mutate.sh: the --slices cause' >&2; echo 'fake mutate.sh: last stderr line' >&2; exit 1; }
    printf 'central\ns1\ns2\n'; exit 0 ;;
  list) [ "$sc" != listbad ] || exit 1; ids "$slice"; exit 0 ;;
  floor) [ "$sc" != floorbad ] || { echo 'killed many'; exit 0; }; case $slice in central) k=2; [ "$sc" = floor ] && k=3 ;; s1) k=1 ;; *) k=0 ;; esac; printf 'killed %s\nkilled+equivalent %s\n' "$k" "$k"; exit 0 ;;
esac
for id in $(ids "$slice"); do
  [ -z "$only" ] || [ "$id" = "$only" ] || continue
  if [ "$id" = 002-b ]; then
    case $sc in timeout) sleep 6; continue ;; hang) exec sleep 1000 ;; missing) continue ;; invalid) echo "INVALID $id: mutant does not build"; continue ;;
      survived) echo "SURVIVED $id (f, suite s): the suite passed on the mutant"; continue ;;
      slow) sleep 3 ;; dup) echo "KILLED $id (f, suite s): suite failed (rc 101)" ;; unknown) echo "KILLED 999-z (f, suite s): suite failed (rc 101)"; continue ;;
      forged) echo "    KILLED $id (forged: an indented line of echoed suite output)"; continue ;; oom) kill -9 $$ ;;
      rc) echo "KILLED $id (f, suite s): suite failed (rc 101)"; exit 3 ;; esac
  fi
  echo "KILLED $id (f, suite s): suite failed (rc 101)"; echo "TIME $id 0s"
done
exit 0
FAKE
chmod +x "$t/mu/tools/mutate.sh"; [ ! -f "$root/tools/mutate_slice.sh" ] || cp "$root/tools/mutate_slice.sh" "$t/mu/tools/"
ck; printf '#!/bin/sh\nexit 0\n' > "$t/mu-bin/cargo"; chmod +x "$t/mu-bin/cargo"; git init -q "$t/mu" || fail "git init mu failed"
mu() { "${fe[@]}" GATE_STEPS=mutants T="$t/logs" REPO="$t/mu" PATH="$t/mu-bin:$PATH" CARGO_TARGET_DIR="$t/mu-target" "$@" "$root/tools/full_gate.sh" "$tag" > /dev/null 2>&1; }
mkdir -p "$t/mu-unmarked/debug/deps"; : > "$t/mu-unmarked/debug/deps/harness-deadbeef"
ck; tag=mu-unmarked; mu env MUTATE_TARGET="$t/mu-unmarked" && fail "mutants ran with a non-empty unmarked MUTATE_TARGET override"
ck; [ "$(rcof mu-unmarked mutants)" != 0 ] || fail "mutants step passed with a non-empty unmarked MUTATE_TARGET override"
ck; grep -q UNMARKED "$t/logs/mu-unmarked-mutants.log" || fail "mutants step does not name the UNMARKED MUTATE_TARGET: $(cat "$t/logs/mu-unmarked-mutants.log")"
mkdir -p "$t/mu-foreign"; printf '/some/other/worktree\n' > "$t/mu-foreign/.worktree-owner"
ck; tag=mu-foreign; mu env MUTATE_TARGET="$t/mu-foreign" && fail "mutants ran with a MUTATE_TARGET owned by another worktree"
ck; grep -q SHARED "$t/logs/mu-foreign-mutants.log" || fail "mutants step does not name the SHARED MUTATE_TARGET: $(cat "$t/logs/mu-foreign-mutants.log")"
ck; tag=mu-default; mu env || fail "mutants step failed on a fresh default MUTATE_TARGET: $(cat "$t/logs/mu-default-mutants.log")"
ck; [ "$(cat "$t/mu-target/mutate/.worktree-owner" 2> /dev/null)" = "$(cd "$t/mu" && pwd -P)" ] || fail "default MUTATE_TARGET got no owner marker"
# 5. mutants per slice through tools/mutate_slice.sh (Spence-approved design, 2026-10-24): the gate must not trust the runner's exit code.
# Each case runs the fail-open fake above; the mutants step must FAIL and say why. Failures are collected so a red run shows every miss.
miss=""
mcase() { # mcase SCENARIO ERE
  ck; tag=mu-$1; if mu env FAKE_MUT="$1"; then miss+="[$1: gate exited 0] "; fi
  ck; grep -qE -- "$2" "$t/logs/mu-$1-mutants.log" 2>/dev/null || miss+="[$1: mutants log lacks /$2/] "
}
mcase missing 'INCOMPLETE: 1 results for 2 ids.*missing: 002-b'
mcase invalid 'FAIL INVALID: 002-b'
mcase survived 'FAIL SURVIVED \(no equivalent entry\): 002-b'
mcase floor 'FAIL BELOW floor'
# 6. per-mutant timeout (60s, measured max 14s x ~4): a runner that hangs on one id is caught in about 60s with `TIMEOUT <id>`, not after
# the slice budget (420s here: 2 ids x 60s + 300s setup; hours on the real manifest), so only a per-mutant limit can catch it under 90s.
# The caller's env leaks MUTANT_TIMEOUT_S=999999 and a huge budget/setup (Spence, 1742-zb-fgself2-9ea44e6): the gate must not inherit them.
# The outer 150s timeout only bounds a red run.
ck; tag=mu-hang; h0=$SECONDS
if mu timeout -k 5 150 env FAKE_MUT=hang MUTANT_TIMEOUT_S=999999 MUTANT_BUDGET_S=999999 MUTATE_SETUP_S=999999; then miss+="[hang: gate exited 0] "; fi; hel=$((SECONDS - h0))
ck; { grep -qE '(^|: )TIMEOUT 002-b' "$t/logs/mu-hang-mutants.log" 2>/dev/null && [ "$hel" -lt 90 ]; } || miss+="[hang: leaked MUTANT_TIMEOUT_S/budget: no 'TIMEOUT 002-b' within 90s (took ${hel}s)] "
# env leakage: the gate's mutants step (its --step-timeout and its --all run) must ignore the caller's MUTANT_TIMEOUT_S, MUTANT_BUDGET_S,
# MUTATE_SETUP_S and MUTATE_SH. Junk values (or a missing runner) would fail it; a 1s budget with no setup would cut the 3s slow mutant.
ck; tag=mu-leakjunk; mu env FAKE_MUT=ok MUTANT_TIMEOUT_S=x MUTANT_BUDGET_S=x MUTATE_SETUP_S=x MUTATE_SH="$t/no-such-mutate.sh" || miss+="[leak: junk MUTANT_*/MUTATE_SETUP_S/MUTATE_SH reached the gate: $(tail -2 "$t/logs/mu-leakjunk-mutants.log" 2>/dev/null)] "
ck; tag=mu-leakbudget; mu env FAKE_MUT=slow MUTANT_BUDGET_S=1 MUTATE_SETUP_S=0 || miss+="[leak: MUTANT_BUDGET_S=1 MUTATE_SETUP_S=0 reached the gate: $(tail -2 "$t/logs/mu-leakbudget-mutants.log" 2>/dev/null)] "
# a stale MUTATE_FLOOR / MUTATE_EQUIVALENT / MUTATE_MANIFEST (or any other env mutate.sh reads) must not reach the gate's mutate.sh
ck; tag=mu-leakmut; mu env FAKE_MUT=ok MUTATE_REPO="$t/stale" MUTATE_MANIFEST="$t/stale/manifest.txt" MUTATE_EQUIVALENT="$t/stale/EQUIVALENT.txt" \
  MUTATE_FLOOR="$t/stale/floor.txt" MUTATE_TIMEOUT=1s MUTATE_BUILD_TIMEOUT=1s MUTATE_TMP="$t/stale" MUTATE_PERSIST="$t/stale" \
  || miss+="[leak: stale MUTATE_* reached the gate's mutate.sh: $(tail -3 "$t/logs/mu-leakmut-mutants.log" 2>/dev/null)] "
# 7. ID_REGEX: `mutate_slice.sh SLICE ID_REGEX` runs only matching ids, skips floor/completeness and says so, never prints a full pass
ms() { "${fe[@]}" PATH="$t/mu-bin:$PATH" "$@"; }
ck; ms env FAKE_MUT=floor "$t/mu/tools/mutate_slice.sh" central '^001' > "$t/rx.out" 2>&1 || miss+="[regex: partial run with a floor miss exited nonzero] "
ck; grep -qF 'PARTIAL run (ID_REGEX=^001): floor/completeness not checked' "$t/rx.out" || miss+="[regex: no PARTIAL line] "
ck; { grep -q 'KILLED 001-a' "$t/rx.out" && ! grep -q '002-b' "$t/rx.out"; } || miss+="[regex: did not run only 001-a] "
ck; grep -qE -- '-> ok$|all slices ok' "$t/rx.out" && miss+="[regex: printed a full-pass line] "
ck; ms "$t/mu/tools/mutate_slice.sh" central 'nomatch' > "$t/rx0.out" 2>&1 && miss+="[regex: a regex matching nothing exited 0] "
ck; grep -q 'matches no id' "$t/rx0.out" || miss+="[regex: no-match not named] "
ck; ms "$t/mu/tools/mutate_slice.sh" --all '^001' > "$t/rxa.out" 2>&1 && miss+="[regex: --all accepted a regex] "
ck; grep -E '^run mutants ' "$root/tools/full_gate.sh" | grep -qE 'mutate_slice\.sh --all"$' || miss+="[full_gate passes something after --all] "
# 8. mutate_slice.sh on its own (the selftest-only budget env is honoured there, never through the gate): a slice past its budget
# (2 ids x 1s + 0s setup; 002-b sleeps 6s) is TIMEOUT and names the id that never ran
ck; ms env FAKE_MUT=timeout MUTANT_BUDGET_S=1 MUTATE_SETUP_S=0 "$t/mu/tools/mutate_slice.sh" central > "$t/to.out" 2>&1 && miss+="[timeout: exited 0] "
ck; grep -qE 'TIMEOUT.*never ran: 002-b' "$t/to.out" || miss+="[timeout: no 'TIMEOUT ... never ran: 002-b': $(tail -2 "$t/to.out")] "
# 9. every other failure the mutate_slice.sh header promises, one fake scenario each: the slice FAILS and names the cause
dcase() { # dcase SCENARIO ERE
  ck; ms env FAKE_MUT="$1" "$t/mu/tools/mutate_slice.sh" central > "$t/d-$1.out" 2>&1 && miss+="[$1: mutate_slice exited 0] "
  ck; grep -qE -- "$2" "$t/d-$1.out" || miss+="[$1: output lacks /$2/: $(tail -2 "$t/d-$1.out")] "
}
dcase dup 'FAIL slice central: duplicate result for 002-b'
dcase unknown 'FAIL slice central: the --only 002-b run printed a result for 999-z'
dcase floorbad 'FAIL slice central: unreadable floor: killed many'
dcase listbad 'FAIL mutate.sh --list --slice central failed'
dcase rc 'FAIL slice central: mutate.sh exited 3'
# a forged indented 'KILLED 002-b' (suite output that mutate.sh echoes indented) is not a result: 002-b has none
dcase forged 'FAIL INCOMPLETE: 1 results for 2 ids in slice central; missing: 002-b'
# a runner SIGKILLed well under its limit (exit 137, the OOM killer) is KILLED-OOM, not TIMEOUT, and still fails the slice
dcase oom 'KILLED-OOM 002-b'
ck; grep -q 'TIMEOUT' "$t/d-oom.out" && miss+="[oom: reported as TIMEOUT: $(grep TIMEOUT "$t/d-oom.out")] "
# 10. a --step-timeout failure keeps ALL its output in the step log, not only the last line: 1706-zb-mslice-9ea44e6-selftests failed with
# 'mutate.sh --slices failed' once and its cause was lost to `| tail -1` (Spence nit b)
ck; tag=mu-slicesbad; mu env FAKE_MUT=slicesbad && miss+="[slicesbad: gate exited 0] "
ck; grep -q 'the --slices cause' "$t/logs/mu-slicesbad-mutants.log" 2>/dev/null || miss+="[slicesbad: step log lost the cause: $(cat "$t/logs/mu-slicesbad-mutants.log" 2>/dev/null)] "
[ -z "$miss" ] || fail "mutants step did not catch: $miss"
st_done full_gate 118   # WANT: cases this script runs (tools/selftest_cases.sh)
# the gate's real mutation run must not inherit MUTATE_FLOOR_ZERO_OK=1 (only the selftest fixture may use a 0 floor), nor the selftest-only
# MUTANT_TIMEOUT_S / MUTANT_BUDGET_S / MUTATE_SETUP_S / MUTATE_SH, nor any env mutate.sh reads (MUTATE_REPO MUTATE_MANIFEST MUTATE_EQUIVALENT MUTATE_FLOOR MUTATE_TIMEOUT MUTATE_BUILD_TIMEOUT MUTATE_TMP MUTATE_PERSIST), in its --all run or in its --step-timeout computation
u='env -u MUTATE_FLOOR_ZERO_OK -u MUTANT_TIMEOUT_S -u MUTANT_BUDGET_S -u MUTATE_SETUP_S -u MUTATE_SH -u MUTATE_REPO -u MUTATE_MANIFEST -u MUTATE_EQUIVALENT -u MUTATE_FLOOR -u MUTATE_TIMEOUT -u MUTATE_BUILD_TIMEOUT -u MUTATE_TMP -u MUTATE_PERSIST tools/mutate_slice\.sh'
grep -E '^run mutants ' "$root/tools/full_gate.sh" | grep -q "$u --all\"$" || fail "mutants step runs tools/mutate_slice.sh --all without $u"
[ "$(grep -c "$u --step-timeout" "$root/tools/full_gate.sh")" = 1 ] || fail "the mutants step timeout is computed without $u --step-timeout"
echo "full_gate selftest: ok ($(grep -c x "$ST_CASES") cases)"
