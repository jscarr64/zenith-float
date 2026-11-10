#!/usr/bin/env bash
# Self-test for the Pop direct-run guard (~/pop-logs/bin/popguard.sh; 11:49 2026-10-24 incident: a bare 58 GB goldens python run
# global-OOMed Pop). Heavy golden generation must refuse (exit 3, clear message) outside a popslot/poprun memory scope:
# tools/goldens_regen.sh before it starts any generator, and scripts/goldlib.py (imported by every mpmath generator) at import.
# Proves: no POPSLOT_SCOPE -> refused, nothing ran; a forged POPSLOT_SCOPE (not our cgroup) -> refused; inside a real
# popslot-*.scope -> runs; an uncapped scope -> refused; POPGUARD pointing nowhere on Pop -> refused exit 2; no popguard.sh
# (CI, other hosts) -> loud WARNING, runs unguarded; the selftest itself never skips on Pop. Uses its own throwaway systemd --user scope, no gate slot.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap popguard_selftest "$@"   # run log (tools/runlog.sh)
set -u
unset ZB_GATE_MODE GOLDENS_REGEN_ARGS   # fixtures run regular mode even inside a merge-gate run (gate 1707: inherited ZB_GATE_MODE=merge refused every fixture)
here="$(cd "$(dirname "$0")" && pwd)"; repo="$(cd "$here/.." && pwd)"
t="$(mktemp -d "$(zb_scratch)/popguard-selftest.XXXXXX")"; runlog_on_exit "rm -rf \"$t\""
. "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
# Skip ONLY on a host that is not Pop (no ~/pop-logs/bin). On Pop (inside full_gate or not) a missing or non-executable guard, or a
# POPGUARD pointing nowhere, FAILS: a silent skip there would hide a disabled guard (Spence, 2026-10-24 block of 7881073).
if [ -n "${POPGUARD:-}" ]; then PG=$POPGUARD; [ -x "$PG" ] || fail "POPGUARD=$PG is missing or not executable (refusing to skip)"
elif [ -x "$HOME/pop-logs/bin/popguard.sh" ]; then PG=$HOME/pop-logs/bin/popguard.sh
elif [ -d "$HOME/pop-logs/bin" ]; then fail "this is Pop ($HOME/pop-logs/bin) but $HOME/pop-logs/bin/popguard.sh is missing or not executable (refusing to skip)"
else echo "popguard selftest: skipped (no $HOME/pop-logs/bin: not Pop)"; exit 0; fi
export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
# fixture: one slice whose generator is a stub interpreter that only leaves a marker
root="$t/r"; d="$root/tests/goldens/slices/x"; mkdir -p "$d"; printf 'gen_x.py\n' > "$d/GENERATORS"; printf '# EMPTY\n' > "$d/SHA256SUMS"
stub="$t/stubpy"; printf '#!/bin/sh\ntouch "%s/ran"\n' "$t" > "$stub"; chmod +x "$stub"
gr() { GOLDENS_ROOT="$root" GOLDENS_PYTHON="$stub" RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1 "$@" "$here/goldens_regen.sh" --slices-only --regen 2>&1; }   # --regen: python runs every call (an unchanged key would skip it, mode=cached)
inscope() { local u="popslot-zbselftest-$$-$RANDOM.scope"; POPSLOT_SCOPE=$u systemd-run --user --scope -q --collect "--unit=$u" -p MemoryMax=256M -p MemorySwapMax=0 -- "$@"; }
ck; rm -f "$t/ran"; out=$(gr env -u POPSLOT_SCOPE); rc=$?; [ "$rc" = 3 ] || fail "bare goldens_regen: want exit 3, got $rc: $out"
ck; printf '%s\n' "$out" | grep -q "REFUSING.*goldens_regen.*popslot" || fail "bare goldens_regen: unclear refusal: $out"
ck; [ ! -e "$t/ran" ] || fail "bare goldens_regen still started a generator"
ck; out=$(gr env POPSLOT_SCOPE=popslot-forged.scope); rc=$?; [ "$rc" = 3 ] || fail "forged POPSLOT_SCOPE accepted (rc $rc): $out"
ck; rm -f "$t/ran"; out=$(gr inscope env); rc=$?; [ "$rc" = 0 ] && [ -e "$t/ran" ] || fail "goldens_regen inside a popslot scope did not run (rc $rc): $out"
# a scope named popslot-*.scope but without a memory cap is refused (popguard.sh checks the scope's memory.max)
ck; rm -f "$t/ran"; u="popslot-zbnocap-$$-$RANDOM.scope"; out=$(gr env POPSLOT_SCOPE=$u systemd-run --user --scope -q --collect "--unit=$u" -p MemorySwapMax=0 -- env); rc=$?
  [ "$rc" = 3 ] && [ ! -e "$t/ran" ] || fail "uncapped popslot-*.scope accepted (rc $rc): $out"
# POPGUARD pointing nowhere must not switch the guard off (Spence, 2026-10-24): sandbox HOMEs simulate the two hosts.
hpg="$t/hpg"; mkdir -p "$hpg/pop-logs/bin"; cp "$PG" "$hpg/pop-logs/bin/popguard.sh"; chmod +x "$hpg/pop-logs/bin/popguard.sh"   # Pop
hno="$t/hno"; mkdir -p "$hno/.local"                                                                                                  # CI/other host
# goldlib resolves its team env from HOME, so a bare sandbox HOME has none and exits 1 at the oracle check (not the guard).
# Stand-in: the sandbox ~/team-env is a link to the real ~/team-env (team env oracle-2026-10-30; never ~/.local/oracle-venv), so the no-guard case's exit code is deterministic (0).
[ -x "$HOME/team-env/oracle-2026-10-30/bin/python" ] || fail "team env $HOME/team-env/oracle-2026-10-30 missing: cannot build the no-guard host stand-in"
ln -s "$HOME/team-env" "$hno/team-env"
sbx() { HOME="$1" ZB_SCRATCH="$t/scr" "${@:2}"; }; mkdir -p "$t/scr"
ck; rm -f "$t/ran"; out=$(gr sbx "$hpg" env -u POPSLOT_SCOPE POPGUARD="$t/none"); rc=$?
  [ "$rc" = 2 ] && [ ! -e "$t/ran" ] && printf '%s\n' "$out" | grep -q "REFUSING.*POPGUARD" || fail "Pop host, POPGUARD=/nonexistent: goldens_regen must refuse exit 2 (rc $rc): $out"
ck; rm -f "$t/ran"; out=$(gr sbx "$hno" env -u POPSLOT_SCOPE -u POPGUARD); rc=$?
  [ "$rc" = 0 ] && [ -e "$t/ran" ] && printf '%s\n' "$out" | grep -q "WARNING.*popguard" || fail "no popguard.sh host: goldens_regen must WARN and run (rc $rc): $out"
ck; out=$(cd "$repo" && sbx "$hpg" env -u POPSLOT_SCOPE POPGUARD="$t/none" PYTHONNOUSERSITE=1 /usr/bin/python3 scripts/goldlib.py --selftest 2>&1); rc=$?
  [ "$rc" = 2 ] && printf '%s\n' "$out" | grep -q "REFUSING.*POPGUARD" || fail "Pop host, POPGUARD=/nonexistent: goldlib must refuse exit 2 (rc $rc): $out"
ck; out=$(cd "$repo" && sbx "$hno" env -u POPSLOT_SCOPE -u POPGUARD PYTHONNOUSERSITE=1 "$hno/team-env/oracle-2026-10-30/bin/python" scripts/goldlib.py --selftest 2>&1); rc=$?
  [ "$rc" = 0 ] && printf '%s\n' "$out" | grep -q "^goldlib: WARNING: no popguard.sh on this host ($hno/pop-logs/bin/popguard.sh): running WITHOUT the Pop memory guard$" \
    && printf '%s\n' "$out" | grep -q "^goldlib selftest: ok" || fail "no popguard.sh host: goldlib must WARN and pass its selftest, exit 0 (rc $rc): $out"
# the selftest itself: on Pop it never skips (a skip inside full_gate would hide a disabled guard); off Pop it skips
st() { ( unset POPLOG_FILE POPLOG_PID POPLOG_STATE_DIR POPLOG_PROGRESS_FILE POPLOG_RUNNING RUNLOG_INNER RUNLOG_ACTIVE
         RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1 "$@" "$here/popguard_selftest.sh" 2>&1 ); }
mkdir -p "$t/hpop/pop-logs/bin"   # Pop, but popguard.sh gone
ck; out=$(st sbx "$t/hpop" env -u POPGUARD); rc=$?; [ "$rc" != 0 ] && printf '%s\n' "$out" | grep -q "SELFTEST FAIL.*refusing to skip" || fail "selftest on Pop without popguard.sh skipped with rc $rc: $out"
ck; out=$(st env POPGUARD="$t/none"); rc=$?; [ "$rc" != 0 ] || fail "selftest with POPGUARD=/nonexistent skipped with rc 0: $out"
ck; out=$(st sbx "$hno" env -u POPGUARD); rc=$?; [ "$rc" = 0 ] && printf '%s\n' "$out" | grep -q skipped || fail "selftest off Pop should skip (rc $rc): $out"
ck; out=$(cd "$repo" && env -u POPSLOT_SCOPE PYTHONNOUSERSITE=1 /usr/bin/python3 scripts/goldlib.py --selftest 2>&1); rc=$?
  [ "$rc" = 3 ] && printf '%s\n' "$out" | grep -q REFUSING || fail "bare goldlib import: want popguard refusal exit 3, got $rc: $out"
ck; out=$(cd "$repo" && inscope env PYTHONNOUSERSITE=1 /usr/bin/python3 scripts/goldlib.py --selftest 2>&1); rc=$?
  printf '%s\n' "$out" | grep -q REFUSING && fail "goldlib refused inside a popslot scope: $out"
echo "popguard selftest: ok"
st_done popguard 15   # WANT: cases this script runs (tools/selftest_cases.sh)
