#!/usr/bin/env bash
# ZF-M17: full_gate.sh logs as the calling agent, never a hard-coded one. The agent is POPLOG_AGENT, else
# POPSLOT_AGENT (set by popslot --agent), else the gate refuses to run (exit 2). The default log dir T is
# $HOME/<agent>-work/logs/zf-gate/... Runs with stub cargo, HOME and POPLOG_ROOT in a selftest scratch dir.
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY POPLOG_AGENT \
  T ZF_GATE_LOG_DIR CARGO_BUILD_JOBS RUST_TEST_THREADS
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$("$root/tools/selftest_scratch.sh" test_gate_agent) || exit 2; trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin"
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/cargo"; cp "$tmp/bin/cargo" "$tmp/bin/rustup"; chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
fail=0
# case NAME POPLOG_AGENT POPSLOT_AGENT WANT_AGENT ("" = must refuse with exit 2)
case_agent() {
  local name=$1 pl=$2 ps=$3 want=$4 home="$tmp/$1" rc=0
  mkdir -p "$home"
  env ${pl:+POPLOG_AGENT=$pl} ${ps:+POPSLOT_AGENT=$ps} HOME="$home" PATH="$tmp/bin:$PATH" POPLOG_ROOT="$home/pop" \
    POPSLOT_PREFIX="$home/slot" CARGO_TARGET_DIR="$home/target" GATE_STEPS=fmt STEP_TIMEOUT=60s \
    bash "$root/tools/full_gate.sh" agenttest > "$home/out" 2>&1 || rc=$?
  if [[ -z $want ]]; then
    [[ $rc -eq 2 ]] || { echo "FAIL[$name]: no agent set: rc=$rc, want 2"; fail=1; }
    [[ -z $(find "$home" -mindepth 1 -maxdepth 1 -name '*-work') ]] || { echo "FAIL[$name]: no agent set, but made $(ls -d "$home"/*-work)"; fail=1; }
    [[ ! -e "$home/pop" ]] || { echo "FAIL[$name]: no agent set, but wrote a poplog"; fail=1; }
    return 0
  fi
  [[ -n $(find "$home/pop/$want" -name '*full_gate-agenttest*.log' 2>/dev/null) ]] || { echo "FAIL[$name]: no poplog under pop/$want: $(find "$home/pop" -name '*.log' 2>/dev/null | tr '\n' ' ')"; fail=1; }
  [[ -s $(find "$home/$want-work/logs/zf-gate" -name agenttest.rc 2>/dev/null | head -1) ]] || { echo "FAIL[$name]: default T not under \$HOME/$want-work/logs/zf-gate: $(find "$home" -name agenttest.rc | tr '\n' ' ')"; fail=1; }
  return 0
}
case_agent A-poplog pa "" pa
case_agent B-popslot "" ps ps
case_agent C-both pa ps pa
case_agent D-none "" "" ""
# E: no hard-coded agent name left in the gate
! grep -n 'rusty' "$root/tools/full_gate.sh" || { echo "FAIL[E]: full_gate.sh still names an agent"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: full_gate.sh agent = POPLOG_AGENT, else POPSLOT_AGENT, else exit 2; T under \$HOME/<agent>-work"
exit "$fail"
