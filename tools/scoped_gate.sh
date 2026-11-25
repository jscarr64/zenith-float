#!/usr/bin/env bash
# scoped_gate.sh: a scoped gate run from the repo itself (Cassie round-2 note (c); Jeff Q6 2026-11-25: a scoped gate
# per PR, the full gate at tier boundaries), so the gate is reproducible from the PR head.
# usage: tools/scoped_gate.sh SCOPE          run the scope's steps
#        tools/scoped_gate.sh --list SCOPE   print the step names
#        tools/scoped_gate.sh --dry-run SCOPE print "name: command" per step, run nothing
# Scopes:
#   rules  for a PR that changes only the rule-check tools (tools/*, .githooks/*, no Rust):
#          full-gate-syntax  bash -n tools/full_gate.sh
#          rules-scan        \ the command of `run NAME "CMD"` in tools/full_gate.sh, verbatim (each must be defined
#          rule-checks       / exactly once there), so this scope cannot drift from the full gate
#          gate-tool-checks  tools/test_gate_tool_checks.sh
#          full-gate-poplog  tools/test_full_gate_poplog.sh
# Each step is one `poprun --agent A --job J --step NAME -- bash -o pipefail -c CMD` run in REPO; wrap the whole run in
# poprun for one log with a PASS/FAIL footer and gate-commit, e.g.
#   PYTHONDONTWRITEBYTECODE=1 ~/pop-logs/bin/poprun --agent spence --job scoped-gate-rules -- bash tools/scoped_gate.sh rules
# Every step runs even after a failure. Env: REPO (default: this checkout), POPRUN (default ~/pop-logs/bin/poprun),
# POPLOG_AGENT (default spence), GATE_JOB (default scoped-gate-SCOPE).
# rc 0 all steps pass, 1 a step failed, 2 usage error, missing poprun, or a step full_gate.sh does not define once.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
repo="${REPO:-$here/..}"
usage() { echo "usage: tools/scoped_gate.sh SCOPE | --list SCOPE | --dry-run SCOPE   (scopes: rules)" >&2; exit 2; }
mode=run
case "${1:-}" in --list) mode=list; shift;; --dry-run) mode=dry; shift;; esac
[ $# = 1 ] || usage
scope=$1
case "$scope" in
  rules) steps=("full-gate-syntax|bash -n tools/full_gate.sh" "rules-scan|@full_gate" "rule-checks|@full_gate"
                "gate-tool-checks|tools/test_gate_tool_checks.sh" "full-gate-poplog|tools/test_full_gate_poplog.sh");;
  *) echo "scoped_gate: unknown scope '$scope'" >&2; usage;;
esac
cd "$repo" || { echo "scoped_gate: cannot cd to REPO $repo" >&2; exit 2; }
names=(); cmds=()
for s in "${steps[@]}"; do
  n=${s%%|*}; c=${s#*|}
  if [ "$c" = @full_gate ]; then
    k=$(grep -c "^run $n \"" tools/full_gate.sh 2>/dev/null)
    [ "$k" = 1 ] || { echo "scoped_gate: tools/full_gate.sh defines step $n ${k:-0} time(s), want exactly 1" >&2; exit 2; }
    c=$(sed -n "s/^run $n \"\(.*\)\"\$/\1/p" tools/full_gate.sh)
    [ -n "$c" ] || { echo "scoped_gate: cannot read the command of step $n in tools/full_gate.sh" >&2; exit 2; }
  fi
  names+=("$n"); cmds+=("$c")
done
case $mode in
  list) printf '%s\n' "${names[@]}"; exit 0;;
  dry) for i in "${!names[@]}"; do echo "${names[$i]}: ${cmds[$i]}"; done; exit 0;;
esac
poprun="${POPRUN:-$HOME/pop-logs/bin/poprun}"
[ -x "$poprun" ] || { echo "scoped_gate: poprun not found or not executable: $poprun" >&2; exit 2; }
agent="${POPLOG_AGENT:-spence}"; job="${GATE_JOB:-scoped-gate-$scope}"
failed=()
for i in "${!names[@]}"; do
  echo "scoped_gate: step ${names[$i]}: ${cmds[$i]}"
  "$poprun" --agent "$agent" --job "$job" --step "${names[$i]}" -- bash -o pipefail -c "${cmds[$i]}"; rc=$?
  if [ $rc = 0 ]; then echo "scoped_gate: ok   ${names[$i]}"; else echo "scoped_gate: FAIL ${names[$i]} (rc $rc)"; failed+=("${names[$i]}"); fi
done
if [ ${#failed[@]} -gt 0 ]; then echo "scoped_gate: FAIL $scope (${#failed[@]} of ${#names[@]} step(s) failed: ${failed[*]})"; exit 1; fi
echo "scoped_gate: PASS $scope (${#names[@]} step(s))"
exit 0
