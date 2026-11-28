#!/usr/bin/env bash
# Selftest for tools/scoped_gate.sh (Cassie round-2 note (c): the scoped gate must be reproducible from the PR head).
# Scope "rules" = full-gate-syntax (bash -n tools/full_gate.sh), rules-scan and rule-checks (each taken verbatim from
# its `run NAME "CMD"` line in tools/full_gate.sh, so the scoped gate cannot drift from the full gate), gate-tool-checks
# (tools/test_gate_tool_checks.sh) and full-gate-poplog (tools/test_full_gate_poplog.sh). Each step is one poprun
# --step (POPRUN), run in REPO; every step runs even after a failure; rc 1 when any step fails; rc 2 on a usage error,
# a step that full_gate.sh does not define exactly once, or a missing poprun.
# poprun is stubbed here; fixtures live under $ZF_SCRATCH (default $HOME/.cache/zenith-float), never /tmp.
# SCOPED_GATE points at another copy (the mutation runner uses this).
set -u
# the GATE_JOB / POPLOG_AGENT defaults are under test: a caller's values (an outer poprun, a gate run with
# GATE_JOB set) must not leak in
unset GATE_JOB POPLOG_AGENT
here="$(cd "$(dirname "$0")" && pwd)"
sg="${SCOPED_GATE:-$here/scoped_gate.sh}"
root="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/scoped-gate-selftest"
mkdir -p "$root" || exit 2
work=$(mktemp -d "$root/run.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "ok   $1"; }
bad() { fail=$((fail+1)); echo "FAIL $1"; }
show() { printf '%s\n' "$1" | sed 's/^/     | /'; }
# stub poprun: records its options, then runs the command after --
cat > "$work/poprun" <<'STUB'
#!/usr/bin/env bash
a=(); while [ "$1" != -- ]; do a+=("$1"); shift; done; shift
echo "${a[*]}" >> "$POPRUN_LOG"
"$@"
STUB
chmod +x "$work/poprun"
# fake REPO with a full_gate.sh; RC_CMD / RS_CMD override its two step commands
mkrepo() {
  local d="$work/$1"; mkdir -p "$d/tools"
  { echo '#!/usr/bin/env bash'; echo 'run() { :; }'
    echo "run rules-scan \"${RS_CMD:-echo RS-ran}\""; echo "run rule-checks \"${RC_CMD:-echo RC-ran}\""; } > "$d/tools/full_gate.sh"
  printf '#!/usr/bin/env bash\necho "GTC-ran in $PWD"\n' > "$d/tools/test_gate_tool_checks.sh"
  printf '#!/usr/bin/env bash\necho FGP-ran\n' > "$d/tools/test_full_gate_poplog.sh"
  chmod +x "$d/tools/"*.sh; echo "$d"
}
# go REPO ARGS... : run the scoped gate against REPO with the stub; sets out, rc, calls
go() { local r=$1; shift; : > "$work/calls"; out=$(cd "$work" && REPO="$r" POPRUN="$work/poprun" POPRUN_LOG="$work/calls" bash "$sg" "$@" 2>&1); rc=$?; calls=$(cat "$work/calls"); }
want_steps="full-gate-syntax rules-scan rule-checks gate-tool-checks full-gate-poplog"
# ---- --list ----
r=$(mkrepo r0); go "$r" --list rules
[ $rc = 0 ] && [ "$(echo $out)" = "$want_steps" ] && ok "--list rules: the 5 steps in order" || { bad "--list rules (rc=$rc)"; show "$out"; }
# ---- a passing run: every step is one poprun --step, in order, in REPO, with the full_gate.sh commands ----
go "$r" rules
steps=$(sed -n 's/.*--step \([^ ]*\).*/\1/p' <<<"$calls" | tr '\n' ' ' | sed 's/ $//')
[ $rc = 0 ] && [ "$steps" = "$want_steps" ] && ok "rules: 5 poprun steps in order, rc 0" || { bad "rules: steps (rc=$rc, steps '$steps')"; show "$out"; show "$calls"; }
grep -q RS-ran <<<"$out" && grep -q RC-ran <<<"$out" && ok "rules-scan and rule-checks run the commands from full_gate.sh" || { bad "full_gate.sh commands not run"; show "$out"; }
grep -q "GTC-ran in $r$" <<<"$out" && grep -q FGP-ran <<<"$out" && ok "gate-tool-checks and full-gate-poplog run, in REPO" || { bad "meta steps not run in REPO"; show "$out"; }
k=$(grep -c -- '--agent spence --job scoped-gate-rules' <<<"$calls")
[ "$k" = 5 ] && ok "default --agent spence --job scoped-gate-rules on every step" || { bad "agent/job defaults ($k of 5)"; show "$calls"; }
: > "$work/calls"; out=$(cd "$work" && REPO="$r" POPRUN="$work/poprun" POPRUN_LOG="$work/calls" POPLOG_AGENT=cassie GATE_JOB=jx bash "$sg" rules 2>&1)
k=$(grep -c -- '--agent cassie --job jx' "$work/calls"); [ "$k" = 5 ] && ok "POPLOG_AGENT and GATE_JOB are passed on" || { bad "POPLOG_AGENT/GATE_JOB ($k of 5)"; show "$(cat "$work/calls")"; }
# ---- failures: every step still runs, rc 1, the failed step is named ----
r=$(RC_CMD="exit 3" mkrepo r1); go "$r" rules
[ $rc = 1 ] && grep -q FGP-ran <<<"$out" && grep -q 'FAIL.*rule-checks' <<<"$out" && ok "a failing rule-checks: rc 1, named, later steps still run" || { bad "failing rule-checks (rc=$rc)"; show "$out"; }
r=$(RS_CMD="exit 1" mkrepo r2); go "$r" rules
[ $rc = 1 ] && grep -q RC-ran <<<"$out" && ok "a failing rules-scan: rc 1, rule-checks still runs" || { bad "failing rules-scan (rc=$rc)"; show "$out"; }
r=$(mkrepo r3); echo 'if then fi' >> "$r/tools/full_gate.sh"; go "$r" rules
[ $rc = 1 ] && grep -q 'FAIL.*full-gate-syntax' <<<"$out" && ok "full_gate.sh syntax error: full-gate-syntax fails, rc 1" || { bad "full_gate.sh syntax error (rc=$rc)"; show "$out"; }
r=$(mkrepo r4); printf '#!/usr/bin/env bash\nexit 1\n' > "$r/tools/test_full_gate_poplog.sh"; go "$r" rules
[ $rc = 1 ] && grep -q 'FAIL.*full-gate-poplog' <<<"$out" && ok "a failing last step: rc 1" || { bad "failing last step (rc=$rc)"; show "$out"; }
# ---- full_gate.sh must define each borrowed step exactly once, checked before anything runs ----
r=$(mkrepo r5); sed -i '/^run rule-checks /d' "$r/tools/full_gate.sh"; go "$r" rules
[ $rc = 2 ] && [ -z "$calls" ] && grep -q rule-checks <<<"$out" && ok "step missing from full_gate.sh: rc 2, nothing run" || { bad "missing step (rc=$rc)"; show "$out"; show "$calls"; }
r=$(mkrepo r6); echo 'run rules-scan "echo again"' >> "$r/tools/full_gate.sh"; go "$r" rules
[ $rc = 2 ] && [ -z "$calls" ] && ok "step defined twice in full_gate.sh: rc 2, nothing run" || { bad "duplicate step (rc=$rc)"; show "$out"; }
# ---- usage ----
r=$(mkrepo r7)
go "$r" bogus;          [ $rc = 2 ] && [ -z "$calls" ] && ok "unknown scope: rc 2" || bad "unknown scope (rc=$rc)"
go "$r";                [ $rc = 2 ] && ok "no scope: rc 2" || bad "no scope (rc=$rc)"
go "$r" --list bogus;   [ $rc = 2 ] && ok "--list of an unknown scope: rc 2" || bad "--list unknown scope (rc=$rc)"
out=$(REPO="$r" POPRUN="$work/no-such-poprun" bash "$sg" rules 2>&1); rc=$?
[ $rc = 2 ] && ok "poprun missing: rc 2" || { bad "poprun missing (rc=$rc)"; show "$out"; }
# ---- the real repo: --dry-run resolves both borrowed steps from the real tools/full_gate.sh ----
out=$(REPO="$here/.." bash "$sg" --dry-run rules 2>&1); rc=$?
[ $rc = 0 ] && grep -q '^rule-checks: .*tools/rule_checks_mutation.sh' <<<"$out" && grep -q '^rules-scan: .*tools/rules_scan.sh' <<<"$out" \
  && ok "--dry-run on this repo: rules-scan and rule-checks resolve from tools/full_gate.sh" || { bad "--dry-run on this repo (rc=$rc)"; show "$out"; }
r=$(mkrepo r8); sed -i 's/^run rule-checks .*/run rule-checks ""/' "$r/tools/full_gate.sh"; go "$r" rules
[ $rc = 2 ] && [ -z "$calls" ] && ok "a borrowed step with an empty command: rc 2, nothing run" || { bad "empty borrowed command (rc=$rc)"; show "$out"; }
# ---- round 3 (Cassie W22): a borrowed step runs under pipefail, so a failure inside a pipeline fails it ----
r=$(RC_CMD="false | true" mkrepo r9); go "$r" rules
[ $rc = 1 ] && grep -q 'FAIL.*rule-checks' <<<"$out" && ok "a borrowed step 'false | true' fails under pipefail: rc 1 (W22)" || { bad "pipefail (want rc=1, got rc=$rc)"; show "$out"; }
echo "scoped_gate_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
