#!/usr/bin/env bash
# rule_checks_daily.sh: the daily full drift sweep (Jeff 2026-11-25 5:20 PM CT: "Yes, add a daily full sweep over
# every worktree that fails on any drift"). The 15-minute check runs drift_check in active mode, which skips idle
# trees and branches; this sweep runs drift_check --all, so long-abandoned unpushed or dirty work FAILs once a day.
# usage: tools/rule_checks_daily.sh ROOT_OR_GLOB... [-- DRIFT_CHECK_ARGS...]
#   ROOT_OR_GLOB  a directory, or a glob of directories (quoted, expanded here; non-dir matches are ignored, a glob
#                 that matches no directory is an error, rc 2). Each root goes to drift_check as --scan ROOT (ROOT and its children).
#   -- ARGS       passed on to drift_check (run-daily.sh passes --only-origin, the selftest passes --now).
# DRIFT_CHECK selects drift_check.py (default: beside this script). Read-only, like drift_check.
# rc 0 no drift, 1 drift (any FAIL in any root), 2 usage error, no roots, or a drift_check / git error (fail closed).
set -u
here="$(cd "$(dirname "$0")" && pwd)"
dc="${DRIFT_CHECK:-$here/drift_check.py}"
pats=(); extra=()
while [ $# -gt 0 ]; do
  case "$1" in
    --) shift; extra=("$@"); break;;
    -*) echo "rule_checks_daily: unknown option $1" >&2; exit 2;;
    *) pats+=("$1"); shift;;
  esac
done
[ -f "$dc" ] || { echo "rule_checks_daily: drift_check not found: $dc" >&2; exit 2; }
roots=()
for p in "${pats[@]+"${pats[@]}"}"; do
  case "$p" in
    *[\*\?\[]*) mapfile -t ms < <(compgen -G "$p" | sort)
                nd=0
                for m in "${ms[@]+"${ms[@]}"}"; do [ -d "$m" ] && { nd=$((nd+1)); roots+=("$m"); }; done
                # a mistyped glob must not silently drop a root (Cassie round-3 N4)
                [ $nd -gt 0 ] || { echo "rule_checks_daily: root glob $p matches no directory" >&2; exit 2; };;
    *) [ -d "$p" ] || { echo "rule_checks_daily: root $p is not a directory" >&2; exit 2; }; roots+=("$p");;
  esac
done
# the same dir named twice (a literal and a glob, ./, a symlink) is swept once, first spelling kept (Cassie r4 N2).
declare -A seen=(); uniq=()
for r in "${roots[@]+"${roots[@]}"}"; do
  k=$(readlink -f -- "$r"); [ -n "${seen[$k]+x}" ] && continue
  seen[$k]=1; uniq+=("$r")
done
roots=("${uniq[@]+"${uniq[@]}"}")
[ ${#roots[@]} -gt 0 ] || { echo "rule_checks_daily: no roots to sweep" >&2; exit 2; }
args=(); for r in "${roots[@]}"; do args+=(--scan "$r"); done
echo "rule_checks_daily: sweeping ${#roots[@]} root(s) with drift_check --all: ${roots[*]}"
out=$(python3 -B "$dc" --all "${args[@]}" "${extra[@]+"${extra[@]}"}" 2>&1); rc=$?
# drift_check's own one-line summary ("drift_check: ok|FAIL (N repo(s) checked, ...)") is folded into the sweep's
# summary, so every "drift_check: FAIL" line left in the log is one finding
sumre='^drift_check: (ok|FAIL) \('
grep -Ev "$sumre" <<<"$out"
dsum=$(grep -E "$sumre" <<<"$out" | tail -1 | sed 's/^drift_check: //')
nf=$(grep -Ev "$sumre" <<<"$out" | grep -c '^drift_check: FAIL '); ne=$(grep -c '^drift_check: ERROR ' <<<"$out")
case $rc in 0) v=PASS;; 1) v=FAIL;; *) v=ERROR; rc=2;; esac
echo "rule_checks_daily: $v (drift_check rc $rc, $nf FAIL line(s), $ne ERROR line(s), ${#roots[@]} root(s); drift_check: ${dsum:-no summary})"
exit $rc
