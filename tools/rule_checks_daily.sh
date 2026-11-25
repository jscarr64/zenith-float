#!/usr/bin/env bash
# rule_checks_daily.sh: the daily full drift sweep (Jeff 2026-11-25 5:20 PM CT: "Yes, add a daily full sweep over
# every worktree that fails on any drift"). The 15-minute check runs drift_check in active mode, which skips idle
# trees and branches; this sweep runs drift_check --all, so long-abandoned unpushed or dirty work FAILs once a day.
# usage: tools/rule_checks_daily.sh [--exclude NAME]... ROOT_OR_GLOB... [-- DRIFT_CHECK_ARGS...]
#   ROOT_OR_GLOB  a directory, or a glob of directories (quoted, expanded here; non-dir matches are ignored, a glob
#                 that matches nothing is noted). Each root goes to drift_check as --scan ROOT (ROOT and its children).
#   --exclude N   leave out any root that has a path component N, as given or after resolving symlinks. MathWorks is
#                 always excluded (Jeff: MathWorks is not touched by this cleanup); --exclude adds names.
#   -- ARGS       passed on to drift_check (the selftest passes --now).
# DRIFT_CHECK selects drift_check.py (default: beside this script). Read-only, like drift_check.
# rc 0 no drift, 1 drift (any FAIL in any root), 2 usage error, no roots, or a drift_check / git error (fail closed).
set -u
here="$(cd "$(dirname "$0")" && pwd)"
dc="${DRIFT_CHECK:-$here/drift_check.py}"
excl=(MathWorks); pats=(); extra=()
while [ $# -gt 0 ]; do
  case "$1" in
    --exclude) [ $# -ge 2 ] || { echo "rule_checks_daily: --exclude needs a NAME" >&2; exit 2; }; excl+=("$2"); shift 2;;
    --) shift; extra=("$@"); break;;
    -*) echo "rule_checks_daily: unknown option $1" >&2; exit 2;;
    *) pats+=("$1"); shift;;
  esac
done
[ -f "$dc" ] || { echo "rule_checks_daily: drift_check not found: $dc" >&2; exit 2; }
excluded() {  # PATH: 0 when PATH, as given or resolved, has an excluded component
  local p x real; real=$(readlink -f -- "$1" 2>/dev/null || echo "$1")
  for p in "$1" "$real"; do for x in "${excl[@]}"; do case "/$p/" in */"$x"/*) return 0;; esac; done; done
  return 1
}
roots=(); nex=0
for p in "${pats[@]+"${pats[@]}"}"; do
  case "$p" in
    *[\*\?\[]*) mapfile -t ms < <(compgen -G "$p" | sort)
                [ ${#ms[@]} -gt 0 ] || echo "rule_checks_daily: note: $p matches nothing";;
    *) [ -d "$p" ] || { echo "rule_checks_daily: root $p is not a directory" >&2; exit 2; }; ms=("$p");;
  esac
  for m in "${ms[@]+"${ms[@]}"}"; do
    [ -d "$m" ] || continue
    if excluded "$m"; then echo "rule_checks_daily: excluded $m"; nex=$((nex+1)); continue; fi
    roots+=("$m")
  done
done
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
echo "rule_checks_daily: $v (drift_check rc $rc, $nf FAIL line(s), $ne ERROR line(s), ${#roots[@]} root(s), $nex excluded; drift_check: ${dsum:-no summary})"
exit $rc
