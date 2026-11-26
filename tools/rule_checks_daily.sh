#!/usr/bin/env bash
# rule_checks_daily.sh: the daily full drift sweep (Jeff 2026-11-25 5:20 PM CT: "Yes, add a daily full sweep over
# every worktree that fails on any drift"). The 15-minute check runs drift_check in active mode, which skips idle
# trees and branches; this sweep runs drift_check --all, so long-abandoned unpushed or dirty work FAILs once a day.
# usage: tools/rule_checks_daily.sh [--exclude NAME]... ROOT_OR_GLOB... [-- DRIFT_CHECK_ARGS...]
#   ROOT_OR_GLOB  a directory, or a glob of directories (quoted, expanded here; non-dir matches are ignored, a glob
#                 that matches no directory is an error, rc 2). Each root goes to drift_check as --scan ROOT (ROOT and its children).
#   --exclude N   leave out any root that has a path component N, as given or after resolving symlinks. MathWorks is
#                 always excluded (Jeff: MathWorks is not touched by this cleanup); --exclude adds names. The names
#                 also go to drift_check --exclude, so a scanned child with such a component is left out too.
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
excluded() {  # PATH: 0 when PATH, as written or followed symlink by symlink, has an excluded component.
  # Decided from names and readlink strings only (Cassie r4): each component is checked BEFORE it is touched
  # ([ -L ] / readlink, never -f), so nothing excluded, and nothing under it, is ever stat'ed or listed.
  local rest=$1 cur= c t x hops=0
  case "$rest" in /*) ;; *) rest="$PWD/$rest";; esac
  while :; do
    while [ "${rest#/}" != "$rest" ]; do rest=${rest#/}; done
    [ -n "$rest" ] || return 1
    c=${rest%%/*}; rest=${rest#"$c"}
    case "$c" in .) continue;; ..) cur=${cur%/*}; continue;; esac
    for x in "${excl[@]}"; do [ "$c" = "$x" ] && return 0; done
    if [ -L "$cur/$c" ]; then
      hops=$((hops+1)); [ $hops -le 40 ] || return 1   # a symlink loop; it fails later as not a directory
      t=$(readlink -- "$cur/$c")
      case "$t" in /*) cur=;; esac
      rest="$t/$rest"
    else
      cur="$cur/$c"
    fi
  done
}
roots=(); nex=0
for p in "${pats[@]+"${pats[@]}"}"; do
  # the pattern's own name decides first, before it is expanded or stat'ed (Cassie r4)
  if excluded "$p"; then echo "rule_checks_daily: excluded $p"; nex=$((nex+1)); continue; fi
  case "$p" in
    *[\*\?\[]*) mapfile -t ms < <(compgen -G "$p" | sort)
                nd=0
                for m in "${ms[@]+"${ms[@]}"}"; do
                  # each match is judged by name before [ -d ] can follow it into MathWorks
                  if excluded "$m"; then echo "rule_checks_daily: excluded $m"; nex=$((nex+1)); nd=$((nd+1)); continue; fi
                  [ -d "$m" ] && { nd=$((nd+1)); roots+=("$m"); }
                done
                # a mistyped glob must not silently drop a root (Cassie round-3 N4)
                [ $nd -gt 0 ] || { echo "rule_checks_daily: root glob $p matches no directory" >&2; exit 2; };;
    *) [ -d "$p" ] || { echo "rule_checks_daily: root $p is not a directory" >&2; exit 2; }; roots+=("$p");;
  esac
done
[ ${#roots[@]} -gt 0 ] || { echo "rule_checks_daily: no roots to sweep" >&2; exit 2; }
args=(); for r in "${roots[@]}"; do args+=(--scan "$r"); done
for x in "${excl[@]}"; do args+=(--exclude "$x"); done   # scanned children too, not only roots
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
