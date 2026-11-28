#!/usr/bin/env bash
# Selftest for tools/rule_checks_daily.sh, the daily full drift sweep (Jeff 2026-11-25 5:20 PM CT: "Yes, add a daily
# full sweep over every worktree that fails on any drift"). It runs drift_check --all (every tree, every branch, no
# active-mode skipping) over every ROOT (a dir or a glob of dirs), and FAILs (rc 1) on any drift; rc 2 on a usage or
# git error.
# Fixtures live under $ZF_SCRATCH (default $HOME/.cache/zenith-float), never /tmp.
# RULE_CHECKS_DAILY / DRIFT_CHECK point at other copies (the mutation runner uses this).
set -u
here="$(cd "$(dirname "$0")" && pwd)"
daily="${RULE_CHECKS_DAILY:-$here/rule_checks_daily.sh}"
export DRIFT_CHECK="${DRIFT_CHECK:-$here/drift_check.py}"
root="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/daily-selftest"
mkdir -p "$root" || exit 2
work=$(mktemp -d "$root/run.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
NOW=$(date +%s)
# fixture identity and config, local to this process (never export these in a shell that makes real commits)
export GIT_AUTHOR_NAME=t GIT_AUTHOR_EMAIL=t@t GIT_COMMITTER_NAME=t GIT_COMMITTER_EMAIL=t@t
export GIT_CONFIG_NOSYSTEM=1 GIT_CONFIG_GLOBAL=/dev/null GIT_CEILING_DIRECTORIES="$work"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "ok   $1"; }
bad() { fail=$((fail+1)); echo "FAIL $1"; }
show() { printf '%s\n' "$1" | sed 's/^/     | /'; }
commit_at() {
  local r=$1 ts=$((NOW-$2)); echo "$3" >> "$r/f.txt"; git -C "$r" add f.txt
  GIT_COMMITTER_DATE="@$ts +0000" GIT_AUTHOR_DATE="@$ts +0000" git -C "$r" commit -q -m "$3"
}
# tree DIR [unpushed] : a work tree at DIR with one pushed commit (origin bare repo beside it under $work/origins),
# idle 2 days (commit, index, logs/HEAD); with "unpushed", 2 more commits 2 days old that are on no remote
tree() {
  local d=$1 o; o="$work/origins/$(echo "$1" | tr / _).git"
  mkdir -p "$(dirname "$d")"; git init -q --bare "$o"; git init -q -b w "$d"; git -C "$d" remote add origin "$o"
  commit_at "$d" 172800 base; git -C "$d" push -q -u origin w 2>/dev/null
  [ "${2:-}" = unpushed ] && { commit_at "$d" 172800 u1; commit_at "$d" 172800 u2; }
  touch -d "@$((NOW-172800))" "$d/.git/index" "$d/.git/logs/HEAD"
}
# sweep ARGS... : run the daily sweep with a fixed --now
# every sweep is capped at 60 s (normal sweeps take a few s), so a sweep that hangs FAILs that case (rc 124)
# instead of hanging the selftest
sweep() { timeout 60 bash "$daily" "$@" -- --now "$NOW"; }
# check NAME WANT_RC [re:REGEX]... [nore:REGEX]... -- ARGS...
check() {
  local name=$1 want=$2; shift 2; local res=() nres=()
  while [ "$1" != -- ]; do case "$1" in re:*) res+=("${1#re:}");; nore:*) nres+=("${1#nore:}");; esac; shift; done; shift
  local out rc good=1; out=$(sweep "$@" 2>&1); rc=$?
  [ "$rc" = "$want" ] || good=0
  local x; for x in "${res[@]+"${res[@]}"}"; do grep -Eq -- "$x" <<<"$out" || good=0; done
  for x in "${nres[@]+"${nres[@]}"}"; do grep -Eq -- "$x" <<<"$out" && good=0; done
  if [ $good = 1 ]; then ok "$name (rc=$rc)"; else bad "$name (want rc=$want, got rc=$rc)"; show "$out"; fi
}
# ---- the sweep is --all: an idle tree with drift FAILs (the 15-minute active mode skips it) ----
tree "$work/r1/idle" unpushed; tree "$work/r1/clean"
check "idle tree with 2 unpushed commits 2 days old: FAIL" 1 re:"FAIL $work/r1/idle " re:'^rule_checks_daily: FAIL' nore:"FAIL $work/r1/clean" -- "$work/r1"
python3 -B "$DRIFT_CHECK" --now "$NOW" --scan "$work/r1" >/dev/null 2>&1; rc=$?
[ $rc = 0 ] && ok "control: active mode skips the same idle tree (rc=0)" || bad "control: active mode (rc=$rc)"
tree "$work/r2/obr"; git -C "$work/r2/obr" switch -q -c old; commit_at "$work/r2/obr" 172800 o1; commit_at "$work/r2/obr" 172800 o2
git -C "$work/r2/obr" switch -q w; touch -d "@$((NOW-172800))" "$work/r2/obr/.git/index" "$work/r2/obr/.git/logs/HEAD"
check "idle branch (not HEAD) 2 ahead in an idle tree: FAIL" 1 re:"FAIL $work/r2/obr old" -- "$work/r2"
check "clean roots: PASS" 0 re:'^rule_checks_daily: PASS' -- "$work/r1/clean"
# ---- roots: dirs and globs of dirs; every root is scanned; any FAIL in any root fails the sweep ----
tree "$work/hm/a-work/t1" unpushed; tree "$work/hm/b-work/t2"; echo notadir > "$work/hm/afile"
check "glob root: each matching dir is scanned, a file match is ignored" 1 re:"FAIL $work/hm/a-work/t1 " re:"ok   $work/hm/b-work/t2" -- "$work/hm/*"
check "bad tree in the 2nd root is not hidden by a clean 1st root" 1 re:"FAIL $work/r1/idle " -- "$work/r1/clean" "$work/r1"
out=$(sweep "$work/r1" "$work/hm/*" 2>&1); rc=$?
k=$(grep -c '^drift_check: FAIL ' <<<"$out")
[ $rc = 1 ] && [ "$k" = 2 ] && grep -q "2 FAIL line(s)" <<<"$out" && ok "summary counts the FAIL lines over all roots (2)" || { bad "summary FAIL count (rc=$rc, $k FAIL lines)"; show "$out"; }
# a mistyped glob must not silently drop its root (Cassie round-3 N4; this case was rc 0 "ignored" until round 3)
check "a glob matching nothing beside a real root fails loudly: rc 2" 2 re:"no-such-\\* matches no directory" -- "$work/no-such-*" "$work/r1/clean"
check "only a glob matching nothing: no roots, rc 2" 2 -- "$work/no-such-*"
check "a literal root that is missing: rc 2" 2 -- "$work/no-such-root" "$work/r1/clean"
check "a literal root that is a file: rc 2" 2 -- "$work/hm/afile"
check "no roots at all: rc 2, says so (Z17)" 2 re:'^rule_checks_daily: no roots to sweep' --
# ---- errors fail closed ----
mkdir -p "$work/e/broken"; echo "gitdir: $work/no-such-gitdir" > "$work/e/broken/.git"
check "a drift_check error (unreadable child): rc 2" 2 re:"$work/e/broken" -- "$work/e"
check "an unknown option: rc 2" 2 -- --bogus "$work/r1"
DRIFT_CHECK="$work/no-such-drift_check.py" sweep "$work/r1/clean" >/dev/null 2>&1; rc=$?
[ $rc = 2 ] && ok "drift_check missing: rc 2" || bad "drift_check missing (want rc=2, got rc=$rc)"
out=$(DRIFT_CHECK="$work/no-such-drift_check.py" bash "$daily" "$work/r1/clean" 2>&1); rc=$?
[ $rc = 2 ] && grep -q 'drift_check not found' <<<"$out" && ok "drift_check missing: rc 2, says so (W12)" || { bad "drift_check missing message (rc=$rc)"; show "$out"; }
# ---- round 3 (Cassie N4): every glob must match at least one directory ----
mkdir -p "$work/fo"; echo f > "$work/fo/onlyfile"
check "a glob that matches only files beside a real root: rc 2" 2 re:'matches no directory' -- "$work/fo/*" "$work/r1/clean"
check "a typo glob after a drifted root still fails rc 2, not 1" 2 re:'matches no directory' -- "$work/r1" "$work/r1-typo*"
err=$(bash "$daily" "$work/no-such-*" "$work/r1/clean" -- --now "$NOW" 2>&1 >/dev/null); rc=$?
[ $rc = 2 ] && grep -q 'matches no directory' <<<"$err" && ok "the 'matches no directory' error goes to stderr (X17)" || { bad "'matches no directory' on stderr (rc=$rc)"; show "$err"; }
check "a bracket glob root is expanded (X19)" 1 re:"FAIL $work/hm/a-work/t1 " re:"ok   $work/hm/b-work/t2" -- "$work/hm/[ab]-work"
# ---- round 4 N2 (Cassie): overlapping roots (the same dir by a literal, a glob, ./ or a symlink) are swept once ----
ln -s "$work/hm/a-work" "$work/alias"
out=$(sweep "$work/hm/a-work" "$work/hm/*" "$work/hm/./b-work" "$work/alias" 2>&1); rc=$?
k=$(grep -c "^drift_check: FAIL $work/hm/a-work/t1 " <<<"$out")
[ $rc = 1 ] && [ "$k" = 1 ] && grep -q "sweeping 2 root(s) with drift_check --all: $work/hm/a-work $work/hm/b-work\$" <<<"$out" && grep -q "1 FAIL line(s), 0 ERROR line(s), 2 root(s);" <<<"$out" && ok "N2: overlapping roots are swept once (2 roots, 1 FAIL)" || { bad "N2: overlapping roots (rc=$rc, $k FAIL lines for t1)"; show "$out"; }
# the sweep itself refuses a literal root that is not a dir (drift_check would too, but later and less plainly)
check "a missing literal root is refused by the sweep itself" 2 re:"^rule_checks_daily: root $work/no-such-root is not a directory" -- "$work/no-such-root" "$work/r1/clean"
# ---- allowlist (Jeff 2026-11-27 4:42 PM CT): run-daily.sh passes --only-origin jscarr64/zenith-float to drift_check,
# and its root is the glob of agent folders ----
tree "$work/al/a-work/acc" unpushed; git -C "$work/al/a-work/acc" remote set-url origin https://github.com/jscarr64/accumath.git
tree "$work/al/b-work/zf" unpushed; git -C "$work/al/b-work/zf" remote set-url origin https://github.com/jscarr64/zenith-float.git
out=$(timeout 60 bash "$daily" "$work/al/*-work" -- --now "$NOW" --only-origin jscarr64/zenith-float 2>&1); rc=$?
[ $rc = 1 ] && grep -q "FAIL $work/al/b-work/zf " <<<"$out" && grep -q "skip $work/al/a-work/acc (" <<<"$out" && ! grep -q "FAIL $work/al/a-work" <<<"$out" && ok "allowlist: the sweep FAILs the idle zenith-float clone and skips the other origin (rc=1)" || { bad "allowlist: sweep with --only-origin (want rc=1, zf FAIL, acc skip; got rc=$rc)"; show "$out"; }
mkdir -p "$work/al0"
check "allowlist: zero agent folders: the root glob matches no directory, rc 2" 2 re:'matches no directory' -- "$work/al0/*-work"
# ---- round 7 N1 (Cassie): a root that is a symlink loop is "not a directory", rc 2 ----
mkdir -p "$work/lp"; ln -s l "$work/lp/l"
check "r7 N1: a root that is a symlink loop: rc 2" 2 re:"^rule_checks_daily: root $work/lp/l is not a directory" -- "$work/lp/l" "$work/r1/clean"
echo "rule_checks_daily_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
