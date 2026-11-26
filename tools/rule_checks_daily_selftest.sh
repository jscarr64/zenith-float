#!/usr/bin/env bash
# Selftest for tools/rule_checks_daily.sh, the daily full drift sweep (Jeff 2026-11-25 5:20 PM CT: "Yes, add a daily
# full sweep over every worktree that fails on any drift"). It runs drift_check --all (every tree, every branch, no
# active-mode skipping) over every ROOT (a dir or a glob of dirs), leaves out any root with a MathWorks path component
# (or another --exclude NAME), and FAILs (rc 1) on any drift; rc 2 on a usage or git error.
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
sweep() { bash "$daily" "$@" -- --now "$NOW"; }
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
check "no roots at all: rc 2" 2 --
# ---- MathWorks is excluded (as a root, by any path component, also through a symlink); --exclude adds names ----
tree "$work/x/MathWorks/mw" unpushed; tree "$work/x/ok"
check "glob root named MathWorks is excluded" 0 re:"excluded $work/x/MathWorks" nore:"FAIL $work/x/MathWorks" -- "$work/x/*"
check "literal MathWorks root only: excluded, no roots, rc 2" 2 re:"excluded $work/x/MathWorks" -- "$work/x/MathWorks"
mkdir -p "$work/lnk"; ln -s "$work/x/MathWorks" "$work/lnk/mwlink"
check "a root that is a symlink into MathWorks is excluded" 0 re:"excluded $work/lnk/mwlink" -- "$work/lnk/*" "$work/x/ok"
check "MathWorks deeper in the path is excluded" 0 re:"excluded $work/x/MathWorks/mw" -- "$work/x/MathWorks/mw" "$work/x/ok"
tree "$work/y/MathWorksOld/t" unpushed
check "a name that only contains MathWorks is not excluded" 1 re:"FAIL $work/y/MathWorksOld/t " -- "$work/y/MathWorksOld"
tree "$work/z/skipme/t" unpushed; tree "$work/z/keep"
check "--exclude NAME also excludes" 0 re:"excluded $work/z/skipme" -- --exclude skipme "$work/z/*"
check "--exclude keeps the MathWorks default" 0 re:"excluded $work/x/MathWorks" -- --exclude skipme "$work/x/*"
# ---- errors fail closed ----
mkdir -p "$work/e/broken"; echo "gitdir: $work/no-such-gitdir" > "$work/e/broken/.git"
check "a drift_check error (unreadable child): rc 2" 2 re:"$work/e/broken" -- "$work/e"
check "an unknown option: rc 2" 2 -- --bogus "$work/r1"
DRIFT_CHECK="$work/no-such-drift_check.py" sweep "$work/r1/clean" >/dev/null 2>&1; rc=$?
[ $rc = 2 ] && ok "drift_check missing: rc 2" || bad "drift_check missing (want rc=2, got rc=$rc)"
check "every root excluded: says no roots (rc 2)" 2 re:'no roots to sweep' -- "$work/x/MathWorks"
# ---- round 3 (Cassie W10, W12) ----
out=$(bash "$daily" --exclude 2>&1); rc=$?
[ $rc = 2 ] && grep -q -- '--exclude needs a NAME' <<<"$out" && ok "--exclude with no NAME: rc 2, says so (W10)" || { bad "--exclude with no NAME (want rc=2, got rc=$rc)"; show "$out"; }
out=$(DRIFT_CHECK="$work/no-such-drift_check.py" bash "$daily" "$work/r1/clean" 2>&1); rc=$?
[ $rc = 2 ] && grep -q 'drift_check not found' <<<"$out" && ok "drift_check missing: rc 2, says so (W12)" || { bad "drift_check missing message (rc=$rc)"; show "$out"; }
# ---- round 3 (Jeff's MathWorks rule, Cassie N3): MathWorks is excluded as a scanned CHILD too, not only as a root ----
tree "$work/c/MathWorks" unpushed; tree "$work/c/ok"
check "a child folder named MathWorks under a root is excluded" 0 re:"excluded $work/c/MathWorks" nore:"FAIL $work/c/MathWorks" re:"ok   $work/c/ok" -- "$work/c"
tree "$work/c2/skipme" unpushed; tree "$work/c2/keep"
check "--exclude NAME also excludes a child folder" 0 re:"excluded $work/c2/skipme" -- --exclude skipme "$work/c2"
# ---- round 3 (Cassie N4): every glob must match at least one directory ----
mkdir -p "$work/fo"; echo f > "$work/fo/onlyfile"
check "a glob that matches only files beside a real root: rc 2" 2 re:'matches no directory' -- "$work/fo/*" "$work/r1/clean"
check "a typo glob after a drifted root still fails rc 2, not 1" 2 re:'matches no directory' -- "$work/r1" "$work/r1-typo*"
# the sweep's own root exclusion resolves symlinks (drift_check --exclude would also hide it, so check the sweep's own line)
check "a lone root symlinked into MathWorks: the sweep excludes it, no roots, rc 2" 2 re:"^rule_checks_daily: excluded $work/lnk/mwlink" -- "$work/lnk/mwlink"
# ---- round 4 (Cassie X16, X17, X19) ----
mkdir -p "$work/h2"; ln -s "$work/x/MathWorks/mw" "$work/h2/l2"; ln -s "$work/h2/l2" "$work/h2/l1"
check "a lone root reached by a 2-hop symlink chain into MathWorks: excluded by the sweep, rc 2 (X16)" 2 re:"^rule_checks_daily: excluded $work/h2/l1" -- "$work/h2/l1"
ln -s "$work/x/MathWorks" "$work/sp"
check "a lone root under a symlinked parent that points into MathWorks: excluded by the sweep, rc 2 (X16)" 2 re:"^rule_checks_daily: excluded $work/sp/mw" -- "$work/sp/mw"
err=$(bash "$daily" "$work/no-such-*" "$work/r1/clean" -- --now "$NOW" 2>&1 >/dev/null); rc=$?
[ $rc = 2 ] && grep -q 'matches no directory' <<<"$err" && ok "the 'matches no directory' error goes to stderr (X17)" || { bad "'matches no directory' on stderr (rc=$rc)"; show "$err"; }
check "a bracket glob root is expanded (X19)" 1 re:"FAIL $work/hm/a-work/t1 " re:"ok   $work/hm/b-work/t2" -- "$work/hm/[ab]-work"
# ---- round 4 (Cassie): roots are judged by readlink strings hop by hop, never by stat'ing a target in MathWorks ----
tree "$work/far/t" unpushed; mkdir -p "$work/h3"; ln -s "$work/far/t" "$work/x/MathWorks/hop2"; ln -s "$work/x/MathWorks/hop2" "$work/h3/l1"
check "a root whose first hop names MathWorks is excluded, whatever the target (no roots, rc 2)" 2 re:"^rule_checks_daily: excluded $work/h3/l1" nore:"FAIL $work" -- "$work/h3/l1"
check "a glob match whose first hop names MathWorks is excluded, whatever the target" 0 re:"^rule_checks_daily: excluded $work/h3/l1" nore:"FAIL $work" -- "$work/h3/*" "$work/x/ok"
if command -v strace >/dev/null; then
  strace -f -qq -e trace=%file -o "$work/trace" bash "$daily" "$work/h3/l1" "$work/h3/*" "$work/x/MathWorks/mw" "$work/sp/mw" "$work/x/ok" -- --now "$NOW" >/dev/null 2>&1; rc=$?
  hits=$(grep -v 'execve(' "$work/trace" | grep -F "$work/x/MathWorks")
  [ $rc = 0 ] && [ -z "$hits" ] && ok "strace shows no syscall inside MathWorks while choosing roots" || { bad "syscalls inside MathWorks while choosing roots (rc=$rc)"; show "$hits"; }
else
  echo "note strace check skipped (no strace)"
fi
# ---- round 4 N1 (Cassie): a root that is a symlink loop stays rc 2 but says so plainly ----
mkdir -p "$work/lp"; ln -s l "$work/lp/l"
check "a root that is a symlink loop: rc 2, 'is a symlink loop'" 2 re:"root $work/lp/l is a symlink loop" -- "$work/lp/l"
check "a glob match that is a symlink loop: rc 2, 'is a symlink loop'" 2 re:"root $work/lp/l is a symlink loop" -- "$work/lp/*"
# ---- round 4 N2 (Cassie): overlapping roots (the same dir by a literal, a glob, ./ or a symlink) are swept once ----
ln -s "$work/hm/a-work" "$work/alias"
out=$(sweep "$work/hm/a-work" "$work/hm/*" "$work/hm/./b-work" "$work/alias" 2>&1); rc=$?
k=$(grep -c "^drift_check: FAIL $work/hm/a-work/t1 " <<<"$out")
[ $rc = 1 ] && [ "$k" = 1 ] && grep -q "sweeping 2 root(s) with drift_check --all: $work/hm/a-work $work/hm/b-work\$" <<<"$out" && grep -q "1 FAIL line(s), 0 ERROR line(s), 2 root(s), 0 excluded" <<<"$out" && ok "N2: overlapping roots are swept once (2 roots, 1 FAIL)" || { bad "N2: overlapping roots (rc=$rc, $k FAIL lines for t1)"; show "$out"; }
# ---- round 4: a glob may have wildcards in its LAST component only. Expanding x/*/mw would stat x/MathWorks/mw
# before any name check could run; a last-component glob only lists its (already name-checked) parent ----
check "a root glob with a wildcard before its last component: refused, rc 2" 2 re:"wildcards are allowed in the last component only" nore:"FAIL $work" -- "$work/x/*/mw" "$work/x/ok"
# ---- round 4: the sweep's lexical walk handles "..", and a relative root is walked from the cwd ----
check "a root with '..' walks back up before following h2/l1 into MathWorks (excluded by the sweep)" 2 re:"^rule_checks_daily: excluded $work/hm/../h2/l1" -- "$work/hm/../h2/l1"
out=$(cd "$work" && bash "$daily" lnk/mwlink -- --now "$NOW" 2>&1); rc=$?
[ $rc = 2 ] && grep -q '^rule_checks_daily: excluded lnk/mwlink' <<<"$out" && ok "a relative root is walked from the cwd (excluded by the sweep)" || { bad "relative root (rc=$rc)"; show "$out"; }
# the sweep itself refuses a literal root that is not a dir (drift_check would too, but later and less plainly)
check "a missing literal root is refused by the sweep itself" 2 re:"^rule_checks_daily: root $work/no-such-root is not a directory" -- "$work/no-such-root" "$work/r1/clean"
# ---- round 5 B2 (Cassie Q17 Q18/Q19 Q20 Q24): '.', the 40-hop limit, '[' mid-path, a double slash ----
mkdir -p "$work/r5x"
check "r5: '.' is skipped, so r5x/./.. walks back before h2/l1 (excluded by the sweep, Q17)" 2 re:"^rule_checks_daily: excluded $work/r5x/./../h2/l1" -- "$work/r5x/./../h2/l1"
for N in 40 41; do
  c="$work/ch$N"; mkdir -p "$c"; i=1; while [ $i -lt $N ]; do ln -s "$c/L$((i+1))" "$c/L$i"; i=$((i+1)); done
  ln -s "$work/hm/a-work" "$c/L$N"
done
check "r5: a 40-hop root chain is followed, not a loop (Q18)" 1 re:"FAIL $work/ch40/L1/t1 " -- "$work/ch40/L1"
check "r5: a 41-hop root chain is a symlink loop (Q19)" 2 re:"root $work/ch41/L1 is a symlink loop" -- "$work/ch41/L1"
check "r5: a '[' wildcard before the last component is refused (Q20)" 2 re:"wildcards are allowed in the last component only" -- "$work/x/[ab]/y*" "$work/x/ok"
check "r5: an empty component before '..' does not shift the walk (r5x//.. then h2/l1: excluded by the sweep, Q24)" 2 re:"^rule_checks_daily: excluded $work/r5x//../h2/l1" -- "$work/r5x//../h2/l1"
out=$(timeout 60 bash "$daily" "$work//r1///clean/" -- --now "$NOW" 2>&1); rc=$?
[ $rc = 0 ] && grep -q "ok   $work//r1///clean/" <<<"$out" && ok "r5: double slashes in a root are walked and finish (Q24)" || { bad "r5: double slashes (rc=$rc)"; show "$out"; }
echo "rule_checks_daily_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
