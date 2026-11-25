#!/usr/bin/env bash
# Selftest for tools/drift_check.py (Jeff 2026-11-25: push and merge as we go; plan v2 C-DRIFT / C-AHEAD).
# Thresholds under test: FAIL when a branch (or a detached HEAD) is more than 1 commit ahead of what is pushed, or
# when unpushed commits or dirty files are older than 30 minutes (exactly 30 min passes, 30 min + 1 s fails).
# Pre-commit mode refuses a commit when 1 commit is already unpushed (the new one would make it 2).
# Every fixture lives under $ZF_SCRATCH (default $HOME/.cache/zenith-float), never /tmp.
# DRIFT_CHECK can point at another copy of drift_check.py (the mutation runner uses this).
set -u
here="$(cd "$(dirname "$0")" && pwd)"
dc="${DRIFT_CHECK:-$here/drift_check.py}"
root="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/drift-selftest"
mkdir -p "$root" || exit 2
work=$(mktemp -d "$root/run.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
NOW=$(date +%s)                     # one fixed "now" for the whole run; every age is relative to it
export GIT_AUTHOR_NAME=t GIT_AUTHOR_EMAIL=t@t GIT_COMMITTER_NAME=t GIT_COMMITTER_EMAIL=t@t
export GIT_CONFIG_NOSYSTEM=1 GIT_CONFIG_GLOBAL=/dev/null
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "ok   $1"; }
bad() { fail=$((fail+1)); echo "FAIL $1"; }
# expect NAME WANT_RC CMD... : run CMD, compare rc
expect() {
  local name=$1 want=$2; shift 2
  local out rc
  out=$("$@" 2>&1); rc=$?
  if [ "$rc" = "$want" ]; then ok "$name (rc=$rc)"; else bad "$name (want rc=$want, got rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; fi
}
# commit_at REPO AGE_SECONDS MSG : one commit whose committer date is NOW-AGE
commit_at() {
  local r=$1 age=$2 msg=$3 ts=$((NOW-$2))
  echo "$msg" >> "$r/f.txt"
  git -C "$r" add f.txt
  GIT_COMMITTER_DATE="@$ts +0000" GIT_AUTHOR_DATE="@$ts +0000" git -C "$r" commit -q --no-verify -m "$msg"
}
# fresh NAME : origin bare repo + clone with one pushed commit on branch w (upstream set)
fresh() {
  local d="$work/$1"
  git init -q --bare "$d/origin.git"
  git init -q -b w "$d/r"
  git -C "$d/r" remote add origin "$d/origin.git"
  commit_at "$d/r" 7200 base
  git -C "$d/r" push -q -u origin w 2>/dev/null
  echo "$d/r"
}
run() { python3 -B "$dc" --now "$NOW" "$@"; }

[ -f "$dc" ] || { echo "drift_check_selftest: $dc missing"; echo "drift_check_selftest: FAIL (0 ok, 1 failed)"; exit 1; }

r=$(fresh a0);  expect "0 ahead, clean"                    0 run --repo "$r"
r=$(fresh a1);  commit_at "$r" 600 c1
                expect "1 ahead, 10 min old"               0 run --repo "$r"
r=$(fresh a2);  commit_at "$r" 600 c1; commit_at "$r" 300 c2
                expect "2 ahead, fresh"                    1 run --repo "$r"
r=$(fresh a1x); commit_at "$r" 1800 c1
                expect "1 ahead, exactly 30 min"           0 run --repo "$r"
r=$(fresh a1y); commit_at "$r" 1801 c1
                expect "1 ahead, 30 min + 1 s"             1 run --repo "$r"
r=$(fresh a1z); commit_at "$r" 1799 c1
                expect "1 ahead, 30 min - 1 s"             0 run --repo "$r"
# custom thresholds are honoured both ways
r=$(fresh t2);  commit_at "$r" 600 c1; commit_at "$r" 300 c2
                expect "2 ahead with --max-ahead 2"        0 run --repo "$r" --max-ahead 2
r=$(fresh t3);  commit_at "$r" 600 c1
                expect "10 min old with --max-age-min 5"   1 run --repo "$r" --max-age-min 5
# no upstream: a local branch never pushed counts commits not on any remote ref
r=$(fresh nu1); git -C "$r" switch -q -c lone; commit_at "$r" 60 l1
                expect "no upstream, 1 commit"             0 run --repo "$r"
r=$(fresh nu2); git -C "$r" switch -q -c lone; commit_at "$r" 60 l1; commit_at "$r" 30 l2
                expect "no upstream, 2 commits"            1 run --repo "$r"
# upstream deleted on the remote: falls back to "not on any remote ref"
r=$(fresh ug);  git -C "$r" switch -q -c gone; commit_at "$r" 60 g1; git -C "$r" push -q -u origin gone 2>/dev/null
                commit_at "$r" 30 g2; commit_at "$r" 20 g3
                git -C "$r" push -q origin --delete gone 2>/dev/null; git -C "$r" fetch -q --prune origin
                expect "upstream gone, 3 unpushed"         1 run --repo "$r"
# a non-current local branch is checked too
r=$(fresh oth); git -C "$r" switch -q -c side; commit_at "$r" 60 s1; commit_at "$r" 50 s2; git -C "$r" switch -q w
                expect "other local branch 2 ahead"        1 run --repo "$r"
# detached HEAD
r=$(fresh dh1); git -C "$r" switch -q --detach; commit_at "$r" 60 d1
                expect "detached HEAD, 1 unpushed"         0 run --repo "$r"
r=$(fresh dh2); git -C "$r" switch -q --detach; commit_at "$r" 60 d1; commit_at "$r" 50 d2
                expect "detached HEAD, 2 unpushed"         1 run --repo "$r"
r=$(fresh dh0); git -C "$r" switch -q --detach
                expect "detached HEAD at pushed commit"    0 run --repo "$r"
# dirty work: oldest dirty file mtime
r=$(fresh dn);  echo x > "$r/new.txt"; touch -d "@$((NOW-60))" "$r/new.txt"
                expect "untracked 1 min old"               0 run --repo "$r"
r=$(fresh dx);  echo x > "$r/new.txt"; touch -d "@$((NOW-1800))" "$r/new.txt"
                expect "untracked exactly 30 min"          0 run --repo "$r"
r=$(fresh do);  echo x > "$r/new.txt"; touch -d "@$((NOW-1801))" "$r/new.txt"
                expect "untracked 30 min + 1 s"            1 run --repo "$r"
r=$(fresh dm);  echo y >> "$r/f.txt"; touch -d "@$((NOW-3600))" "$r/f.txt"
                expect "modified tracked file 1 h old"     1 run --repo "$r"
r=$(fresh di);  echo '*.log' > "$r/.gitignore"; git -C "$r" add .gitignore; GIT_COMMITTER_DATE="@$((NOW-7000)) +0000" git -C "$r" commit -q --no-verify -m ign
                git -C "$r" push -q 2>/dev/null; echo x > "$r/old.log"; touch -d "@$((NOW-9999))" "$r/old.log"
                expect "ignored file 2.7 h old"            0 run --repo "$r"
r=$(fresh dsub); mkdir -p "$r/d/e"; echo x > "$r/d/e/deep.txt"; touch -d "@$((NOW-4000))" "$r/d/e/deep.txt"
                expect "untracked file in a new subdir, old" 1 run --repo "$r"
r=$(fresh ddn); rm -f "$r/f.txt"
                expect "deleted tracked file, index fresh" 0 run --repo "$r"
r=$(fresh ddo); rm -f "$r/f.txt"; touch -d "@$((NOW-4000))" "$r/.git/index"
                expect "deleted tracked file, age from index 66 min" 1 run --repo "$r"
# pre-commit mode
r=$(fresh pc0); expect "pre-commit, 0 unpushed"            0 run --pre-commit --repo "$r"
r=$(fresh pc1); commit_at "$r" 10 p1
                expect "pre-commit, 1 unpushed (refuse)"   1 run --pre-commit --repo "$r"
r=$(fresh pcd); git -C "$r" switch -q --detach; commit_at "$r" 10 p1
                expect "pre-commit, detached, 1 unpushed"  1 run --pre-commit --repo "$r"
# the real hook: .githooks/pre-commit refuses the 2nd unpushed commit and allows it after a push
hook="${PRE_COMMIT_HOOK:-$here/../.githooks/pre-commit}"
if [ -f "$hook" ]; then
  r=$(fresh hk); mkdir -p "$r/.githooks" "$r/tools"; cp "$hook" "$r/.githooks/pre-commit"; cp "$dc" "$r/tools/drift_check.py"
  git -C "$r" config core.hooksPath .githooks
  git -C "$r" add .githooks tools; git -C "$r" commit -q -m hooks && ok "hook: 1st unpushed commit allowed" || bad "hook: 1st unpushed commit allowed"
  echo z >> "$r/f.txt"; git -C "$r" add f.txt
  if git -C "$r" commit -q -m second 2>/dev/null; then bad "hook: 2nd unpushed commit refused"; else ok "hook: 2nd unpushed commit refused"; fi
  git -C "$r" push -q 2>/dev/null
  git -C "$r" commit -q -m second 2>/dev/null && ok "hook: commit allowed again after push" || bad "hook: commit allowed again after push"
else bad "hook: $hook missing"; fi
# scan mode: one good and one bad repo under a directory; the bad one is named
s="$work/scan"; mkdir -p "$s"
g=$(fresh sg); b=$(fresh sb); commit_at "$b" 10 b1; commit_at "$b" 5 b2
mv "$(dirname "$g")" "$s/good"; mv "$(dirname "$b")" "$s/bad"
out=$(run --scan "$s" 2>&1); rc=$?
if [ $rc = 1 ] && grep -q "bad/r" <<<"$out" && ! grep -q "FAIL.*good/r" <<<"$out"; then ok "scan: bad repo named, good not (rc=1)"; else bad "scan (rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; fi
expect "scan: only the good repo" 0 run --scan "$s/good"
mkdir -p "$work/emptydir"; expect "scan: empty dir" 0 run --scan "$work/emptydir"
# linked worktrees of one repo are reported once per worktree, branches once per repo
r=$(fresh wt); git -C "$r" worktree add -q -b wtb "$work/wt-linked" 2>/dev/null; commit_at "$work/wt-linked" 30 x1; commit_at "$work/wt-linked" 20 x2
                expect "linked worktree branch 2 ahead"    1 run --repo "$work/wt-linked"
# read-only: the check must not touch the index (it runs git --no-optional-locks)
r=$(fresh ro); echo y >> "$r/f.txt"; touch -d "@$((NOW-60))" "$r/f.txt"; before=$(stat -c %Y.%s "$r/.git/index")
sleep 1; run --repo "$r" >/dev/null 2>&1; after=$(stat -c %Y.%s "$r/.git/index")
[ "$before" = "$after" ] && ok "index untouched" || bad "index untouched ($before -> $after)"
# usage errors fail closed with rc 2
expect "not a git repo" 2 run --repo "$work"
expect "no target given" 2 run
echo "drift_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
