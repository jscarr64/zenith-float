#!/usr/bin/env bash
# Selftest for tools/drift_check.py and .githooks/pre-commit (Jeff 2026-11-25: push and merge as we go; plan v2
# C-DRIFT / C-AHEAD). Each section below landed as its own red commit before the code that makes it pass.
# Every fixture lives under $ZF_SCRATCH (default $HOME/.cache/zenith-float), never /tmp.
# DRIFT_CHECK / PRE_COMMIT_HOOK point at other copies (the mutation runner uses this).
set -u
here="$(cd "$(dirname "$0")" && pwd)"
dc="${DRIFT_CHECK:-$here/drift_check.py}"
root="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/drift-selftest"
mkdir -p "$root" || exit 2
work=$(mktemp -d "$root/run.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
NOW=$(date +%s)                     # one fixed "now" for the whole run; every age is relative to it
# fixture identity and config, local to this process (never export these in a shell that makes real commits)
export GIT_AUTHOR_NAME=t GIT_AUTHOR_EMAIL=t@t GIT_COMMITTER_NAME=t GIT_COMMITTER_EMAIL=t@t
export GIT_CONFIG_NOSYSTEM=1 GIT_CONFIG_GLOBAL=/dev/null GIT_CEILING_DIRECTORIES="$work"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "ok   $1"; }
bad() { fail=$((fail+1)); echo "FAIL $1"; }
show() { printf '%s\n' "$1" | sed 's/^/     | /'; }
# expect NAME WANT_RC CMD... : run CMD, compare rc
expect() {
  local name=$1 want=$2; shift 2
  local out rc
  out=$("$@" 2>&1); rc=$?
  if [ "$rc" = "$want" ]; then ok "$name (rc=$rc)"; else bad "$name (want rc=$want, got rc=$rc)"; show "$out"; fi
}
# commit_at REPO AGE_SECONDS MSG : one commit whose committer date is NOW-AGE
commit_at() {
  local r=$1 age=$2 msg=$3 ts=$((NOW-$2))
  echo "$msg" >> "$r/f.txt"
  git -C "$r" add f.txt
  GIT_COMMITTER_DATE="@$ts +0000" GIT_AUTHOR_DATE="@$ts +0000" git -C "$r" commit -q -m "$msg"
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
# ---- ahead cap: FAIL when a branch or detached HEAD is more than 1 commit ahead of what is pushed ----
r=$(fresh a0);  expect "0 ahead, clean"                    0 run --repo "$r"
r=$(fresh a1);  commit_at "$r" 600 c1
                expect "1 ahead"                           0 run --repo "$r"
r=$(fresh a2);  commit_at "$r" 600 c1; commit_at "$r" 300 c2
                expect "2 ahead"                           1 run --repo "$r"
r=$(fresh t2);  commit_at "$r" 600 c1; commit_at "$r" 300 c2
                expect "2 ahead with --max-ahead 2"        0 run --repo "$r" --max-ahead 2
r=$(fresh t3);  commit_at "$r" 600 c1; commit_at "$r" 500 c2; commit_at "$r" 300 c3
                expect "3 ahead with --max-ahead 2"        1 run --repo "$r" --max-ahead 2
# no upstream: a local branch never pushed counts commits not on any remote ref
r=$(fresh nu1); git -C "$r" switch -q -c lone; commit_at "$r" 60 l1
                expect "no upstream, 1 commit"             0 run --repo "$r"
r=$(fresh nu2); git -C "$r" switch -q -c lone; commit_at "$r" 60 l1; commit_at "$r" 30 l2
                expect "no upstream, 2 commits"            1 run --repo "$r"
# upstream deleted on the remote: falls back to "not on any remote ref"
r=$(fresh ug);  git -C "$r" switch -q -c gone; commit_at "$r" 60 g1; git -C "$r" push -q -u origin gone 2>/dev/null
                commit_at "$r" 30 g2; commit_at "$r" 20 g3
                git -C "$r" push -q origin --delete gone 2>/dev/null; git -C "$r" fetch -q --prune origin
                expect "upstream gone, 2 unpushed"         1 run --repo "$r"
# upstream is another LOCAL branch (branch.feat.remote=.): count against the remotes, not against w (Cassie B4)
r=$(fresh lup); commit_at "$r" 60 w1; git -C "$r" switch -q -c feat --track w 2>/dev/null; commit_at "$r" 30 f1
                expect "local upstream: 2 commits on no remote ref" 1 run --repo "$r"
r=$(fresh lu1); git -C "$r" switch -q -c feat --track w 2>/dev/null; commit_at "$r" 30 f1
                expect "local upstream: 1 commit on no remote ref"  0 run --repo "$r"
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
# detached on a branch tip that is 2 ahead: reported once, as the branch (Cassie N6 X08)
r=$(fresh dtp); git -C "$r" switch -q -c side; commit_at "$r" 60 s1; commit_at "$r" 50 s2; git -C "$r" switch -q --detach side
out=$(run --repo "$r" 2>&1); rc=$?; k=$(grep -c '^drift_check: FAIL /' <<<"$out")
[ $rc = 1 ] && [ "$k" = 1 ] && ok "detached on a 2-ahead branch tip: reported once" || { bad "detached on a 2-ahead branch tip: reported once (rc=$rc, $k reports)"; show "$out"; }
# two --repo: a bad first repo is not hidden by a good second one
b=$(fresh mb); commit_at "$b" 60 x1; commit_at "$b" 50 x2; g=$(fresh mg)
                expect "bad repo then good repo"           1 run --repo "$b" --repo "$g"
# a GIT_DIR exported by a hook (git does this) must not redirect the check to another repo
r=$(fresh ge);  commit_at "$r" 60 e1; commit_at "$r" 50 e2; o=$(fresh geo)
                expect "GIT_DIR of another repo in the env is ignored" 1 env GIT_DIR="$o/.git" python3 -B "$dc" --now "$NOW" --repo "$r"
# read-only: the check must not touch the index (it runs git --no-optional-locks)
r=$(fresh ro); echo y >> "$r/f.txt"; touch -d "@$((NOW-60))" "$r/f.txt"; before=$(stat -c %Y.%s "$r/.git/index")
sleep 1; run --repo "$r" >/dev/null 2>&1; after=$(stat -c %Y.%s "$r/.git/index")
[ "$before" = "$after" ] && ok "index untouched" || bad "index untouched ($before -> $after)"
# usage and git errors fail closed with rc 2 (Cassie B7)
mkdir -p "$work/notgit"
expect "not a git repo"                                  2 run --repo "$work/notgit"
r=$(fresh sd); mkdir -p "$r/sub"
expect "a subdir of a work tree, not its top"            2 run --repo "$r/sub"
expect "no target given"                                 2 run
r=$(fresh gx); commit_at "$r" 60 x1; commit_at "$r" 50 x2; h=$(git -C "$r" rev-parse HEAD); rm -f "$r/.git/objects/${h:0:2}/${h:2}"
expect "missing commit object: git error"                2 run --repo "$r"
echo "drift_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
