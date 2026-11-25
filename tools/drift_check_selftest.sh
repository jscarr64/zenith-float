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
# ---- unpushed age: FAIL when the OLDEST unpushed commit is more than 30 min old (exactly 30 min passes) ----
r=$(fresh a1x); commit_at "$r" 1800 c1
                expect "1 ahead, exactly 30 min"           0 run --repo "$r"
r=$(fresh a1y); commit_at "$r" 1801 c1
                expect "1 ahead, 30 min + 1 s"             1 run --repo "$r"
r=$(fresh a1z); commit_at "$r" 1799 c1
                expect "1 ahead, 30 min - 1 s"             0 run --repo "$r"
r=$(fresh am5); commit_at "$r" 600 c1
                expect "10 min old with --max-age-min 5"   1 run --repo "$r" --max-age-min 5
# the OLDEST commit decides, not the newest (Cassie B5a)
r=$(fresh oa1); commit_at "$r" 1801 c1; commit_at "$r" 60 c2
                expect "oldest of 2 is 30 min + 1 s (--max-ahead 2)" 1 run --repo "$r" --max-ahead 2
r=$(fresh oa0); commit_at "$r" 1800 c1; commit_at "$r" 60 c2
                expect "oldest of 2 is exactly 30 min (--max-ahead 2)" 0 run --repo "$r" --max-ahead 2
r=$(fresh dha); git -C "$r" switch -q --detach; commit_at "$r" 1800 d1
                expect "detached HEAD, 1 unpushed exactly 30 min" 0 run --repo "$r"
r=$(fresh dhb); git -C "$r" switch -q --detach; commit_at "$r" 1801 d1
                expect "detached HEAD, 1 unpushed 30 min + 1 s"   1 run --repo "$r"
# ---- dirty work: FAIL when the OLDEST dirty path (modified, deleted, untracked, not ignored) is > 30 min old ----
r=$(fresh dn);  echo x > "$r/new.txt"; touch -d "@$((NOW-60))" "$r/new.txt"
                expect "untracked 1 min old"               0 run --repo "$r"
r=$(fresh dx);  echo x > "$r/new.txt"; touch -d "@$((NOW-1800))" "$r/new.txt"
                expect "untracked exactly 30 min"          0 run --repo "$r"
r=$(fresh do);  echo x > "$r/new.txt"; touch -d "@$((NOW-1801))" "$r/new.txt"
                expect "untracked 30 min + 1 s"            1 run --repo "$r"
r=$(fresh dm);  echo y >> "$r/f.txt"; touch -d "@$((NOW-3600))" "$r/f.txt"
                expect "modified tracked file 1 h old"     1 run --repo "$r"
r=$(fresh di);  echo '*.log' > "$r/.gitignore"; git -C "$r" add .gitignore; GIT_COMMITTER_DATE="@$((NOW-7000)) +0000" git -C "$r" commit -q -m ign
                git -C "$r" push -q 2>/dev/null; echo x > "$r/old.log"; touch -d "@$((NOW-9999))" "$r/old.log"
                expect "ignored file 2.7 h old"            0 run --repo "$r"
r=$(fresh dsub); mkdir -p "$r/d/e"; echo x > "$r/d/e/deep.txt"; touch -d "@$((NOW-4000))" "$r/d/e/deep.txt"
                expect "untracked file in a new subdir, old" 1 run --repo "$r"
r=$(fresh ddn); rm -f "$r/f.txt"
                expect "deleted tracked file, index fresh" 0 run --repo "$r"
r=$(fresh ddo); rm -f "$r/f.txt"; touch -d "@$((NOW-4000))" "$r/.git/index"
                expect "deleted tracked file, age from index 66 min" 1 run --repo "$r"
r=$(fresh drn); git -C "$r" mv f.txt g.txt; touch -d "@$((NOW-60))" "$r/g.txt"; touch -d "@$((NOW-4000))" "$r/.git/index"
                expect "staged rename: age from the new file, not the index" 0 run --repo "$r"
# the OLDEST dirty file decides, not the newest (Cassie B5b)
r=$(fresh d2o); echo x > "$r/a.txt"; echo y > "$r/b.txt"; touch -d "@$((NOW-1801))" "$r/a.txt"; touch -d "@$((NOW-5))" "$r/b.txt"
                expect "two dirty files, oldest 30 min + 1 s" 1 run --repo "$r"
r=$(fresh d2x); echo x > "$r/a.txt"; echo y > "$r/b.txt"; touch -d "@$((NOW-1800))" "$r/a.txt"; touch -d "@$((NOW-5))" "$r/b.txt"
                expect "two dirty files, oldest exactly 30 min" 0 run --repo "$r"
# empty repo (no commit yet): dirty work still counts (Cassie B6)
e="$work/emp1"; git init -q -b w "$e"; echo x > "$e/u.txt"; touch -d "@$((NOW-1800))" "$e/u.txt"
                expect "empty repo, untracked exactly 30 min" 0 run --repo "$e"
e="$work/emp2"; git init -q -b w "$e"; echo x > "$e/u.txt"; touch -d "@$((NOW-1801))" "$e/u.txt"
                expect "empty repo, untracked 30 min + 1 s"   1 run --repo "$e"
# a submodule with old uncommitted changes is dirty work too (Cassie N6 X09)
r=$(fresh dsm); s=$(fresh dsms)
git -C "$r" -c protocol.file.allow=always submodule -q add "$s" sm 2>/dev/null; git -C "$r" commit -q -m sm; git -C "$r" push -q 2>/dev/null
echo z >> "$r/sm/f.txt"; touch -d "@$((NOW-4000))" "$r/sm/f.txt" "$r/sm"
                expect "submodule with 66 min old changes"   1 run --repo "$r"
# ---- pre-commit cap: refuse a new commit while 1 commit is already unpushed (.githooks/pre-commit) ----
r=$(fresh pc0); expect "pre-commit, 0 unpushed"            0 run --pre-commit --repo "$r"
r=$(fresh pc1); commit_at "$r" 10 p1
                expect "pre-commit, 1 unpushed (refuse)"   1 run --pre-commit --repo "$r"
r=$(fresh pcd); git -C "$r" switch -q --detach; commit_at "$r" 10 p1
                expect "pre-commit, detached, 1 unpushed"  1 run --pre-commit --repo "$r"
r=$(fresh pcl); git -C "$r" switch -q -c lone; commit_at "$r" 10 p1
                expect "pre-commit, no upstream, 1 unpushed" 1 run --pre-commit --repo "$r"
r=$(fresh pco); commit_at "$r" 60 c1; commit_at "$r" 50 c2; git -C "$r" switch -q -c other origin/w 2>/dev/null
                expect "pre-commit judges the current branch only" 0 run --pre-commit --repo "$r"
e="$work/pcu"; git init -q -b w "$e"
                expect "pre-commit, unborn branch (first commit of a new repo)" 0 run --pre-commit --repo "$e"
# the real hook
hook="${PRE_COMMIT_HOOK:-$here/../.githooks/pre-commit}"
r=$(fresh hk); mkdir -p "$r/.githooks" "$r/tools"; cp "$hook" "$r/.githooks/pre-commit"; chmod +x "$r/.githooks/pre-commit"; cp "$dc" "$r/tools/drift_check.py"
git -C "$r" config core.hooksPath .githooks
git -C "$r" add .githooks tools; git -C "$r" commit -q -m hooks 2>/dev/null && ok "hook: 1st unpushed commit allowed" || bad "hook: 1st unpushed commit allowed"
echo z >> "$r/f.txt"; git -C "$r" add f.txt
if git -C "$r" commit -q -m second 2>/dev/null; then bad "hook: 2nd unpushed commit refused"; else ok "hook: 2nd unpushed commit refused"; fi
git -C "$r" push -q 2>/dev/null
git -C "$r" commit -q -m second 2>/dev/null && ok "hook: commit allowed again after push" || bad "hook: commit allowed again after push"
e="$work/hku"; git init -q -b w "$e"; mkdir -p "$e/.githooks" "$e/tools"; cp "$hook" "$e/.githooks/pre-commit"; chmod +x "$e/.githooks/pre-commit"; cp "$dc" "$e/tools/drift_check.py"
git -C "$e" config core.hooksPath .githooks; git -C "$e" add -A
git -C "$e" commit -q -m first 2>/dev/null && ok "hook: first commit of a new repo allowed" || bad "hook: first commit of a new repo allowed"
# outside a work tree the hook fails closed with its own message (Cassie N6 X16)
out=$(cd "$work/notgit" && bash "$hook" 2>&1); rc=$?
[ $rc = 2 ] && grep -q 'pre-commit: not in a git work tree' <<<"$out" && ok "hook outside a work tree: rc 2" || { bad "hook outside a work tree (rc=$rc)"; show "$out"; }
echo "drift_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
