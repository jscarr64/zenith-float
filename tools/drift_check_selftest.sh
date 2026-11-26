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
# ---- --scan DIR: DIR and its direct children; a child with a .git it cannot read fails closed (Cassie B3) ----
s="$work/scan"; mkdir -p "$s"
g=$(fresh sg); b=$(fresh sb); commit_at "$b" 10 b1; commit_at "$b" 5 b2
mv "$g" "$s/good"; mv "$b" "$s/bad"; mkdir -p "$s/plain"; git init -q --bare "$s/bare.git"
out=$(run --all --scan "$s" 2>&1); rc=$?
if [ $rc = 1 ] && grep -q "FAIL $s/bad" <<<"$out" && ! grep -q "FAIL $s/good" <<<"$out"; then ok "scan: bad repo named, good not, plain dir and bare repo skipped (rc=1)"; else bad "scan (rc=$rc)"; show "$out"; fi
expect "scan: only the good repo"                        0 run --all --scan "$s/good"
mkdir -p "$work/emptydir"; expect "scan: empty dir"      0 run --all --scan "$work/emptydir"
expect "scan: not a directory (Cassie B7 X05)"           2 run --all --scan "$work/no-such-dir"
s="$work/scan-broken"; mkdir -p "$s/broken"; echo "gitdir: $work/no-such-gitdir" > "$s/broken/.git"
out=$(run --all --scan "$s" 2>&1); rc=$?
[ $rc = 2 ] && grep -q "$s/broken" <<<"$out" && ok "scan: child whose .git points at a missing gitdir: rc 2, named" || { bad "scan: broken .git child (rc=$rc)"; show "$out"; }
s="$work/scan-head"; mkdir -p "$s"; c=$(fresh sh); commit_at "$c" 60 h1; commit_at "$c" 50 h2; mv "$c" "$s/corrupt"; echo garbage > "$s/corrupt/.git/HEAD"
out=$(run --all --scan "$s" 2>&1); rc=$?
[ $rc = 2 ] && grep -q "$s/corrupt" <<<"$out" && ok "scan: child with a corrupt HEAD and 2 unpushed: rc 2, named" || { bad "scan: corrupt HEAD child (rc=$rc)"; show "$out"; }
# linked worktrees share branches: a bad branch is reported once, not once per worktree
r=$(fresh wt); git -C "$r" switch -q -c side2; commit_at "$r" 60 a1; commit_at "$r" 50 a2; git -C "$r" switch -q w
git -C "$r" worktree add -q --detach "$r-l1" w 2>/dev/null; git -C "$r" worktree add -q --detach "$r-l2" w 2>/dev/null
out=$(run --repo "$r" --repo "$r-l1" --repo "$r-l2" 2>&1); rc=$?
k=$(grep -c 'side2: 2 commits ahead' <<<"$out")
[ $rc = 1 ] && [ "$k" = 1 ] && ok "branch shared by 3 worktrees reported once" || { bad "branch shared by 3 worktrees reported once (rc=$rc, reported $k times)"; show "$out"; }
out=$(run --all --scan "$work/wt" 2>&1); rc=$?
k=$(grep -c 'side2: 2 commits ahead' <<<"$out")
[ $rc = 1 ] && [ "$k" = 1 ] && ok "scan of 3 linked worktrees: branch reported once" || { bad "scan of 3 linked worktrees (rc=$rc, reported $k times)"; show "$out"; }
r=$(fresh wtl); git -C "$r" worktree add -q -b wtb "$work/wtl-linked" 2>/dev/null; commit_at "$work/wtl-linked" 30 x1; commit_at "$work/wtl-linked" 20 x2
                expect "linked worktree branch 2 ahead"    1 run --repo "$work/wtl-linked"
# ---- active mode (Jeff 2026-11-25 3:55 PM CT: "limit the check to active work trees") ----
# --scan checks only ACTIVE work trees: newest change (HEAD commit time, index mtime, logs/HEAD mtime, newest dirty
# file mtime) within the window (default 24 h; exactly 24 h is active), or dirty work less than 30 min old, or listed
# with --include / --include-file. Unknown activity counts as active (fail closed). --all checks every work tree.
# Inside an active tree: the HEAD branch, plus other branches whose tip changed within the window.
W=86400
# age REPO SECS : index and logs/HEAD last touched SECS ago (git itself touched them at real now)
age() { touch -d "@$((NOW-$2))" "$1/.git/index" "$1/.git/logs/HEAD"; }
# stale NAME COMMIT_AGE : repo 2 commits ahead, both commits COMMIT_AGE old, index/logs/HEAD 2 days old
stale() { local r; r=$(fresh "$1"); commit_at "$r" "$2" s1; commit_at "$r" "$2" s2; age "$r" 172800; echo "$r"; }
r=$(stale ina 172800)
expect "active: tree idle 2 days is skipped"                 0 run --scan "$work/ina"
expect "active: --all checks the idle tree"                  1 run --all --scan "$work/ina"
expect "active: --include checks the idle tree"              1 run --scan "$work/ina" --include "$r"
printf '# comment\n%s\n' "$r" > "$work/ina-inc.txt"
expect "active: --include-file checks the idle tree"         1 run --scan "$work/ina" --include-file "$work/ina-inc.txt"
out=$(run --scan "$work/ina" 2>&1); grep -q "skip $r" <<<"$out" && ok "active: skipped tree is named" || { bad "active: skipped tree is named"; show "$out"; }
r=$(stale cwe $W)
expect "active: HEAD commit exactly 24 h old is active"      1 run --scan "$work/cwe"
r=$(stale cwo $((W+1)))
expect "active: HEAD commit 24 h + 1 s old is idle"          0 run --scan "$work/cwo"
r=$(stale ixr 172800); touch -d "@$((NOW-100))" "$r/.git/index"
expect "active: recent index alone makes it active"          1 run --scan "$work/ixr"
r=$(stale rlr 172800); touch -d "@$((NOW-100))" "$r/.git/logs/HEAD"
expect "active: recent logs/HEAD alone makes it active"      1 run --scan "$work/rlr"
r=$(stale dwr 172800); echo u > "$r/u.txt"; touch -d "@$((NOW-7200))" "$r/u.txt"
expect "active: dirty file 2 h old (inside window) is active" 1 run --scan "$work/dwr"
r=$(stale d29 172800); echo u > "$r/u.txt"; touch -d "@$((NOW-1799))" "$r/u.txt"
expect "active: dirty < 30 min is active even with a 5 min window" 1 run --scan "$work/d29" --active-window-min 5
r=$(stale d30 172800); echo u > "$r/u.txt"; touch -d "@$((NOW-1800))" "$r/u.txt"
expect "active: dirty exactly 30 min, 5 min window: idle"    0 run --scan "$work/d30" --active-window-min 5
r=$(stale brk 172800); echo 0123456789abcdef0123456789abcdef01234567 > "$r/.git/refs/heads/w"
out=$(run --scan "$work/brk" 2>&1); rc=$?
[ $rc != 0 ] && ok "active: unreadable activity fails closed (checked, rc=$rc)" || { bad "active: unreadable activity fails closed (rc=0)"; show "$out"; }
# branches inside an active tree
r=$(fresh obr); git -C "$r" switch -q -c old; commit_at "$r" 172800 o1; commit_at "$r" 172800 o2; git -C "$r" switch -q w
expect "active: other branch idle 2 days is skipped"         0 run --scan "$work/obr"
expect "active: --all checks the idle branch"                1 run --all --scan "$work/obr"
r=$(fresh obe); git -C "$r" switch -q -c old; commit_at "$r" 172800 o1; commit_at "$r" $W o2; git -C "$r" switch -q w
expect "active: other branch tip exactly 24 h old is checked" 1 run --scan "$work/obe"
r=$(fresh obo); git -C "$r" switch -q -c old; commit_at "$r" 172800 o1; commit_at "$r" $((W+1)) o2; git -C "$r" switch -q w
expect "active: other branch tip 24 h + 1 s old is skipped"  0 run --scan "$work/obo"
r=$(fresh obk); echo 0123456789abcdef0123456789abcdef01234567 > "$r/.git/refs/heads/bogus"
expect "active: other branch with unreadable tip is checked (git error, rc 2)" 2 run --scan "$work/obk"
r=$(stale ind 172800); out=$(run --scan "$work/ind" --include "$r" 2>&1)
! grep -q "skip $r" <<<"$out" && ok "active: an --include tree that is also scanned is checked once, not skipped" || { bad "active: --include tree also skipped"; show "$out"; }
expect "active: --repo always checks every branch"           1 run --repo "$r"
r=$(fresh obt); git -C "$r" hash-object -w --stdin <<<blob > "$r/.git/refs/heads/blobref"
out=$(run --scan "$work/obt" 2>&1); rc=$?
[ $rc = 0 ] && ! grep -q "skip.*blobref" <<<"$out" && ok "active: branch with no commit date (points at a blob) is checked, not skipped" || { bad "active: branch with no commit date (rc=$rc)"; show "$out"; }
ra=$(stale ia 172800); rb=$(stale ib 172800); echo "$ra" > "$work/ia-inc.txt"
out=$(run --scan "$work/ia" --include "$rb" --include-file "$work/ia-inc.txt" 2>&1)
grep -q "FAIL $ra " <<<"$out" && grep -q "FAIL $rb " <<<"$out" && ok "active: --include and --include-file together both count" || { bad "active: --include and --include-file together"; show "$out"; }
expect "active: bad window value"                            2 run --scan "$work/obo" --active-window-min -1
e="$work/act-empty"; mkdir -p "$e"; git init -q -b w "$e/r"; echo x > "$e/r/u.txt"; touch -d "@$((NOW-1801))" "$e/r/u.txt"
expect "active: empty repo (no HEAD commit) counts as active, old dirty file fails" 1 run --scan "$e"
# ---- round 2 (Cassie R2-B1): every unreadable scan entry fails closed; the scan dir itself is checked ----
# a scan child that cannot be read (mode 000) and holds 2 unpushed commits: rc 2 and named, never skipped (R2-B1a)
s="$work/scan-perm"; mkdir -p "$s"; c=$(fresh sp); commit_at "$c" 60 p1; commit_at "$c" 50 p2; mv "$c" "$s/locked"; chmod 000 "$s/locked"
out=$(run --all --scan "$s" 2>&1); rc=$?; chmod 755 "$s/locked"
[ $rc = 2 ] && grep -q "$s/locked" <<<"$out" && ok "scan: unreadable child (mode 000) with 2 unpushed: rc 2, named" || { bad "scan: unreadable child (mode 000) (want rc=2 and named, got rc=$rc)"; show "$out"; }
# an unreadable --scan DIR: rc 2 with a message, no traceback (R2-B1b)
s="$work/scan-perm-dir"; mkdir -p "$s/x"; chmod 000 "$s"
out=$(run --all --scan "$s" 2>&1); rc=$?; chmod 755 "$s"
[ $rc = 2 ] && ! grep -q Traceback <<<"$out" && grep -q "$s" <<<"$out" && ok "scan: unreadable --scan dir: rc 2, named, no traceback" || { bad "scan: unreadable --scan dir (want rc=2 without a traceback, got rc=$rc)"; show "$out"; }
# the scan DIR itself is checked: --scan of a repo with 2 unpushed commits fails (Cassie Y16)
r=$(fresh sself); commit_at "$r" 60 q1; commit_at "$r" 50 q2
expect "scan: DIR itself is a repo with 2 unpushed"          1 run --scan "$r"
# --include-file that cannot be read: rc 2 (Cassie Y05)
expect "active: --include-file missing"                      2 run --scan "$work/emptydir" --include-file "$work/no-such-include.txt"
# activity that cannot be parsed (HEAD commit time not a number, via a git stub) counts as active: checked, ok (Cassie Y01)
realgit=$(command -v git); mkdir -p "$work/gitstub"
cat > "$work/gitstub/git" <<STUB
#!/usr/bin/env bash
case " \$* " in *" log -1 --format=%ct HEAD "*) [ -n "\${STUB_CT_BAD:-}" ] && { echo not-a-number; exit 0; };; esac
exec "$realgit" "\$@"
STUB
chmod +x "$work/gitstub/git"
r=$(stale yct 172800); git -C "$r" push -q -f origin w 2>/dev/null
out=$(STUB_CT_BAD=1 PATH="$work/gitstub:$PATH" python3 -B "$dc" --now "$NOW" --scan "$work/yct" 2>&1); rc=$?
[ $rc = 0 ] && grep -q "ok   $r" <<<"$out" && ! grep -q "skip $r" <<<"$out" && ok "active: HEAD time not a number: activity unknown, checked, ok (rc=0)" || { bad "active: HEAD time not a number (want rc=0 and checked, got rc=$rc)"; show "$out"; }
# ---- round 2 (Cassie R2-B2): numeric options must be finite and >= 0, else rc 2 ----
r=$(fresh nmo); commit_at "$r" 7000 n1          # 1 unpushed commit, ~2 h old: a real age violation
for v in nan inf -inf -1 -0.0001; do
  expect "--max-age-min $v is a usage error"                 2 run --repo "$r" --max-age-min "$v"
done
c=$(fresh nmc)
expect "--max-age-min 0 is valid (clean tree)"               0 run --repo "$c" --max-age-min 0
for v in nan inf -inf; do
  expect "--active-window-min $v is a usage error"           2 run --scan "$work/obo" --active-window-min "$v"
done
expect "--active-window-min 0 is valid (clean tree, Cassie Y10)" 0 run --repo "$c" --active-window-min 0
expect "--active-window-min -0.0001 is a usage error"        2 run --repo "$c" --active-window-min -0.0001
# ---- round 2 (Cassie R2-B2d): dirty work older than max-age is a violation, so it is never skipped as idle ----
r=$(stale d31 172800); echo u > "$r/u.txt"; touch -d "@$((NOW-1801))" "$r/u.txt"
expect "active: dirty 30 min + 1 s, 5 min window: checked, FAIL" 1 run --scan "$work/d31" --active-window-min 5
r=$(stale dold 172800); echo u > "$r/u.txt"; touch -d "@$((NOW-172800))" "$r/u.txt"
expect "active: dirty 2 days old in a tree idle 2 days: checked, FAIL" 1 run --scan "$work/dold"
r="$work/d30/r"; out=$(run --scan "$work/d30" --active-window-min 5 2>&1)
grep -q "skip $r" <<<"$out" && ok "active: dirty exactly 30 min, 5 min window: still skipped (named)" || { bad "active: dirty exactly 30 min should still be skipped"; show "$out"; }
# ---- round 2 sweep: decisions the round-2 fixes made reachable only through these cases (A06, X05) ----
r=$(stale dx30 172800); git -C "$r" push -q origin w 2>/dev/null; touch -d "@$((NOW-172800))" "$r/.git/index" "$r/.git/logs/HEAD"
echo u > "$r/u.txt"; touch -d "@$((NOW-1800))" "$r/u.txt"
out=$(run --scan "$work/dx30" 2>&1); rc=$?
[ $rc = 0 ] && grep -q "ok   $r" <<<"$out" && ! grep -q "skip $r" <<<"$out" && ok "active: dirty exactly 30 min in a tree idle 2 days: inside the window, checked (A06)" || { bad "active: dirty exactly 30 min, default window (rc=$rc)"; show "$out"; }
out=$(run --all --scan "$work/no-such-dir" 2>&1); rc=$?
[ $rc = 2 ] && grep -q "no-such-dir is not a directory" <<<"$out" && ok "scan: a missing dir is named as not a directory (X05)" || { bad "scan: missing dir message (rc=$rc)"; show "$out"; }
# ---- round 3 (Cassie W05): dirty files of two ages in an idle tree with a small window: the fresh one makes it active ----
r=$(stale w05 172800); git -C "$r" push -q origin w 2>/dev/null; touch -d "@$((NOW-172800))" "$r/.git/index" "$r/.git/logs/HEAD"
echo a > "$r/u1.txt"; echo b > "$r/u2.txt"; touch -d "@$((NOW-60))" "$r/u1.txt"; touch -d "@$((NOW-1801))" "$r/u2.txt"
out=$(run --scan "$work/w05" --active-window-min 5 2>&1); rc=$?
[ $rc = 1 ] && ! grep -q "skip $r" <<<"$out" && ok "active: dirty 60 s and 1801 s in an idle tree, 5 min window: checked, FAIL (W05)" || { bad "active: dirty 60 s + 1801 s, 5 min window (want rc=1, not skipped, got rc=$rc)"; show "$out"; }
# ---- round 3 (Jeff's MathWorks rule, Cassie N3): --exclude NAME leaves out scanned entries with that path component ----
s="$work/xs"; mkdir -p "$s"; m=$(fresh xmw); commit_at "$m" 60 m1; commit_at "$m" 50 m2; mv "$m" "$s/MathWorks"; k=$(fresh xok); mv "$k" "$s/ok"
expect "exclude: control, a MathWorks child is scanned without --exclude" 1 run --all --scan "$s"
out=$(run --all --scan "$s" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $s/MathWorks" <<<"$out" && ! grep -q "FAIL $s/MathWorks" <<<"$out" && grep -q "ok   $s/ok" <<<"$out" && ok "exclude: a MathWorks child is left out and named, the other child checked (rc=0)" || { bad "exclude: MathWorks child (want rc=0 and excluded, got rc=$rc)"; show "$out"; }
t="$work/xt"; mkdir -p "$t/MathWorks"; m=$(fresh xmw2); commit_at "$m" 60 m1; commit_at "$m" 50 m2; mv "$m" "$t/MathWorks/r"; mkdir -p "$work/xl"; ln -s "$t/MathWorks/r" "$work/xl/link"
out=$(run --all --scan "$work/xl" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $work/xl/link" <<<"$out" && ok "exclude: a child that is a symlink into MathWorks is left out (rc=0)" || { bad "exclude: symlink child (want rc=0, got rc=$rc)"; show "$out"; }
out=$(run --all --scan "$t/MathWorks/r" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $t/MathWorks/r" <<<"$out" && ! grep -q FAIL <<<"$out" && ok "exclude: a scan DIR inside MathWorks is left out itself (rc=0)" || { bad "exclude: scan DIR inside MathWorks (want rc=0, got rc=$rc)"; show "$out"; }
u="$work/xu"; mkdir -p "$u"; m=$(fresh xmo); commit_at "$m" 60 m1; commit_at "$m" 50 m2; mv "$m" "$u/MathWorksOld"
expect "exclude: a name that only contains MathWorks is still checked" 1 run --all --scan "$u" --exclude MathWorks
v="$work/xv"; mkdir -p "$v"; m=$(fresh xa); commit_at "$m" 60 m1; commit_at "$m" 50 m2; mv "$m" "$v/a"; m=$(fresh xb); commit_at "$m" 60 m1; commit_at "$m" 50 m2; mv "$m" "$v/b"
out=$(run --all --scan "$v" --exclude a --exclude b 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $v/a" <<<"$out" && grep -q "excluded $v/b" <<<"$out" && ok "exclude: repeatable (rc=0)" || { bad "exclude: repeatable (want rc=0, got rc=$rc)"; show "$out"; }
expect "exclude: --repo is never excluded"                   1 run --repo "$t/MathWorks/r" --exclude MathWorks
# ---- round 4 (Cassie X01, X05): the name decides as given, and matching is exact (case-sensitive) ----
s="$work/x01"; mkdir -p "$s"; m=$(fresh x01r); commit_at "$m" 60 m1; commit_at "$m" 50 m2; ln -s "$m" "$s/MathWorks"
out=$(run --all --scan "$s" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $s/MathWorks" <<<"$out" && ok "exclude: a symlink NAMED MathWorks that points elsewhere is left out (X01)" || { bad "exclude: symlink named MathWorks (want rc=0, got rc=$rc)"; show "$out"; }
s="$work/x05"; mkdir -p "$s"; for n in mathworks MATHWORKS; do m=$(fresh "x05$n"); commit_at "$m" 60 m1; commit_at "$m" 50 m2; mv "$m" "$s/$n"; done
out=$(run --all --scan "$s" --exclude MathWorks 2>&1); rc=$?
[ $rc = 1 ] && grep -q "FAIL $s/mathworks " <<<"$out" && grep -q "FAIL $s/MATHWORKS " <<<"$out" && ok "exclude: mathworks and MATHWORKS are still checked (exact match, X05)" || { bad "exclude: case variants (want rc=1 and both checked, got rc=$rc)"; show "$out"; }
# ---- round 4 B1 (Cassie): an excluded --scan folder is decided by name BEFORE it is stat'ed or listed ----
# an unlistable (mode 000) excluded folder, and one whose parent cannot be searched, must still just be "excluded"
t="$work/b1"; mkdir -p "$t/a/MathWorks/locked" "$t/b/MathWorks/r"; chmod 000 "$t/a/MathWorks/locked" "$t/b/MathWorks"
out=$(run --all --scan "$t/a/MathWorks/locked" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $t/a/MathWorks/locked" <<<"$out" && ok "B1: an unlistable excluded --scan folder is not listed (rc 0, excluded)" || { bad "B1: unlistable excluded --scan folder (want rc=0 excluded, got rc=$rc)"; show "$out"; }
out=$(run --all --scan "$t/b/MathWorks/r" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $t/b/MathWorks/r" <<<"$out" && ok "B1: an excluded --scan folder is not stat'ed (unsearchable parent: rc 0, excluded)" || { bad "B1: excluded --scan folder stat'ed (want rc=0 excluded, got rc=$rc)"; show "$out"; }
chmod 755 "$t/a/MathWorks/locked" "$t/b/MathWorks"
# and, when strace is installed, no file syscall names the excluded folder at all
if command -v strace >/dev/null; then
  mkdir -p "$t/c/MathWorks"; m=$(fresh b1r); mv "$m" "$t/c/MathWorks/r"
  strace -f -qq -e trace=%file,%desc -o "$t/trace" python3 -B "$dc" --now "$NOW" --all --scan "$t/c/MathWorks" --scan "$t/c/MathWorks/r" --exclude MathWorks >/dev/null 2>&1; rc=$?
  hits=$(grep -v 'execve(' "$t/trace" | grep -F "$t/c/MathWorks")
  [ $rc = 0 ] && [ -z "$hits" ] && ok "B1: strace shows no syscall on an excluded --scan folder" || { bad "B1: syscalls on an excluded --scan folder (rc=$rc)"; show "$hits"; }
else
  echo "note B1 strace check skipped (no strace)"
fi
# ---- round 4 (Cassie): a symlink is judged by its readlink string, hop by hop, never by stat'ing its target ----
# link -> .../MathWorks/hop, and hop -> a drifted repo elsewhere: the first hop already names MathWorks, so excluded
t="$work/lx"; mkdir -p "$t/MathWorks" "$t/s" "$t/u"; m=$(fresh lxr); commit_at "$m" 60 m1; commit_at "$m" 50 m2
ln -s "$m" "$t/MathWorks/hop"; ln -s "$t/MathWorks/hop" "$t/s/link"; ln -s "$t/MathWorks" "$t/u/lk"
out=$(run --all --scan "$t/s" --exclude MathWorks 2>&1); rc=$?
[ $rc = 0 ] && grep -q "excluded $t/s/link" <<<"$out" && ok "lexical: a child whose first hop names MathWorks is excluded, whatever the target (rc 0)" || { bad "lexical: child linked through MathWorks (want rc=0 excluded, got rc=$rc)"; show "$out"; }
if command -v strace >/dev/null; then
  strace -f -qq -e trace=%file,%desc -o "$t/trace" python3 -B "$dc" --now "$NOW" --all --scan "$t/s" --scan "$t/u/lk/hop" --exclude MathWorks >/dev/null 2>&1; rc=$?
  hits=$(grep -v 'execve(' "$t/trace" | grep -F "$t/MathWorks")
  [ $rc = 0 ] && [ -z "$hits" ] && ok "lexical: strace shows no syscall on a symlink target inside MathWorks" || { bad "lexical: syscalls inside MathWorks through a symlink (rc=$rc)"; show "$hits"; }
else
  echo "note lexical strace check skipped (no strace)"
fi
echo "drift_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
