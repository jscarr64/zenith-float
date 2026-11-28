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
# every run is capped at 60 s (normal runs take about 1 s), so a run that hangs FAILs that case (rc 124)
# instead of hanging the selftest
run() { timeout 60 python3 -B "$dc" --now "$NOW" "$@"; }
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
# The hook also runs #40's library rules (tools/rules_scan.sh, cargo fmt, cargo clippy; merge of #42 and #40). They are
# stubs here (exit 0), so this part tests the unpushed cap only; the rules part is tools/test_rules_enforcement.sh part 4.
stubs="$work/hook-stubs"; mkdir -p "$stubs"; printf '#!/usr/bin/env bash\nexit 0\n' > "$stubs/cargo"; chmod +x "$stubs/cargo"
hookfx() { mkdir -p "$1/.githooks" "$1/tools"; cp "$hook" "$1/.githooks/pre-commit"; chmod +x "$1/.githooks/pre-commit"; cp "$dc" "$1/tools/drift_check.py"
  printf '#!/usr/bin/env bash\nexit 0\n' > "$1/tools/rules_scan.sh"; chmod +x "$1/tools/rules_scan.sh"; }
opath=$PATH; PATH="$stubs:$PATH"
r=$(fresh hk); hookfx "$r"
git -C "$r" config core.hooksPath .githooks
git -C "$r" add .githooks tools; git -C "$r" commit -q -m hooks 2>/dev/null && ok "hook: 1st unpushed commit allowed" || bad "hook: 1st unpushed commit allowed"
echo z >> "$r/f.txt"; git -C "$r" add f.txt
if git -C "$r" commit -q -m second 2>/dev/null; then bad "hook: 2nd unpushed commit refused"; else ok "hook: 2nd unpushed commit refused"; fi
git -C "$r" push -q 2>/dev/null
git -C "$r" commit -q -m second 2>/dev/null && ok "hook: commit allowed again after push" || bad "hook: commit allowed again after push"
e="$work/hku"; git init -q -b w "$e"; hookfx "$e"
git -C "$e" config core.hooksPath .githooks; git -C "$e" add -A
git -C "$e" commit -q -m first 2>/dev/null && ok "hook: first commit of a new repo allowed" || bad "hook: first commit of a new repo allowed"
# outside a work tree the hook fails closed with its own message (Cassie N6 X16)
out=$(cd "$work/notgit" && bash "$hook" 2>&1); rc=$?
[ $rc = 2 ] && grep -q 'pre-commit: not in a git work tree' <<<"$out" && ok "hook outside a work tree: rc 2" || { bad "hook outside a work tree (rc=$rc)"; show "$out"; }
PATH=$opath
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
[ $rc = 2 ] && ok "active: unreadable activity fails closed (checked, rc=$rc)" || { bad "active: unreadable activity fails closed (want rc=2, got rc=$rc)"; show "$out"; }
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
# ---- allowlist (Jeff 2026-11-27 4:42 PM CT): --only-origin OWNER/NAME makes --scan check only work trees whose
# remote.origin.url is GitHub OWNER/NAME (https, scp or ssh form, any case, with or without .git). Other origins,
# local-path origins and no origin are skipped and named. The roots are agent folders: each --scan DIR and its
# direct children, so a clone two levels down is not seen. ----
# zf PATH URL [clean|dirty|ahead] : a work tree at PATH with one pushed commit, then its origin URL set to URL
# ("-" keeps the local bare origin, "none" removes the remote)
zf() {
  local p=$1 o; o="$work/zforigins/$(echo "$1" | tr / _).git"
  mkdir -p "$(dirname "$p")" "$work/zforigins"; git init -q --bare "$o"; git init -q -b w "$p"
  git -C "$p" remote add origin "$o"; commit_at "$p" 7200 base; git -C "$p" push -q -u origin w 2>/dev/null
  case "$2" in -) ;; none) git -C "$p" remote remove origin;; *) git -C "$p" remote set-url origin "$2";; esac
  case "${3:-clean}" in
    dirty) echo d > "$p/old.txt"; touch -d "@$((NOW-7200))" "$p/old.txt";;
    ahead) commit_at "$p" 60 a1; commit_at "$p" 60 a2;;
  esac
}
T="$work/team"; ZF=https://github.com/jscarr64/zenith-float.git; oo=(--only-origin jscarr64/zenith-float)
# allow NAME WANT_RC [re:REGEX]... [nore:REGEX]... -- ARGS... : run drift_check --only-origin jscarr64/zenith-float ARGS
allow() {
  local name=$1 want=$2; shift 2; local res=() nres=()
  while [ "$1" != -- ]; do case "$1" in re:*) res+=("${1#re:}");; nore:*) nres+=("${1#nore:}");; esac; shift; done; shift
  local out rc good=1; out=$(run "${oo[@]}" "$@" 2>&1); rc=$?
  [ "$rc" = "$want" ] || good=0
  local x; for x in "${res[@]+"${res[@]}"}"; do grep -Eq -- "$x" <<<"$out" || good=0; done
  for x in "${nres[@]+"${nres[@]}"}"; do grep -Eq -- "$x" <<<"$out" && good=0; done
  if [ $good = 1 ]; then ok "allowlist: $name (rc=$rc)"; else bad "allowlist: $name (want rc=$want, got rc=$rc)"; show "$out"; fi
}
# a clean zenith-float clone in its own agent folder: a skip case scans it too, so the allowlist matches something
zf "$T/anchor-work/zf" "$ZF"; A=(--scan "$T/anchor-work")
allow "zero roots: usage error" 2 re:'give --repo or --scan' --
mkdir -p "$T/empty-work"
allow "an empty agent folder: nothing checked" 0 re:'\(0 repo\(s\) checked' -- --scan "$T/empty-work"
zf "$T/a-work/acc" https://github.com/jscarr64/accumath.git dirty
allow "a non-zenith-float origin is skipped, even with old dirty work" 0 re:"skip $T/a-work/acc \(origin " nore:FAIL -- --scan "$T/a-work" "${A[@]}"
zf "$T/b-work/zb" https://github.com/jscarr64/zenith-bound.git ahead
allow "a zenith-bound origin is skipped, even 2 ahead" 0 re:"skip $T/b-work/zb \(origin " nore:FAIL -- --scan "$T/b-work" "${A[@]}"
zf "$T/c-work/loc" - dirty
allow "a local-path origin is skipped" 0 re:"skip $T/c-work/loc \(origin " nore:FAIL -- --scan "$T/c-work" "${A[@]}"
zf "$T/d-work/sub/deep" "$ZF" dirty
allow "a zenith-float clone two levels down is not seen" 0 re:'\(0 repo\(s\) checked' nore:deep -- --scan "$T/d-work"
zf "$T/e-real/zf" "$ZF" dirty; ln -s "$T/e-real" "$T/e-work"
allow "a symlinked agent folder is followed: its dirty clone FAILs" 1 re:"FAIL $T/e-work/zf: dirty work" -- --scan "$T/e-work"
zf "$T/f-work/zf" "$ZF" dirty
allow "a dirty zenith-float clone FAILs" 1 re:"FAIL $T/f-work/zf: dirty work" -- --scan "$T/f-work"
zf "$T/g-work/zf" "$ZF" ahead
allow "a zenith-float clone 2 ahead FAILs" 1 re:"FAIL $T/g-work/zf w: 2 commits ahead" -- --scan "$T/g-work"
zf "$T/h-work/zf" "$ZF"
allow "a clean zenith-float clone is checked: ok" 0 re:"ok   $T/h-work/zf\$" re:'\(1 repo\(s\) checked' -- --scan "$T/h-work"
# every GitHub URL form of jscarr64/zenith-float counts; near names, other owners, other hosts and no origin do not
i=0
for u in https://github.com/jscarr64/zenith-float https://github.com/jscarr64/zenith-float.git \
         git@github.com:jscarr64/zenith-float.git ssh://git@github.com/jscarr64/zenith-float \
         https://GitHub.com/JScarr64/Zenith-Float.git; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "origin $u is zenith-float: checked" 1 re:"FAIL $T/u$i-work/zf w: 2 commits ahead" -- --scan "$T/u$i-work"
done
for u in https://github.com/jscarr64/zenith-float-old.git https://github.com/other/zenith-float.git \
         https://gitlab.com/jscarr64/zenith-float.git https://github.com/jscarr64/zenith-float/x.git none; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "origin $u is not zenith-float: skipped" 0 re:"skip $T/u$i-work/zf \(" nore:FAIL -- --scan "$T/u$i-work" "${A[@]}"
done
# the agent folder itself is filtered the same way, and --repo is never filtered (the pre-commit hook uses it)
zf "$T/acc-as-root" https://github.com/jscarr64/accumath.git ahead
allow "a non-zenith-float --scan DIR itself is skipped" 0 re:"skip $T/acc-as-root \(" nore:FAIL -- --scan "$T/acc-as-root" "${A[@]}"
allow "--repo is checked whatever its origin" 1 re:"FAIL $T/acc-as-root w: 2 commits ahead" -- --repo "$T/acc-as-root"
# ---- round 7 B1 (Cassie): --only-origin is trimmed and lower-cased, must be OWNER/NAME (else rc 2), and matching
# none of the work trees under --scan is an ERROR (rc 2), not a silent ok ----
for o in JScarr64/Zenith-Float " jscarr64/zenith-float "; do
  out=$(run --only-origin "$o" --scan "$T/g-work" 2>&1); rc=$?
  [ $rc = 1 ] && grep -q "FAIL $T/g-work/zf w: 2 commits ahead" <<<"$out" && ok "r7 B1: --only-origin '$o' is normalized: checked (rc=1)" || { bad "r7 B1: --only-origin '$o' (want rc=1 FAIL, got rc=$rc)"; show "$out"; }
done
for o in https://github.com/jscarr64/zenith-float jscarr64/zenith-float.git jscarr64/ jscarr64 /zenith-float a/b/c "" "jscarr64/zenith float"; do
  out=$(run --only-origin "$o" --scan "$T/g-work" 2>&1); rc=$?
  [ $rc = 2 ] && grep -q -- "--only-origin must be OWNER/NAME" <<<"$out" && ok "r7 B1: malformed --only-origin '$o': rc 2" || { bad "r7 B1: malformed --only-origin '$o' (want rc=2, got rc=$rc)"; show "$out"; }
done
allow "r7 B1: work trees under --scan but none matches: ERROR rc 2" 2 re:"ERROR --only-origin jscarr64/zenith-float matched none of 1 work tree" -- --scan "$T/a-work"
allow "r7 B1: an empty agent folder beside a zf one is fine" 0 -- --scan "$T/empty-work" --scan "$T/h-work"
# ---- round 7 B2 (Cassie): every real URL form of our repo is CHECKED (origin read with git ls-remote --get-url, so
# insteadOf counts); credentials in the URL are never printed ----
for u in ssh://git@github.com:22/jscarr64/zenith-float.git git@github.com:jscarr64/zenith-float.git/ \
         https://github.com/jscarr64/zenith-float/ https://www.github.com/jscarr64/zenith-float \
         https://user:s3cret@github.com/jscarr64/zenith-float.git https://x-access-token:ghs_s3cret@github.com/jscarr64/zenith-float \
         http://github.com/jscarr64/zenith-float git://github.com/jscarr64/zenith-float.git ssh://github.com/jscarr64/zenith-float \
         github.com:jscarr64/zenith-float HTTPS://GITHUB.COM/JSCARR64/ZENITH-FLOAT.GIT; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "r7 B2: origin $u is zenith-float: checked" 1 re:"FAIL $T/u$i-work/zf w: 2 commits ahead" nore:s3cret -- --scan "$T/u$i-work"
done
i=$((i+1)); zf "$T/u$i-work/zf" gh:jscarr64/zenith-float ahead; git -C "$T/u$i-work/zf" config url.https://github.com/.insteadOf gh:
allow "r7 B2: an insteadOf alias that resolves to zenith-float: checked" 1 re:"FAIL $T/u$i-work/zf w: 2 commits ahead" -- --scan "$T/u$i-work"
# clearly another repo: skipped (look-alike hosts, github.com mid-URL, extra or empty parts, near names, a bare path)
for u in https://github.com.evil.com/jscarr64/zenith-float https://evilgithub.com/jscarr64/zenith-float \
         https://evil.com/github.com/jscarr64/zenith-float https://github.com/jscarr64/zenith.git-float \
         https://github.com/jscarr64/zenith-floatgit https://github.com/jscarr64/zenith-float.git.git \
         https://github.com//zenith-float https://github.com/jscarr64/ jscarr64/zenith-float file:///x/jscarr64/zenith-float; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "r7 B2: origin $u is another repo: skipped" 0 re:"skip $T/u$i-work/zf \(" nore:FAIL -- --scan "$T/u$i-work" "${A[@]}"
done
# what git itself treats as a local path (no ":" before the first "/") is skipped, even with a ":" later on;
# a near name that only ends in "git" is another repo
for u in /x/y:z/jscarr64/zenith-float ./a:b rel/a:b https://github.com/jscarr64/zenith-floatxgit; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "r7 B2: origin $u is another repo: skipped" 0 re:"skip $T/u$i-work/zf \(" nore:FAIL -- --scan "$T/u$i-work" "${A[@]}"
done
# an origin that cannot be read or parsed is an ERROR (rc 2), never a skip, and its URL is not printed
for u in ftp://u:s3cret@github.com/jscarr64/zenith-float https:///jscarr64/zenith-float u@v@github.com:jscarr64/zenith-float; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "r7 B2: origin $u cannot be parsed: ERROR" 2 re:"ERROR $T/u$i-work/zf: origin " nore:s3cret nore:"skip $T/u$i-work" -- --scan "$T/u$i-work" "${A[@]}"
done
i=$((i+1)); zf "$T/u$i-work/zf" "$ZF" ahead; git -C "$T/u$i-work/zf" config --add remote.origin.url https://github.com/other/x.git
allow "r7 B2: two origin url values: ERROR" 2 re:"ERROR $T/u$i-work/zf: origin " nore:"skip $T/u$i-work" -- --scan "$T/u$i-work" "${A[@]}"
i=$((i+1)); zf "$T/u$i-work/zf" "$ZF" ahead; git -C "$T/u$i-work/zf" config --unset remote.origin.url
allow "r7 B2: an origin remote with no url: ERROR" 2 re:"ERROR $T/u$i-work/zf: origin " nore:"skip $T/u$i-work" -- --scan "$T/u$i-work" "${A[@]}"
# ---- round 7 B3 CR06 (Cassie): the work-tree check runs before the origin check, so a broken .git is the work-tree ERROR ----
mkdir -p "$T/k-work/broken"; echo "gitdir: $work/no-such-gitdir" > "$T/k-work/broken/.git"
allow "r7 CR06: a broken .git under --only-origin: work-tree ERROR, rc 2" 2 re:"ERROR $T/k-work/broken has a .git but is not a readable work tree" nore:"skip $T/k-work/broken" -- --scan "$T/k-work" "${A[@]}"
# ---- round 7 N1 (Cassie): a child that is a symlink loop fails closed (rc 2) with the unreadable-work-tree ERROR;
# a --scan DIR that is a loop is "not a directory" ----
mkdir -p "$T/lp-work"; ln -s l "$T/lp-work/l"
allow "r7 N1: a symlink-loop child: rc 2, ERROR line" 2 re:"ERROR $T/lp-work/l has a .git but is not a readable work tree" -- --scan "$T/lp-work" "${A[@]}"
allow "r7 N1: a --scan DIR that is a symlink loop: rc 2" 2 re:"lp-work/l is not a directory" -- --scan "$T/lp-work/l"
# ---- round 8 (Cassie): tests that kill M03, M04, M06, M07, M08 and M12. No code change; they pin today's
# behavior. M11 (a $ on OWNER_NAME) is equivalent to fullmatch after strip, not a row. ----
# M03/M04: a valid zenith-float URL with trailing junk cannot be parsed (origin_repo turns that into
# ERROR rc 2). A space after the URL is eaten by the path and would be skipped, and git stores a
# newline in the url as two values (already ERROR, on the count, for real code and for the mutant).
# So this calls github_path() in $dc directly: a newline leftover makes fullmatch fail (None);
# match()/search() accept the prefix and would return jscarr64/zenith-float.
got=$(python3 -B -c 'import importlib.util, sys
spec = importlib.util.spec_from_file_location("dc", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
print(m.github_path("https://github.com/jscarr64/zenith-float.git\ntrailing"))
' "$dc")
[ "$got" = "None" ] && ok "r8 M03/M04: a zenith-float URL with trailing junk cannot be parsed" || { bad "r8 M03/M04: trailing junk (want None, got $got)"; }
# M06: a remote whose name contains the substring origin, with no origin remote, is skipped as no origin
# (the companion allowlisted clone keeps the zero-match ERROR from firing)
zf "$T/m06-work/zf" "$ZF" ahead; git -C "$T/m06-work/zf" remote rename origin original
allow "r8 M06: a remote named original with no origin is skipped as no-origin" 0 re:"skip $T/m06-work/zf \(" nore:FAIL -- --scan "$T/m06-work" "${A[@]}"
# M07: host github.com. is not github.com (today: skip)
zf "$T/m07-work/zf" https://github.com./jscarr64/zenith-float.git ahead
allow "r8 M07: origin https://github.com./jscarr64/zenith-float.git is skipped" 0 re:"skip $T/m07-work/zf \(" nore:FAIL -- --scan "$T/m07-work" "${A[@]}"
# M08 / M12: #frag and ?x=1 stay in the path (today: skip). The # case also pins M12 (split on # before parse).
for u in https://github.com/jscarr64/zenith-float.git#frag https://github.com/jscarr64/zenith-float.git?x=1; do
  i=$((i+1)); zf "$T/u$i-work/zf" "$u" ahead
  allow "r8 M08/M12: origin $u is skipped" 0 re:"skip $T/u$i-work/zf \(" nore:FAIL -- --scan "$T/u$i-work" "${A[@]}"
done
echo "drift_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
