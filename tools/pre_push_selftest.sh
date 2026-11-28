#!/usr/bin/env bash
# Selftest for .githooks/pre-push. Two parts, both must pass (merge of #42 and #40; Jeff 2026-11-27 10:42 PM CT):
#   part 1 (wiring, #42): real pushes to a local bare origin with stub scanners; the hook judges the PUSHED commit
#          in a clean --shared clone, never the pusher's working tree.
#   part 2 (real scanners, #40): the repo's own rules_scan.sh and no_suppress_scan.py on fixture commits of a
#          --shared clone of this repo.
# usage: tools/pre_push_selftest.sh               both parts (rc 0 only when both pass)
#        tools/pre_push_selftest.sh --part 1|2    one part
# PRE_PUSH_HOOK points both parts at another hook (the mutation runner uses this). Scratch under $ZF_SCRATCH, not /tmp.
set -u
if [ "${1:-}" != --part ]; then
  bash "$0" --part 1; r1=$?
  bash "$0" --part 2; r2=$?
  echo "pre_push_selftest: $([ $r1 = 0 ] && [ $r2 = 0 ] && echo PASS || echo FAIL) (part 1 rc=$r1, part 2 rc=$r2)"
  [ $r1 = 0 ] && [ $r2 = 0 ]; exit
fi
case "${2:-}" in 1|2) ;; *) echo "usage: tools/pre_push_selftest.sh [--part 1|2]" >&2; exit 2;; esac
if [ "$2" = 1 ]; then
# ---------------------------------------------------------------- part 1 (#42) ----------------------------------------
# Selftest for .githooks/pre-push (plan v2 C-HOOK). The hook must judge the PUSHED commit, never the pusher's working
# tree: every scanner runs inside a clean --shared clone checked out at the pushed sha. Root cause of the 2026-11-25
# 2:11 PM CT refusal of d41d938 (pushed from a detached LINKED worktree): git exports GIT_DIR=<that worktree's
# gitdir> to hooks, so the #40 hook's `git -C "$t" checkout` acted on the pusher's repo and left the scratch clone
# empty ("could not check out ... or it has no tools/no_suppress_scan.py"); and rules_scan ran from the pusher's
# worktree, where its count checks read the working-tree files (rules_scan.py ROOT), so dirty files can fail (or
# pass) a push of a commit that is clean on its own. The hook must unset GIT_DIR/GIT_WORK_TREE/GIT_INDEX_FILE. A deliberate red-first commit (tests that fail, no fix yet) must push.
# Real pushes to a local bare origin with core.hooksPath set; the scanners are stubs that fail on marker files, so
# this tests the hook's wiring, not the scanners (they have their own selftests). Scratch under $ZF_SCRATCH, not /tmp.
# PRE_PUSH_HOOK can point at another hook (the mutation runner, and the red run against #40's hook, use this).
# The fixture identity below is local to this process; never export it in a shell that makes real commits.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
hook="${PRE_PUSH_HOOK:-$here/../.githooks/pre-push}"
root="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/pre-push-selftest"
mkdir -p "$root" || exit 2
work=$(mktemp -d "$root/run.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
export GIT_AUTHOR_NAME=t GIT_AUTHOR_EMAIL=t@t GIT_COMMITTER_NAME=t GIT_COMMITTER_EMAIL=t@t
export GIT_CONFIG_NOSYSTEM=1 GIT_CONFIG_GLOBAL=/dev/null
export ZF_SCRATCH="$work/scratch"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "ok   $1"; }
bad() { fail=$((fail+1)); echo "FAIL $1"; }

# fixture: origin bare + clone with the hook and two stub scanners
fixture() {
  local d="$work/$1"
  git init -q --bare "$d/origin.git"; git init -q -b w "$d/r"; local r="$d/r"
  git -C "$r" remote add origin "$d/origin.git"
  mkdir -p "$r/.githooks" "$r/tools"
  cp "$hook" "$r/.githooks/pre-push"; chmod +x "$r/.githooks/pre-push"
  cat > "$r/tools/rules_scan.sh" <<'S'
#!/usr/bin/env bash
# stub: record where and how it ran; fail when the scanned tree has a VIOLATION file
echo "rules_scan cwd=$(pwd) args=$*" >> "$STUBLOG"
[ -e VIOLATION ] && { echo "stub rules_scan: VIOLATION"; exit 1; }
exit 0
S
  cat > "$r/tools/no_suppress_scan.py" <<'S'
import os, sys
root = sys.argv[sys.argv.index("--root") + 1] if "--root" in sys.argv else "."
open(os.environ["STUBLOG"], "a").write("no_suppress root=%s\n" % os.path.abspath(root))
sys.exit(1 if os.path.exists(os.path.join(root, "SUPPRESS")) else 0)
S
  chmod +x "$r/tools/rules_scan.sh"
  echo base > "$r/f.txt"; git -C "$r" add -A; git -C "$r" commit -q -m base
  git -C "$r" push -q origin w 2>/dev/null          # first push before the hook is enabled
  git -C "$r" config core.hooksPath .githooks
  echo "$r"
}
# fx NAME : build the fixture in this shell, so r and STUBLOG are set here (not in a $(...) subshell)
fx() { r=$(fixture "$1"); export STUBLOG="$work/$1/stub.log"; : > "$STUBLOG"; }
# push NAME WANT(0 ok / 1 refused) REPO [refspec]
push() {
  local name=$1 want=$2 r=$3 spec=${4:-w}
  local out rc
  out=$(cd "$r" && RATCHET_BASE="$(git rev-list --max-parents=0 HEAD | tail -1)" git push origin $spec 2>&1); rc=$?
  [ $rc -ne 0 ] && rc=1
  if [ "$rc" = "$want" ]; then ok "$name"; else bad "$name (want $( [ $want = 0 ] && echo pushed || echo refused ), got $( [ $rc = 0 ] && echo pushed || echo refused ))"; printf '%s\n' "$out" | sed 's/^/     | /'; fi
}
c() { local r=$1; shift; git -C "$r" "$@"; }

fx h1; echo a >> "$r/f.txt"; c "$r" commit -qam a
push "clean commit, clean tree: pushed" 0 "$r"
grep -q "rules_scan cwd=$r " "$STUBLOG" && bad "rules_scan ran in the pusher's worktree" || ok "rules_scan did not run in the pusher's worktree"
grep -q "rules_scan cwd=$ZF_SCRATCH/" "$STUBLOG" && ok "rules_scan ran inside the scratch clone" || { bad "rules_scan ran inside the scratch clone"; sed 's/^/     | /' "$STUBLOG"; }
grep -q "no_suppress root=$ZF_SCRATCH/" "$STUBLOG" && ok "no_suppress ran on the scratch clone" || bad "no_suppress ran on the scratch clone"

fx h2; touch "$r/VIOLATION"; c "$r" add VIOLATION; c "$r" commit -qm v
push "commit carries a rules violation: refused" 1 "$r"

fx h3; echo a >> "$r/f.txt"; c "$r" commit -qam a; touch "$r/VIOLATION"
push "clean commit, violation only in the dirty tree (untracked): pushed  [2:11 PM root cause]" 0 "$r"

fx h3b; echo a >> "$r/f.txt"; c "$r" commit -qam a; touch "$r/SUPPRESS"
push "clean commit, suppression only in the dirty tree: pushed" 0 "$r"

fx h3c; touch "$r/VIOLATION"; c "$r" add VIOLATION; c "$r" commit -qm v; rm "$r/VIOLATION"
push "violating commit, violation deleted in the dirty tree only: refused" 1 "$r"

fx h4; touch "$r/SUPPRESS"; c "$r" add SUPPRESS; c "$r" commit -qm s
push "commit carries a suppression: refused" 1 "$r"

fx h5; mkdir -p "$r/tests"; echo 'assert!(false) // red first' > "$r/tests/red.rs"; c "$r" add tests; c "$r" commit -qm "red: test first"
push "red-first commit (failing test, no fix yet) alone: pushed" 0 "$r"
echo 'assert!(true)' > "$r/tests/red.rs"; c "$r" commit -qam "fix"
push "red-first then fix, range pushed: pushed" 0 "$r"

fx h6; echo a >> "$r/f.txt"; c "$r" commit -qam a; old=$(c "$r" rev-parse HEAD)
touch "$r/VIOLATION"; c "$r" add VIOLATION; c "$r" commit -qm v
push "push an older clean sha while HEAD violates: pushed" 0 "$r" "$old:refs/heads/w"
push "then the violating HEAD: refused" 1 "$r" w

fx h7; c "$r" push -q origin w:refs/heads/tmpb >/dev/null 2>&1
push "branch deletion: pushed" 0 "$r" ":refs/heads/tmpb"

fx h8; c "$r" rm -q tools/no_suppress_scan.py; c "$r" commit -qm "no no_suppress"
push "tree without no_suppress_scan.py: refused (fail closed once #40 lands)" 1 "$r"
fx h9; c "$r" rm -q tools/rules_scan.sh; c "$r" commit -qm "no rules_scan"
push "tree without tools/rules_scan.sh: refused (fail closed)" 1 "$r"

fx h10; echo a >> "$r/f.txt"; c "$r" commit -qam a; c "$r" switch -q --detach
push "detached HEAD push of HEAD:w: pushed" 0 "$r" "HEAD:refs/heads/w"

# linked worktree (Paul's setup: worktree of a shared repo, detached): the clone must still check out the sha
fx h12; echo a >> "$r/f.txt"; c "$r" commit -qam a; c "$r" branch w2
push "same sha to two refs: pushed" 0 "$r" "w w2"
fx h13; echo a >> "$r/f.txt"; c "$r" commit -qam a; mkdir -p "$work/h13-bin"
printf '#!/usr/bin/env bash\n[ "${1:-}" = clone ] && { echo "shim: clone refused" >&2; exit 1; }\nexec %s "$@"\n' "$(command -v git)" > "$work/h13-bin/git"; chmod +x "$work/h13-bin/git"
# git puts its exec dir first on a hook's PATH, so drive the hook directly with the shim first
sha=$(c "$r" rev-parse HEAD); out=$(cd "$r" && printf 'refs/heads/w %s refs/heads/w %s\n' "$sha" 0000000000000000000000000000000000000000 \
  | PATH="$work/h13-bin:$PATH" .githooks/pre-push origin "$work/h13/origin.git" 2>&1); rc=$?
[ $rc -ne 0 ] && echo "$out" | grep -q 'FAIL could not check out' && ok "clean clone cannot be made: refused (fail closed)" \
  || { bad "clean clone cannot be made: refused (fail closed) (rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; }
fx h11; c "$r" worktree add -q --detach "$work/h11-wt" HEAD 2>/dev/null
echo a >> "$work/h11-wt/f.txt"; c "$work/h11-wt" commit -qam a
push "push from a detached linked worktree: pushed" 0 "$work/h11-wt" "HEAD:refs/heads/w"

# edges (Cassie N6 X13, X14; N7)
fx h14; echo a >> "$r/f.txt"; c "$r" commit -qam a; sha=$(c "$r" rev-parse HEAD)
out=$(cd "$r" && printf '\nrefs/heads/w %s refs/heads/w %s\n' "$sha" 0000000000000000000000000000000000000000 | .githooks/pre-push origin "$work/h14/origin.git" 2>&1); rc=$?
[ $rc = 0 ] && ok "an empty stdin line is ignored" || { bad "an empty stdin line is ignored (rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; }
grep -q "args=c74e834^\.\.$sha" "$STUBLOG" && ok "default RATCHET_BASE is c74e834^ (CI's range)" || { bad "default RATCHET_BASE is c74e834^"; sed 's/^/     | /' "$STUBLOG"; }
fx h15; echo a >> "$r/f.txt"; c "$r" commit -qam a; c "$r" push -q origin w >/dev/null 2>&1; c "$r" reset -q --hard HEAD~1; echo b >> "$r/f.txt"; c "$r" commit -qam b
push "force push (non-ff) of a clean commit: pushed" 0 "$r" "+w"
touch "$r/VIOLATION"; c "$r" add VIOLATION; c "$r" commit -q --amend -m v
push "force push of a violating commit: refused" 1 "$r" "+w"
fx h16; echo a >> "$r/f.txt"; c "$r" commit -qam a; c "$r" tag -a -m t1 t1; : > "$STUBLOG"
push "annotated tag push: pushed" 0 "$r" "refs/tags/t1"
grep -q "rules_scan cwd=$ZF_SCRATCH/" "$STUBLOG" && ok "annotated tag push: the tagged commit was scanned" || { bad "annotated tag push was not scanned"; sed 's/^/     | /' "$STUBLOG"; }
fx h17; touch "$r/VIOLATION"; c "$r" add VIOLATION; c "$r" commit -qm v; c "$r" tag -a -m t2 t2
push "annotated tag on a violating commit: refused" 1 "$r" "refs/tags/t2"

# scratch: nothing left behind, nothing in /tmp
left=$(find "$ZF_SCRATCH" -mindepth 2 -maxdepth 2 2>/dev/null | wc -l)
[ "$left" = 0 ] && ok "scratch clones removed after every push" || { bad "scratch left $left entries"; find "$ZF_SCRATCH" -maxdepth 2 | sed 's/^/     | /'; }
if sed 's/#.*//' "$hook" | grep -q '/tmp'; then bad "hook uses /tmp"; else ok "hook does not use /tmp"; fi
echo "pre_push_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
exit $?
fi
# ---------------------------------------------------------------- part 2 (#40) ----------------------------------------
# Selftest for .githooks/pre-push (Jeff 2026-11-24 1:17 PM CT: "pre-push hook"). The hook runs, for every pushed
# commit, the CI rules-scan command (tools/rules_scan.sh "${RATCHET_BASE:-c74e834^}..<sha>", .github/workflows/
# rules-scan.yml) and tools/no_suppress_scan.py on that commit's tree. Fixture: a --shared clone of this repo's HEAD
# with one base commit that drops any clippy.toml threshold keys (empty since clippy.toml is at defaults, 2e50ed1);
# RATCHET_BASE = that base, so each case scans only its own fixture commit.
# usage: tools/pre_push_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
hook="${PRE_PUSH_HOOK:-$root/.githooks/pre-push}"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
scratch="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/pre-push-selftest"  # job scratch, not /tmp (Spence m6); like fold_pr1_selftest
mkdir -p "$scratch" || exit 2
w=$(mktemp -d "$scratch/run.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
if [ ! -x "$hook" ]; then
  echo "FAIL: $hook missing or not executable"; echo "pre_push_selftest: FAIL"; exit 1
fi
c="$w/clone"
git clone -q --shared "$root" "$c" || exit 2
g() { git -C "$c" -c user.name=t -c user.email=t@t "$@"; }
sed -i -E '/^[a-z-]*(threshold|max-[a-z-]*|-limit) *=/d' "$c/clippy.toml"
printf 'pre_push_selftest fixture base\n' > "$c/.pre-push-fixture-base"  # a base commit even when clippy.toml had no keys
g add -A; g commit -qm "fixture base: no clippy thresholds"; base=$(g rev-parse HEAD)
mkc() { # mkc FILE CONTENT -> sha of a commit on base adding FILE
  g checkout -q --detach "$base"; mkdir -p "$c/$(dirname "$1")"; printf '%b' "$2" > "$c/$1"
  g add -A; g commit -qm "case $1"; g rev-parse HEAD
}
zero=0000000000000000000000000000000000000000
# case WANT LABEL SHA [OUTPUT-REGEX]
case_() {
  local want=$1 label=$2 sha=$3 rx=${4:-}; n=$((n + 1))
  local out rc
  out=$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$sha" "$zero" \
        | RATCHET_BASE="$base" bash "$hook" origin "$w/remote" 2>&1); rc=$?
  if [ "$want" = pass ] && [ $rc -ne 0 ]; then bad "$label: expected pass, rc=$rc"; echo "$out" | tail -5
  elif [ "$want" = fail ] && [ $rc -eq 0 ]; then bad "$label: expected fail, rc=0"
  elif [ -n "$rx" ] && ! grep -Eq "$rx" <<<"$out"; then bad "$label: output lacks /$rx/"; echo "$out" | tail -5
  else echo "ok $n $label (rc=$rc)"; fi
}
clean=$(mkc zenith-float-num/src/zz_fixture.rs '/// fixture\npub fn zz(x: u32) -> u32 {\n    x + 1\n}\n')
clone=$(mkc zenith-float-num/src/zz_fixture.rs 'pub fn zz(x: &Vec<u8>) -> Vec<u8> {\n    x.clone()\n}\n')
tallow=$(mkc tests/zz_fixture.rs '#[allow(dead_code)]\nfn zz() {}\n')
lallow=$(mkc zenith-float-num/src/zz_fixture.rs '#[allow(clippy::too_many_arguments)]\npub fn zz() {}\n')
relaxed=$(g checkout -q --detach "$base" && printf 'too-many-arguments-threshold = 10\n' >> "$c/clippy.toml" \
  && g commit -qam "case relaxed clippy.toml" && g rev-parse HEAD)
unsafe_=$(mkc zenith-float-num/src/zz_fixture.rs 'pub fn zz(p: *const u8) -> u8 {\n    unsafe { *p }\n}\n')
# Spence m3 on 0a1cb6c: a tracked file that matches .gitignore is still pushed (and CI scans it); the hook must see it
ignored=$(g checkout -q --detach "$base" && printf 'zz_ignored.rs\n' >> "$c/.gitignore" \
  && printf '#[allow(dead_code)]\nfn zz() {}\n' > "$c/tests/zz_ignored.rs" && g add -A && g add -f tests/zz_ignored.rs \
  && g commit -qm "case tracked file matching .gitignore" && g rev-parse HEAD)
g checkout -q --detach "$tallow"   # the working tree is NOT the pushed commit: the hook must scan the pushed tree
case_ pass "clean commit passes both scans" "$clean" 'pre-push: ok'
case_ fail "added .clone() in library code fails rules-scan" "$clone" 'FAIL rules-scan'
case_ fail "#[allow] in tests/ passes rules-scan but fails no-suppress" "$tallow" 'ok rules-scan.*|FAIL no-suppress'
case_ fail "#[allow] in library code fails" "$lallow" 'FAIL (rules-scan|no-suppress)'
case_ fail "unsafe in library code fails no-suppress" "$unsafe_" 'FAIL no-suppress'
case_ fail "a relaxed clippy.toml threshold fails no-suppress" "$relaxed" '\[clippy-config\] too-many-arguments-threshold'
case_ fail "#[allow] in a tracked file that matches .gitignore fails no-suppress" "$ignored" 'tests/zz_ignored.rs:1 \[allow\]'
n=$((n + 1)); out=$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$zero" "$clean" | RATCHET_BASE="$base" bash "$hook" origin "$w/remote" 2>&1)
if [ $? -eq 0 ]; then echo "ok $n branch delete (nothing to scan) passes"; else bad "branch delete: $out"; fi
grep -q 'FAIL no-suppress' <<<"$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$tallow" "$zero" | RATCHET_BASE="$base" bash "$hook" origin x 2>&1)" \
  && grep -q 'ok rules-scan' <<<"$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$tallow" "$zero" | RATCHET_BASE="$base" bash "$hook" origin x 2>&1)" \
  && echo "ok $((n + 1)) tests/ #[allow]: rules-scan ok, no-suppress FAIL (why the hook runs both)" || bad "tests/ allow split"
[ $fail -eq 0 ] && echo "pre_push_selftest: ok ($n cases + split)" || echo "pre_push_selftest: FAIL"
exit $fail
