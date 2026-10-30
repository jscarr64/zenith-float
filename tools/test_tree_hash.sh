#!/usr/bin/env bash
# Red test for the gate's tree-clean check (tools/tree_hash.sh + full_gate.sh tree-clean): the check compares
# content hashes (and mode, symlink target) of every tracked and untracked non-ignored file before and after,
# so a change to an ALREADY-dirty file is caught (git status would read ` M` both times). Scratch git repo.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
r="$tmp/repo"; mkdir -p "$r/sub dir"; cd "$r"; git init -q
printf 'a\n' > a.txt; printf 'b\n' > "sub dir/b c.txt"; printf 'x\n' > x.sh; ln -s a.txt link; printf 'target/\n' > .gitignore
git add -A; git -c user.name=s -c user.email=s@invalid commit -qm init
printf 'dirty\n' >> a.txt   # pre-dirty tracked file
fail=0
# case NAME EXPECT(same|diff) MUTATION: hash, mutate, hash again, restore the scratch repo
case_() {
  local n=$1 want=$2 m=$3 pre post got
  pre=$(bash "$root/tools/tree_hash.sh") || { echo "FAIL[$n]: tree_hash.sh failed"; fail=1; return; }
  [[ -n $pre ]] || { echo "FAIL[$n]: empty hash list"; fail=1; }
  eval "$m"
  post=$(bash "$root/tools/tree_hash.sh")
  got=same; [[ "$pre" == "$post" ]] || got=diff
  [[ $got == "$want" ]] || { echo "FAIL[$n]: expected $want, got $got"; fail=1; }
  git checkout -q -- . 2>/dev/null; git clean -qfd; printf 'dirty\n' >> a.txt; rm -rf target; ln -sfn a.txt link; chmod 644 x.sh
  git checkout -q -- a.txt; printf 'dirty\n' >> a.txt
}
case_ untouched same ':'
case_ predirty-rewrite diff 'printf "more\n" >> a.txt'
case_ clean-rewrite diff 'printf "y\n" > x.sh'
case_ delete diff 'rm "sub dir/b c.txt"'
case_ new-untracked diff 'printf "n\n" > new.txt'
case_ mode diff 'chmod 755 x.sh'
case_ symlink diff 'ln -sfn x.sh link'
case_ space-path diff 'printf "z\n" > "sub dir/b c.txt"'
case_ ignored same 'mkdir -p target; printf "o\n" > target/out'
# the gate step: an empty/missing pre-hash file is a FAIL, never a PASS on nothing
grep -q 'tools/tree_hash.sh > ' "$root/tools/full_gate.sh" && grep -q '^run tree-clean .*\[\[ -s .*tools/tree_hash.sh | diff -u' "$root/tools/full_gate.sh" \
  || { echo "FAIL[gate]: full_gate.sh tree-clean does not diff tree_hash.sh against a non-empty pre-hash"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: tree_hash.sh catches content/mode/link/add/delete changes, incl. on pre-dirty files"
exit "$fail"
