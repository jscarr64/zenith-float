#!/usr/bin/env bash
# Selftest for tools/no_tmp_scan.sh (Spence n3 on #40; Jeff 2026-11-24 4:40 PM CT item 8): each case is a fixture git repo with
# one script. A /tmp work path, a bare mktemp, a TMPDIR mktemp template or a Python tempfile call without dir= must FAIL by
# form; the job-scratch pattern, comments, the shared pop lock files, $dir/tmp paths and tempfile with dir= must PASS.
# usage: tools/no_tmp_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
scan="$root/tools/no_tmp_scan.sh"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/no-tmp-selftest"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp
w=$(mktemp -d "$zf_tmp/run.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
# case WANT FORM FILE CONTENT
case_() {
  local want=$1 form=$2 file=$3 body=$4
  n=$((n + 1)); local r="$w/r$n"; mkdir -p "$r/$(dirname "$file")"
  printf '%s\n' "$body" > "$r/$file"
  git -C "$r" init -q && git -C "$r" add -A && git -C "$r" -c user.name=t -c user.email=t@t commit -qm f
  "$scan" --root "$r" > "$w/out$n" 2>&1; local rc=$?
  if [[ $want == pass ]]; then
    [[ $rc -eq 0 ]] || bad "case $n ($file) want PASS, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-200)"
  else
    [[ $rc -eq 1 ]] && grep -q "^$file:[0-9]* \[$form\]" "$w/out$n" && grep -q "^no_tmp_scan: FAIL [0-9]* finding" "$w/out$n" || bad "case $n ($file) want FAIL [$form], rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-200)"
  fi
}
case_ fail mktemp-bare tools/a.sh 't=$(mktemp -d); echo "$t"'
case_ fail mktemp-bare tools/b.sh 'f="$(mktemp)"'
case_ fail mktemp-bare scripts/c.sh 'if true; then x=$(mktemp -d) || exit 2; fi'
case_ fail tmp-path tools/d.sh 'x=$(mktemp -d /tmp/foo.XXXXXX)'
case_ fail tmp-path tools/e.sh 'y=$(mktemp -d "${TMPDIR:-/tmp}/foo.XXXXXX")'
case_ fail mktemp-tmpdir tools/f.sh 'y=$(mktemp -d "${TMPDIR:-$HOME/.cache}/foo.XXXXXX")'
case_ fail tmp-path tools/g.sh 'cp x /tmp/'
case_ fail py-tempfile tools/h.py 'import tempfile
with tempfile.TemporaryDirectory() as t:
    pass'
case_ fail py-tempfile tools/i.py 'import tempfile
d = tempfile.mkdtemp(prefix="x")'
case_ fail tmp-path tools/j.py 'd = os.environ.get("TMPDIR", "/tmp")'
case_ fail mktemp-bare .githooks/pre-commit 't=$(mktemp -d)'
case_ fail mktemp-bare zenith-float-bound/tools/k.sh 'out=$(mktemp)'
case_ pass - tools/ok1.sh 'zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/job"; mkdir -p "$zf_tmp"
t=$(mktemp -d "$zf_tmp/run.XXXXXX")'
case_ pass - tools/ok2.sh '# a comment naming /tmp and mktemp -d is fine
d=$CARGO_TARGET_DIR/tmp/x  # $dir/tmp is not /tmp'
case_ pass - tools/ok3.sh 'exec 9> /tmp/pop-deep.lock; l=/tmp/pop-gate.lock; p=${P:-/tmp/pop-gate.slot}'
case_ pass - tools/ok4.sh 't=$(mktemp -d "$(zb_scratch)"/x.XXXXXX); real=$(command -v mktemp)'
case_ pass - tools/ok5.py 'import tempfile
with tempfile.TemporaryDirectory(
    dir=scratch()) as t:
    pass'
case_ pass - docs/notes.md 'mktemp -d and /tmp in prose are not scripts'
want=18
[[ $n -eq $want ]] || bad "ran $n cases, want $want"
# the real tree is clean
"$scan" > "$w/tree" 2>&1 || bad "the tree has /tmp work files: $(tr '\n' '|' < "$w/tree" | cut -c1-300)"
if [[ $fail -eq 0 ]]; then echo "no_tmp_selftest: ok ($n cases, tree clean)"; else echo "no_tmp_selftest: FAIL"; fi
exit $fail
