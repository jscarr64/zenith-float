#!/usr/bin/env bash
# Selftest for tools/no_debug_assert_scan.sh (Jeff 2026-11-24 4:40 PM CT item 5): each case is a fixture git repo with one
# Rust file. debug_assert!, debug_assert_eq! and debug_assert_ne! in code must FAIL; assert!, comments, the word inside a
# string or doc comment and cfg!(debug_assertions) must PASS. Then the real tree must be clean.
# usage: tools/no_debug_assert_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
scan="$root/tools/no_debug_assert_scan.sh"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/no-debug-assert-selftest"; mkdir -p "$zf_tmp" || exit 2  # job scratch
w=$(mktemp -d "$zf_tmp/run.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
case_() {
  local want=$1 body=$2
  n=$((n + 1)); local r="$w/r$n"; mkdir -p "$r/src"
  printf '%s\n' "$body" > "$r/src/lib.rs"
  git -C "$r" init -q && git -C "$r" add -A && git -C "$r" -c user.name=t -c user.email=t@t commit -qm f
  "$scan" --root "$r" > "$w/out$n" 2>&1; local rc=$?
  if [[ $want == pass ]]; then
    [[ $rc -eq 0 ]] || bad "case $n want PASS, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-200)"
  else
    [[ $rc -eq 1 ]] && grep -q '^src/lib.rs:[0-9]* \[debug-assert\]' "$w/out$n" || bad "case $n want FAIL, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-200)"
  fi
}
case_ fail 'fn f(x: u8) { debug_assert!(x > 0); }'
case_ fail 'fn f(x: u8) { debug_assert_eq!(x, 1); }'
case_ fail 'fn f(x: u8) { debug_assert_ne!(x, 1); }'
case_ fail 'fn f(x: u8) {
    let y = x; debug_assert! (y < 9); // trailing comment
}'
case_ pass 'fn f(x: u8) { assert!(x > 0); assert_eq!(x, 1); }'
case_ pass '// debug_assert!(x) in a comment is fine
/// so is debug_assert! in a doc comment
fn f() -> bool { cfg!(debug_assertions) }'
case_ pass 'fn f(x: u8) { assert!(x > 0); } // was debug_assert!(x > 0)'
case_ pass 'const S: &str = "debug_assert!(x)";'
"$scan" --root "$root" > "$w/tree" 2>&1 || bad "real tree: $(tr '\n' '|' < "$w/tree" | cut -c1-300)"
[[ $fail -eq 0 ]] && echo "no_debug_assert_selftest: ok ($n cases + tree)" || exit 1
