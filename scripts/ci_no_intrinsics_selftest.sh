#!/usr/bin/env bash
# Self-test for scripts/ci_no_intrinsics.sh's unsafe scan. `unsafe` inside a comment (//, ///, //!, /* */, nested and multi-line
# block comments) must not flag. Real `unsafe {`, `unsafe fn` and `unsafe impl` in code must flag, including code that follows a
# block comment or a string containing "//" on the same line. Files on the allowlist stay permitted.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
t="$(mktemp -d "${HOME}/.cache/zf-ci-no-intrinsics-selftest.XXXXXX")" || { echo "selftest: mktemp failed"; exit 2; }
trap 'rm -rf "$t"' EXIT
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
n=0
mkrepo() { # mkrepo <file content>  -> fresh git repo at $t/r with scripts/ and src/x.rs
  rm -rf "$t/r"; mkdir -p "$t/r/scripts" "$t/r/src" "$t/r/allowed"
  cp "$here/ci_no_intrinsics.sh" "$t/r/scripts/"; [ -f "$here/rust_code_only.py" ] && cp "$here/rust_code_only.py" "$t/r/scripts/"
  printf 'allowed/\n' > "$t/r/scripts/no_intrinsics_allowlist.txt"
  printf '%b' "$1" > "$t/r/src/x.rs"; printf 'unsafe fn ok() {}\n' > "$t/r/allowed/a.rs"
  git -C "$t/r" init -q && git -C "$t/r" add -A
}
pass_case() { n=$((n+1)); mkrepo "$2"; out=$(bash "$t/r/scripts/ci_no_intrinsics.sh" 2>&1) || fail "$1: flagged but must pass: $out"; }
flag_case() { n=$((n+1)); mkrepo "$2"; out=$(bash "$t/r/scripts/ci_no_intrinsics.sh" 2>&1) && fail "$1: not flagged: $out"
  echo "$out" | grep -q 'unsafe outside the allowlist: src/x.rs' || fail "$1: wrong report: $out"; }
pass_case "clean"                 'fn f() {}\n'
pass_case "inner doc comment"     '//! Environment (read through `std::env::var`, no unsafe; ...)\nfn f() {}\n'
pass_case "outer doc comment"     '/// unsafe is never used here\nfn f() {}\n'
pass_case "line comment after code" 'fn f() {} // no unsafe\n'
pass_case "block comment"         '/* unsafe */ fn f() {}\n'
pass_case "multi-line block"      '/*\n * unsafe { }\n */\nfn f() {}\n'
pass_case "nested block"          '/* a /* inner */ unsafe */ fn f() {}\n'
pass_case "allowlisted file"      'fn f() {}\n'
flag_case "unsafe block"          'fn f() { unsafe { g() } }\n'
flag_case "unsafe fn"             'unsafe fn f() {}\n'
flag_case "unsafe impl"           'unsafe impl Send for T {}\n'
flag_case "code after block comment" '/* note */ unsafe { g() }\n'
flag_case "code after multi-line block" '/* a\n b */ unsafe fn f() {}\n'
flag_case "string with // then unsafe" 'fn f() { let s = "//"; unsafe { g() } }\n'
echo "ci_no_intrinsics selftest: ok ($n cases)"
