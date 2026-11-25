#!/usr/bin/env bash
# Selftest for tools/no_held_scan.sh (Jeff 2026-11-25 3:36 AM CT brief, item 3): each case is a fixture git repo with one
# file. A held marker or counter, a None-gap mode list, a known-wrong list or "today it returns" wording, an ignored or
# xfail test and a mode `continue` must FAIL; a HOLD file, past-tense history, `mod held_*;` and a mode branch that is
# not a skip must PASS. Then the real tree must be clean.
# usage: tools/no_held_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
scan="$root/tools/no_held_scan.sh"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/no-held-selftest"; mkdir -p "$zf_tmp" || exit 2  # job scratch
w=$(mktemp -d "$zf_tmp/run.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
case_() {
  local want=$1 file=$2 kind=$3 body=$4
  n=$((n + 1)); local r="$w/r$n"; mkdir -p "$r/$(dirname "$file")"
  printf '%s\n' "$body" > "$r/$file"
  git -C "$r" init -q && git -C "$r" add -A && git -C "$r" -c user.name=t -c user.email=t@t commit -qm f
  "$scan" --root "$r" > "$w/out$n" 2>&1; local rc=$?
  if [[ $want == pass ]]; then
    [[ $rc -eq 0 ]] || bad "case $n want PASS, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-200)"
  else
    [[ $rc -eq 1 ]] && grep -q "^$file:[0-9]* \[$kind\]" "$w/out$n" || bad "case $n want FAIL [$kind], rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-200)"
  fi
}
case_ fail src/lib.rs held '// M27-F12 (HELD): pow(1, NaN) is ruled 1'
case_ fail src/lib.rs held 'fn f() { let mut held = 0; }'
case_ fail src/lib.rs held 'fn f(h: &mut usize) { held_cells += 1; }'
case_ fail src/lib.rs held 'fn f(x: u64) { dd_held(x); }'
case_ fail src/lib.rs held 'const F: &[usize] = &include!("grid_held/f25.in");'
case_ fail tests/t.rs none-gap 'const NONE_GAP_MODES: [RoundingMode; 6] = [RoundingMode::Up; 6];'
case_ fail tests/t.rs none-gap '    // NONE-GAP F-B: None overflow gives +-Inf'
case_ fail tests/t.rs none-gap '/// KNOWN GAP (ZF-M27 F-A): None is not swept'
case_ fail tools/x.sh none-gap 'KNOWN_GAPS="F-B F-C"'
case_ fail tests/t.rs inverted 'const KNOWN_WRONG: &[&str] = &["pow(1, NaN)"];'
case_ fail tests/t.rs inverted '/// zenith-float today returns NaN for pow(1, NaN)'
case_ fail tests/t.rs inverted '/// pow(1, NaN) returns NaN today'
case_ fail tests/t.py inverted 'KNOWN_WRONG_CELLS = {"a"}  # known-wrong'
case_ fail tests/t.rs xfail '#[ignore = "red: not implemented"]
#[test]
fn t() {}'
case_ fail tests/t.py xfail '@pytest.mark.xfail
def test_a(): pass'
case_ fail tests/t.rs mode-skip 'fn f() {
    for rm in MODES {
        if rm == RoundingMode::None {
            continue;
        }
        g(rm);
    }
}'
case_ fail tests/t.rs mode-skip 'fn f() { for rm in MODES { if matches!(rm, RoundingMode::None | RoundingMode::Work) { continue; } g(rm); } }'
case_ pass src/held.rs - '// HOLD H99: held, not compiled (no `mod` line names this file). Index and use-or-delete rule: HOLD.md.
// // HELD: let mut held = 0; held_cells += 1;
// #[ignore = "red"]'
case_ pass src/lib.rs - '// Re-pinned (ZF-M27 F-B): was 1_196_293, with the None overflow cells in (NONE-GAP F-B).
mod held_mul_div_reds;
// HOLD H12: mod held_logdet_reds;
/// The bound held at every row; the function returns 1.
fn f(held_out: usize) -> bool { held_out == 3 }'
case_ pass src/lib.rs - 'fn f(rm: RoundingMode) -> u8 { if rm == RoundingMode::None { return 0; } 1 }'
case_ pass src/lib.rs - '/// It was NaN before the fix (F-12); it is 1 now, as IEEE 754 rules.
fn ignore_case() {}'
"$scan" --root "$root" > "$w/tree" 2>&1 || bad "real tree: $(tr '\n' '|' < "$w/tree" | cut -c1-300)"
[[ $fail -eq 0 ]] && echo "no_held_selftest: ok ($n cases + tree)" || exit 1
