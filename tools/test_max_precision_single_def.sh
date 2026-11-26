#!/usr/bin/env bash
# Red test (Paul's order 3, 2026-11-10; P6 64-bit only): there is ONE definition of Precision, MIN_PRECISION and MAX_PRECISION,
# in zenith-float-num/src/defs.rs. MAX_PRECISION is the fixed value 2^60 - 1, with no usize-width cap: it is not derived from
# usize::MAX and uses no `as`. zenith-float-bound/src/precision.rs only re-exports the zf items. The bound mutants for
# them (117-123) are retired: the bound tool cannot mutate a file outside its crate (git apply skips a ../ path), so no row may
# name src/precision.rs, those ids, or a ../ path.
# Each check also runs on a planted copy and must catch it.   usage: tools/test_max_precision_single_def.sh
set -u
root=$(cd "$(dirname "$0")/.." && pwd)
fail=0; bad() { echo "FAIL: $*"; fail=1; }
DEFS=zenith-float-num/src/defs.rs; BP=zenith-float-bound/src/precision.rs; MAN=zenith-float-bound/tools/mutants/manifest.txt
check() {  # DIR -> 0 when single-definition
  local d=$1 r=0 blk
  blk=$(awk '/^pub const MAX_PRECISION/{p=1} p{print} p&&/;[[:space:]]*$/{exit}' "$d/$DEFS")
  [ -n "$blk" ] || { echo "  no MAX_PRECISION in $DEFS"; r=1; }
  grep -qE 'usize::MAX|[[:space:]]as[[:space:]]' <<< "$blk" && { echo "  MAX_PRECISION derived from usize width: $blk"; r=1; }
  grep -qE '^pub const MAX_PRECISION: Precision = \(1 << 60\) - 1;$' <<< "$blk" || { echo "  MAX_PRECISION is not the fixed (1 << 60) - 1: $blk"; r=1; }
  grep -nE '^[[:space:]]*(pub )?(const|type|fn|static)[[:space:]]+(MAX_PRECISION|MIN_PRECISION|Precision|max_precision)\b' "$d/$BP" \
    && { echo "  $BP defines its own precision item"; r=1; }
  grep -qE '^pub use zenith_float_num::\{(Precision, MAX_PRECISION, MIN_PRECISION|MAX_PRECISION, MIN_PRECISION, Precision)\};$' "$d/$BP" \
    || { echo "  $BP does not re-export zenith_float_num's Precision, MIN_PRECISION and MAX_PRECISION"; r=1; }
  grep -nP '^(?!#).*(\t1(1[7-9]|2[0-3])-P-|\.\./)|^src/precision\.rs:' "$d/$MAN" && { echo "  a bound mutant row names the retired precision ids, src/precision.rs or a ../ path (git apply skips it)"; r=1; }
  return $r
}
check "$root" || bad "precision items are not one zf definition"
t=$(mktemp -d "${TMPDIR:-$HOME/.cache}/maxprec.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
plant() {  # FILE OLD-REGEX NEW
  mkdir -p "$t/$(dirname "$DEFS")" "$t/$(dirname "$BP")" "$t/$(dirname "$MAN")"
  cp "$root/$DEFS" "$t/$DEFS"; cp "$root/$BP" "$t/$BP"; cp "$root/$MAN" "$t/$MAN"
  sed -i -E "s/$2/$3/" "$t/$1"
  cmp -s "$root/$1" "$t/$1" && { bad "plant on $1 changed nothing ($2)"; return; }
  if check "$t" > /dev/null; then bad "plant on $1 ($2 -> $3) was not caught"; fi
}
plant "$DEFS" '^pub const MAX_PRECISION: Precision = .*;$' 'pub const MAX_PRECISION: Precision = (usize::MAX \/ 4) - 1;'
plant "$BP" '^pub use zenith_float_num::.*;$' 'pub const MAX_PRECISION: Precision = (1 << 60) - 1;'
plant "$MAN" '^(src\/mag\/mod\.rs:Mag::finite\t001-F2)' '..\/zenith-float-num\/src\/defs.rs:MAX_PRECISION\t119-P-maxsub\t119-P-maxsub.patch\tball\n\1'
[ $fail = 0 ] && { echo "test_max_precision_single_def: ok"; exit 0; }
exit 1
