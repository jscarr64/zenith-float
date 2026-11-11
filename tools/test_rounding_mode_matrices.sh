#!/usr/bin/env bash
# Red test (Spence #39 F9, ruling C13; P1 "every rounding mode (all 7)"): in zenith-float-num every named rounding-mode
# matrix (`[RoundingMode; N]`) lists all 7 modes (None, Up, Down, ToZero, FromZero, ToEven, ToOdd), each once, and no
# text in the crate says the matrix has 6 modes (Spence's F9 CHECK: '\b(6|six) (rounding )?modes').
# The one exception is the KNOWN-GAP list below (Paul/Jeff/Spence 2026-11-11, #39 final rule): a family that fails the
# strict None oracle (MPFR toward zero at rp(p,64), FIXLIST ZF-M27) keeps its None rows out, through an array named
# NONE_GAP_MODES holding exactly the 6 other modes. Such an array may exist only in a listed file. Its doc comment names
# every finding id listed for that file (`NONE-GAP F-A`...). Every use (`in NONE_GAP_MODES`, `&NONE_GAP_MODES`) carries a
# `NONE-GAP F-x:` comment within the 3 lines above it. A listed file without the array is a stale entry and fails, so the
# ZF-M27 conformance PR has to delete each entry when it fixes the family. Inline mode loops and golden-backed mode sets
# are outside this check. The check also runs on a planted copy and must catch it.
# usage: tools/test_rounding_mode_matrices.sh
set -u
root=$(cd "$(dirname "$0")/.." && pwd)
fail=0; bad() { echo "FAIL: $*"; fail=1; }
command -v python3 > /dev/null || { echo "FAIL: python3 not found"; exit 2; }
check() {  # DIR (a git checkout or a planted copy with the same layout) -> 0 when every named matrix has all 7 modes
  local d=$1 r=0
  python3 - "$d" <<'PY' || r=1
import os, re, sys
d = sys.argv[1]
want = {"None", "Up", "Down", "ToZero", "FromZero", "ToEven", "ToOdd"}
# KNOWN-GAP list: file -> finding ids (fold-pr1-147; acceptance: ~/pop-logs/spence/acceptance/
# zenith-float-rusty-zf-m27-none-conformance.md). F-A: an exact pow/nth_root result comes back as its predecessor,
# flagged inexact. F-B: None overflow gives +-Inf, and Jeff ruled +-MAX. F-C: log with a base near a grid point is
# 1 ulp below the toward-zero cut.
KNOWN_GAPS = {
    "zenith-float-num/tests/mpfr/exact_bit_functions.rs": ["F-A", "F-B"],
    "zenith-float-num/tests/mpfr/wide_primitives.rs": ["F-A", "F-B"],
    "zenith-float-num/tests/mpfr/retry_budget.rs": ["F-A", "F-C"],
    "zenith-float-num/tests/mpfr/subnormal_ops.rs": ["F-B"],
    "zenith-float-num/tests/exponent_range.rs": ["F-B"],
    "zenith-float-num/tests/r13_props.rs": ["F-B"],
    "zenith-float-num/tests/sinh_cosh_overflow.rs": ["F-B"],
}
GAP = "NONE_GAP_MODES"
arr = re.compile(r"((?:[ \t]*///[^\n]*\n)*)[^\n]*?\b(\w+): \[RoundingMode;\s*(\d+)\]\s*=\s*\[(.*?)\];", re.S)
text = re.compile(r"\b(6|six) (rounding )?modes", re.I)
use = re.compile(r"(\bin\s+|&)" + GAP + r"\b")
label = re.compile(r"NONE-GAP (F-[ABC]): \S")
bad = 0
seen = set()
for top, _dirs, files in os.walk(os.path.join(d, "zenith-float-num")):
    parts = os.path.relpath(top, d).split(os.sep)
    if any(x == "target" or x.startswith(".") for x in parts):
        continue
    for name in files:
        if not name.endswith((".rs", ".md")):
            continue
        path = os.path.join(top, name)
        rel = os.path.relpath(path, d)
        s = open(path, encoding="utf-8").read()
        if name.endswith(".rs"):
            lines = s.splitlines()
            for m in arr.finditer(s):
                line = s.count("\n", 0, m.start(2)) + 1
                doc, ident, n, body = m.group(1), m.group(2), int(m.group(3)), m.group(4)
                modes = re.findall(r"RoundingMode::(\w+)", body)
                if ident == GAP:
                    ids = KNOWN_GAPS.get(rel)
                    if ids is None:
                        print(f"  {rel}:{line}: {GAP} outside the KNOWN-GAP list")
                        bad = 1
                        continue
                    seen.add(rel)
                    if n != 6 or len(modes) != 6 or set(modes) != want - {"None"}:
                        print(f"  {rel}:{line}: {GAP} must hold exactly the 6 modes other than None, has {sorted(modes)}")
                        bad = 1
                    for i in ids:
                        if f"NONE-GAP {i}" not in doc:
                            print(f"  {rel}:{line}: {GAP} doc does not name its finding {i} (NONE-GAP {i})")
                            bad = 1
                    continue
                if n != 7 or len(modes) != 7 or set(modes) != want:
                    print(f"  {rel}:{line}: [RoundingMode; {n}] has {sorted(modes)}, missing {sorted(want - set(modes))}")
                    bad = 1
            for m in use.finditer(s):
                line = s.count("\n", 0, m.start()) + 1
                above = "\n".join(lines[max(0, line - 4):line])
                if not label.search(above):
                    print(f"  {rel}:{line}: use of {GAP} without a 'NONE-GAP F-x: reason' label in the 3 lines above")
                    bad = 1
        for i, l in enumerate(s.splitlines(), 1):
            if text.search(l):
                print(f"  {rel}:{i}: says 6 modes: {l.strip()[:120]}")
                bad = 1
if os.path.isdir(os.path.join(d, ".git")) or os.path.isfile(os.path.join(d, ".git")):
    for rel in sorted(set(KNOWN_GAPS) - seen):
        print(f"  {rel}: stale KNOWN-GAP entry (no {GAP} array): delete the entry")
        bad = 1
sys.exit(bad)
PY
  return $r
}
check "$root" || bad "a rounding-mode matrix in zenith-float-num does not cover all 7 modes (or breaks the KNOWN-GAP rules)"
t=$(mktemp -d "${TMPDIR:-$HOME/.cache}/rmmat.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/zenith-float-num/tests/mpfr"
printf 'const M: [RoundingMode; 7] = [\n    RoundingMode::None,\n    RoundingMode::Up,\n    RoundingMode::Down,\n    RoundingMode::ToZero,\n    RoundingMode::FromZero,\n    RoundingMode::ToEven,\n    RoundingMode::ToOdd,\n];\n' > "$t/zenith-float-num/tests/good.rs"
check "$t" > /dev/null || bad "the 7-mode control copy was flagged"
sed 's/    RoundingMode::None,//; s/; 7\]/; 6]/' "$t/zenith-float-num/tests/good.rs" > "$t/zenith-float-num/tests/six.rs"
check "$t" > /dev/null && bad "a planted 6-mode array was not caught"
rm "$t/zenith-float-num/tests/six.rs"
sed 's/RoundingMode::ToOdd/RoundingMode::ToEven/' "$t/zenith-float-num/tests/good.rs" > "$t/zenith-float-num/tests/dup.rs"
check "$t" > /dev/null && bad "a planted 7-entry array with a repeated mode was not caught"
rm "$t/zenith-float-num/tests/dup.rs"
echo '// checked in all six modes' > "$t/zenith-float-num/tests/text.rs"
check "$t" > /dev/null && bad "planted 'six modes' text was not caught"
rm "$t/zenith-float-num/tests/text.rs"
g="$t/zenith-float-num/tests/sinh_cosh_overflow.rs"
gap() {  # DOC LABEL -> a planted gap file with that doc line and use label
  printf '%s\nconst NONE_GAP_MODES: [RoundingMode; 6] = [\n    RoundingMode::Up,\n    RoundingMode::Down,\n    RoundingMode::ToZero,\n    RoundingMode::FromZero,\n    RoundingMode::ToEven,\n    RoundingMode::ToOdd,\n];\nfn f() {\n    %s\n    for rm in NONE_GAP_MODES {\n        g(rm);\n    }\n}\n' "$1" "$2" > "$g"
}
gap '/// NONE-GAP F-B: None overflow gives +-Inf' '// NONE-GAP F-B: sinh/cosh overflow under None'
check "$t" > /dev/null || bad "a well-formed listed gap was flagged"
gap '/// NONE-GAP F-B: None overflow gives +-Inf' '// no label here'
check "$t" > /dev/null && bad "an unlabelled use of NONE_GAP_MODES was not caught"
gap '/// a gap with no finding id' '// NONE-GAP F-B: sinh/cosh overflow under None'
check "$t" > /dev/null && bad "a gap array whose doc lacks its finding id was not caught"
gap '/// NONE-GAP F-B: None overflow gives +-Inf' '// NONE-GAP F-B: sinh/cosh overflow under None'
sed -i 's/RoundingMode::ToOdd,/RoundingMode::None,/' "$g"
check "$t" > /dev/null && bad "a gap array that is not exactly the 6 non-None modes was not caught"
mv "$g" "$t/zenith-float-num/tests/unlisted.rs"
check "$t" > /dev/null && bad "a NONE_GAP_MODES array in an unlisted file was not caught"
rm "$t/zenith-float-num/tests/unlisted.rs"
mkdir "$t/.git"
check "$t" > /dev/null && bad "a stale KNOWN-GAP entry (listed file without the array) was not caught"
[ $fail = 0 ] && { echo "test_rounding_mode_matrices: ok"; exit 0; }
exit 1
