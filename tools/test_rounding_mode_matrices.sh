#!/usr/bin/env bash
# Red test (Spence #39 F9, ruling C13; P1 "every rounding mode (all 7)"): in zenith-float-num every named rounding-mode
# matrix (`[RoundingMode; N]`) lists all 7 modes (None, Up, Down, ToZero, FromZero, ToEven, ToOdd), each once, and no
# text in the crate says the matrix has 6 modes (Spence's F9 CHECK: '\b(6|six) (rounding )?modes').
# There is no exception: the ZF-M27 None gaps (F-A, F-B, F-C) are fixed, and the KNOWN-GAP list and its NONE_GAP_MODES
# arrays are gone (b20). An array or a use named NONE_GAP_MODES, or a `NONE-GAP F-x:` label, anywhere fails. Every
# inline mode loop (`in [RoundingMode::A, ...]`) lists all 7 modes or carries, in the 3 lines above it, a label
# `MODE-SUBSET mode-specific: <why>` (the expected value holds only for these modes; Spence's C13 loop ruling); the old
# `MODE-SUBSET gap F-x` label is refused. Golden-backed mode sets (generators) are outside this check.
# The check also runs on a planted copy and must catch it.
# usage: tools/test_rounding_mode_matrices.sh
set -u
root=$(cd "$(dirname "$0")/.." && pwd)
fail=0; bad() { echo "FAIL: $*"; fail=1; }
command -v python3 > /dev/null || { echo "FAIL: python3 not found"; exit 2; }
check() {  # DIR (a git checkout or a planted copy with the same layout) -> 0 when every named matrix has all 7 modes
  local d=$1 r=0
  python3 -B - "$d" <<'PY' || r=1
import os, re, sys
d = sys.argv[1]
want = {"None", "Up", "Down", "ToZero", "FromZero", "ToEven", "ToOdd"}
GAP = "NONE_GAP_MODES"
arr = re.compile(r"((?:[ \t]*///[^\n]*\n)*)[^\n]*?\b(\w+): \[RoundingMode;\s*(\d+)\]\s*=\s*\[(.*?)\];", re.S)
text = re.compile(r"\b(6|six) (rounding )?modes", re.I)
use = re.compile(r"\b" + GAP + r"\b")
label = re.compile(r"NONE-GAP F-[ABC]:")
inline = re.compile(r"\bin\s*\[\s*((?:RoundingMode::\w+\s*,?\s*)+)\]")
subset = re.compile(r"MODE-SUBSET mode-specific: \S")
bad = 0
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
                if n != 7 or len(modes) != 7 or set(modes) != want:
                    print(f"  {rel}:{line}: [RoundingMode; {n}] has {sorted(modes)}, missing {sorted(want - set(modes))}")
                    bad = 1
            for m in inline.finditer(s):
                modes = re.findall(r"RoundingMode::(\w+)", m.group(1))
                if len(modes) == 7 and set(modes) == want:
                    continue
                line = s.count("\n", 0, m.start()) + 1
                above = "\n".join(lines[max(0, line - 4):line])
                if not subset.search(above):
                    print(f"  {rel}:{line}: inline loop over {modes} has no 'MODE-SUBSET mode-specific: why' label")
                    bad = 1
            for m in use.finditer(s):
                line = s.count("\n", 0, m.start()) + 1
                print(f"  {rel}:{line}: {GAP}: the None gaps are fixed (ZF-M27); check every mode (MODES)")
                bad = 1
        for i, l in enumerate(s.splitlines(), 1):
            if label.search(l) or "MODE-SUBSET gap" in l:
                print(f"  {rel}:{i}: a None-gap label: {l.strip()[:120]}")
                bad = 1
        for i, l in enumerate(s.splitlines(), 1):
            if text.search(l):
                print(f"  {rel}:{i}: says 6 modes: {l.strip()[:120]}")
                bad = 1
sys.exit(bad)
PY
  return $r
}
check "$root" || bad "a rounding-mode matrix in zenith-float-num does not cover all 7 modes (or keeps a None gap)"
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/test-rounding-mode-matrices"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp (Spence n3)
t=$(mktemp -d "$zf_tmp/rmmat.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
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
check "$t" > /dev/null && bad "a planted NONE_GAP_MODES gap (the old well-formed form) was not caught"
rm "$g"
printf 'fn f(rm: RoundingMode) {\n    // NONE-GAP F-A: exact pow results\n    g(rm);\n}\n' > "$g"
check "$t" > /dev/null && bad "a planted NONE-GAP label was not caught"
printf 'fn f() {\n    // MODE-SUBSET gap F-B: None overflow\n    for rm in [RoundingMode::Up, RoundingMode::Down] {\n        g(rm);\n    }\n}\n' > "$g"
check "$t" > /dev/null && bad "a planted 'MODE-SUBSET gap' loop label was not caught"
rm "$g"
printf 'fn f() {\n    for rm in [RoundingMode::Up, RoundingMode::Down] {\n        g(rm);\n    }\n}\n' > "$t/zenith-float-num/tests/loop.rs"
check "$t" > /dev/null && bad "an unlabelled 2-mode inline loop was not caught"
printf 'fn f() {\n    // MODE-SUBSET mode-specific: Up and Down bracket the value\n    for rm in [RoundingMode::Up, RoundingMode::Down] {\n        g(rm);\n    }\n}\n' > "$t/zenith-float-num/tests/loop.rs"
check "$t" > /dev/null || bad "a labelled inline loop was flagged"
sed -i 's/MODE-SUBSET mode-specific:/MODE-SUBSET because/' "$t/zenith-float-num/tests/loop.rs"
check "$t" > /dev/null && bad "an inline loop with a malformed label was not caught"
rm "$t/zenith-float-num/tests/loop.rs"
[ $fail = 0 ] && { echo "test_rounding_mode_matrices: ok"; exit 0; }
exit 1
