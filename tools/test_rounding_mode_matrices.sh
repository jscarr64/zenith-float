#!/usr/bin/env bash
# Red test (Spence #39 F9, ruling C13; P1 "every rounding mode (all 7)"): in zenith-float-num every named rounding-mode
# matrix (`[RoundingMode; N]`) lists all 7 modes (None, Up, Down, ToZero, FromZero, ToEven, ToOdd), each once, and no
# text in the crate says the matrix has 6 modes (Spence's F9 CHECK: '\b(6|six) (rounding )?modes').
# Inline mode loops and golden-backed mode sets are outside this check (rulings pending, fold-pr1-139).
# The check also runs on a planted copy and must catch it.   usage: tools/test_rounding_mode_matrices.sh
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
arr = re.compile(r"\[RoundingMode;\s*(\d+)\]\s*=\s*\[(.*?)\];", re.S)
text = re.compile(r"\b(6|six) (rounding )?modes", re.I)
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
            for m in arr.finditer(s):
                line = s.count("\n", 0, m.start()) + 1
                modes = re.findall(r"RoundingMode::(\w+)", m.group(2))
                if int(m.group(1)) != 7 or len(modes) != 7 or set(modes) != want:
                    print(f"  {rel}:{line}: [RoundingMode; {m.group(1)}] has {sorted(modes)}, missing {sorted(want - set(modes))}")
                    bad = 1
        for i, l in enumerate(s.splitlines(), 1):
            if text.search(l):
                print(f"  {rel}:{i}: says 6 modes: {l.strip()[:120]}")
                bad = 1
sys.exit(bad)
PY
  return $r
}
check "$root" || bad "a rounding-mode matrix in zenith-float-num does not cover all 7 modes"
t=$(mktemp -d "${TMPDIR:-$HOME/.cache}/rmmat.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/zenith-float-num/tests"
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
[ $fail = 0 ] && { echo "test_rounding_mode_matrices: ok"; exit 0; }
exit 1
