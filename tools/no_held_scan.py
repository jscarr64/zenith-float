#!/usr/bin/env python3
"""no_held_scan: no held, xfail, skip or inverted-expectation list in the repo (Jeff 2026-11-25 3:36 AM CT brief, item 3).

A known-wrong result is fixed red-first, never parked: a test does not leave cells out of its digest or its oracle check
(a held counter or list), does not run every mode but one (a None-gap mode array, a mode `continue`), does not pin a
wrong answer as expected (a KNOWN_WRONG / known-wrong list, "today it gives" wording), and is not ignored or xfail.

Scans every git-tracked *.rs, *.py and *.sh file, comments included (held markers are comments). Skipped:
  - HOLD files (first line `// HOLD Hnn:` or `# HOLD Hnn:`): held code that is not compiled, under HOLD.md (Jeff's
    HOLD method); the index itself is markdown and not scanned;
  - this scan, its selftest and tools/test_rounding_mode_matrices.sh (they name the patterns they refuse).
FAILS (file:line [kind] text) on:
  held        a HELD marker (the word in capitals), a held_* counter or list (`held_x +=`, `held_x[..]`, `let mut held`),
              dd_held, grid_held
  none-gap    NONE_GAP_MODES, a `NONE-GAP F-x:` label, a `MODE-SUBSET gap` label, KNOWN GAP / KNOWN-GAP
  inverted    KNOWN_WRONG, known-wrong, "today ... gives/returns" or "currently ... gives/returns" wording
  xfail       xfail, an `#[ignore]` attribute, `#[should_panic]`
  mode-skip   `if rm == RoundingMode::X { continue; }` (or matches!(rm, ..)): a mode left out of a sweep
Past-tense history (`was ...`, `(NONE-GAP F-B)` in a re-pin note) does not match.
usage: tools/no_held_scan.sh [--root DIR]   (exit 0 clean, 1 findings, 2 usage)
"""
import re
import subprocess
import sys

LINE_RULES = [
    ("held", re.compile(r"\bHELD\b")),
    ("held", re.compile(r"\bheld_[a-z0-9_]+\s*(\+=|-=|\[|=[^=>])")),
    ("held", re.compile(r"\blet\s+(mut\s+)?held\b")),
    ("held", re.compile(r"\bdd_held\b|\bgrid_held\b")),
    ("none-gap", re.compile(r"\bNONE_GAP_MODES\b|NONE-GAP F-\w+:|MODE-SUBSET gap|\bKNOWN[ _-]GAPS?\b")),
    ("inverted", re.compile(r"\bKNOWN_WRONG\w*|(?i:\bknown-wrong\b)")),
    ("inverted", re.compile(r"(?i)\b(today|currently)\b[^.;]{0,60}\b(gives|returns)\b")),
    ("inverted", re.compile(r"(?i)\b(gives|returns)\b[^.;]{0,60}\b(today|currently)\b")),
    ("xfail", re.compile(r"(?i)\bxfail\b|^\s*#\[(ignore|should_panic)\b")),
]
MODE_SKIP = re.compile(
    r"if\s+(\w+\s*==\s*RoundingMode::\w+|matches!\(\s*\w+\s*,\s*RoundingMode::[^)]*\))\s*\{\s*continue;\s*\}")
# the scan, its selftest and the mode-matrix guard name the patterns they refuse
SELF = {"tools/no_held_scan.py", "tools/no_held_scan.sh", "tools/no_held_selftest.sh", "tools/test_rounding_mode_matrices.sh"}
HOLD = re.compile(r"^(//|#) HOLD H\d+:")


def main(argv):
    root = "."
    if argv[:1] == ["--root"] and len(argv) == 2:
        root = argv[1]
    elif argv:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    out = subprocess.run(["git", "-C", root, "ls-files", "-z", "*.rs", "*.py", "*.sh"], capture_output=True, check=True)
    found = 0
    for f in out.stdout.decode().split("\0"):
        if not f or f in SELF:
            continue
        try:
            text = open(f"{root}/{f}", encoding="utf-8").read()
        except (FileNotFoundError, UnicodeDecodeError):
            continue
        if HOLD.match(text):
            continue
        for n, line in enumerate(text.split("\n"), 1):
            for kind, rx in LINE_RULES:
                if rx.search(line):
                    print(f"{f}:{n} [{kind}] {line.strip()[:160]}")
                    found += 1
                    break
        for m in MODE_SKIP.finditer(text):
            n = text.count("\n", 0, m.start()) + 1
            print(f"{f}:{n} [mode-skip] {' '.join(m.group(0).split())[:160]}")
            found += 1
    if found:
        print(f"no_held_scan: FAIL {found} finding(s); fix each red-first (see the header of tools/no_held_scan.py)")
        return 1
    print("no_held_scan: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
