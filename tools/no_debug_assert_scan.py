#!/usr/bin/env python3
"""no_debug_assert_scan: no Rust code in the repo uses debug_assert!, debug_assert_eq! or debug_assert_ne! (Jeff 2026-11-24
4:40 PM CT item 5: a debug_assert checks only in a debug build, so a release build checked less than a debug build did).
Shipped code refuses with an always-on Result (P2), proves the case unreachable, or makes it unrepresentable by type; test
code uses assert!/assert_eq!, which run in every build.

Scans every git-tracked *.rs file. `//` comments (line and trailing) are skipped.
FAILS (file:line [debug-assert] text) on any debug_assert!/debug_assert_eq!/debug_assert_ne! macro call.
usage: tools/no_debug_assert_scan.sh [--root DIR]   (exit 0 clean, 1 findings, 2 usage)
"""
import re
import subprocess
import sys

MACRO = re.compile(r"\bdebug_assert(_eq|_ne)?\s*!")


def code_part(line):
    """The line without string literal contents and without a trailing // comment."""
    out, in_str, i = [], False, 0
    while i < len(line):
        c = line[i]
        if c == '"' and (i == 0 or line[i - 1] != "\\"):
            in_str = not in_str
        if not in_str and line.startswith("//", i):
            break
        if not in_str or c == '"':
            out.append(c)
        i += 1
    return "".join(out)


def main(argv):
    root = "."
    if argv[:1] == ["--root"] and len(argv) == 2:
        root = argv[1]
    elif argv:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    files = subprocess.run(["git", "-C", root, "ls-files", "-z", "*.rs"], capture_output=True, check=True).stdout
    found = 0
    for f in files.decode().split("\0"):
        if not f:
            continue
        try:
            text = open(f"{root}/{f}", encoding="utf-8").read()
        except (FileNotFoundError, UnicodeDecodeError):
            continue
        for n, line in enumerate(text.split("\n"), 1):
            if MACRO.search(code_part(line)):
                print(f"{f}:{n} [debug-assert] {line.strip()}")
                found += 1
    if found:
        print(f"no_debug_assert_scan: FAIL {found} finding(s); use an always-on check (see the header of tools/no_debug_assert_scan.py)")
        return 1
    print("no_debug_assert_scan: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
