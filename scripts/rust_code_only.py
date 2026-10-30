#!/usr/bin/env python3
"""Print Rust source with comments removed, as `path:line:code`, one line per source line (scripts/ci_no_intrinsics.sh).

Removes // line comments (including /// and //!) and /* */ block comments (nested, multi-line). String literals ("...", with
escapes, and raw r"..." / r#"..."#) and char literals are kept as code, so "//" inside a string does not start a comment.
Usage: rust_code_only.py FILE...   Exit 2 when a file cannot be read.
"""
import re
import sys

CHAR = re.compile(r"'(\\(x[0-9a-fA-F]{2}|u\{[0-9a-fA-F]{1,6}\}|.)|[^\\'\n])'")
RAW = re.compile(r'r(#*)"')


def strip(src):
    out, i, n, depth = [], 0, len(src), 0
    while i < n:
        c = src[i]
        if depth:
            if src.startswith("/*", i):
                depth, i = depth + 1, i + 2
            elif src.startswith("*/", i):
                depth, i = depth - 1, i + 2
            else:
                out.append("\n" if c == "\n" else " ")
                i += 1
            continue
        if src.startswith("//", i):
            j = src.find("\n", i)
            i = n if j < 0 else j
            continue
        if src.startswith("/*", i):
            depth, i = 1, i + 2
            out.append("  ")
            continue
        m = RAW.match(src, i) if c == "r" and (i == 0 or not (src[i - 1].isalnum() or src[i - 1] == "_")) else None
        if m:
            end = src.find('"' + m.group(1), m.end())
            j = n if end < 0 else end + 1 + len(m.group(1))
            out.append(src[i:j]); i = j
            continue
        if c == '"':
            j = i + 1
            while j < n and src[j] != '"':
                j += 2 if src[j] == "\\" else 1
            out.append(src[i:j + 1]); i = j + 1
            continue
        if c == "'":
            m = CHAR.match(src, i)
            if m:
                out.append(m.group(0)); i = m.end()
                continue
        out.append(c)
        i += 1
    return "".join(out)


def main(paths):
    rc = 0
    for p in paths:
        try:
            with open(p, encoding="utf-8") as fh:
                src = fh.read()
        except OSError as e:
            print("rust_code_only: cannot read %s: %s" % (p, e), file=sys.stderr)
            rc = 2
            continue
        for k, line in enumerate(strip(src).split("\n"), 1):
            if line.strip():
                print("%s:%d:%s" % (p, k, line))
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
