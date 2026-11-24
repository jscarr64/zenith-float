#!/usr/bin/env python3
"""no_tmp_scan: no script in the repo puts work files in /tmp (Jeff 2026-10-23/24: a reboot wiped /tmp and with it work dirs,
a dependency pin and a commit; Spence n3 on #40 e7c6d1a; Jeff 2026-11-24 4:40 PM CT item 8). Temp files and dirs go under
the job scratch dir "${ZF_SCRATCH:-$HOME/.cache/zenith-float}/JOB" (root scripts) or "$(zb_scratch)"
(zenith-float-bound/tools/scratch.sh).

Scans every git-tracked shell script (*.sh, *.bash, .githooks/*) and Python file (*.py) in the tree. Comment lines are skipped.
FAILS (file:line [form] text) on:
  [tmp-path]      a /tmp path in code (`/tmp`, `/tmp/x`, `${TMPDIR:-/tmp}`); only the shared cross-agent lock files
                  /tmp/pop-deep.lock, /tmp/pop-gate.lock and /tmp/pop-gate.slot* are allowed (a lock, not a work dir: moving
                  one would split the lock between agents)
  [mktemp-bare]   mktemp with no template path (`mktemp`, `mktemp -d`): it defaults to $TMPDIR or /tmp
  [mktemp-tmpdir] a mktemp template under $TMPDIR (`"${TMPDIR:-...}/x.XXXXXX"`): TMPDIR is /tmp on most machines
  [py-tempfile]   a Python tempfile call (mkdtemp, mkstemp, TemporaryDirectory, TemporaryFile, NamedTemporaryFile,
                  SpooledTemporaryFile) without dir=
usage: tools/no_tmp_scan.sh [--root DIR]   (exit 0 clean, 1 findings, 2 usage)
"""
import os
import re
import subprocess
import sys

TMP_PATH = re.compile(r"""(^|[\s"'=:(\-])/tmp([/\s"'$}):;]|$)""")
LOCKS = re.compile(r"/tmp/pop-(deep\.lock|gate\.lock|gate\.slot)")
MKTEMP = re.compile(r"(\$\(|^|[;&|]|\bthen|\bdo|\belse)\s*mktemp\b(?:\s+-[A-Za-z]+)*\s*(\)|;|&&|\|\||$)")
MKTEMP_TMPDIR = re.compile(r"(\$\(|^|[;&|])\s*mktemp\b(?:\s+-[A-Za-z]+)*\s+\"?\$\{?TMPDIR")
PY_TEMP = re.compile(r"\btempfile\.(mkdtemp|mkstemp|TemporaryDirectory|TemporaryFile|NamedTemporaryFile|SpooledTemporaryFile)\(")
SELF = {"tools/no_tmp_scan.py", "tools/no_tmp_selftest.sh", "zenith-float-bound/tools/no_tmp_selftest.sh"}


def files(root):
    out = subprocess.run(["git", "-C", root, "ls-files", "-z"], capture_output=True, check=True).stdout.decode()
    for f in out.split("\0"):
        if not f or f in SELF:
            continue
        if f.endswith((".sh", ".bash", ".py")) or f.startswith(".githooks/"):
            yield f


def py_call_has_dir(lines, i, start):
    """True when the tempfile call opened at lines[i][start:] passes dir= before its closing paren."""
    depth, j, k = 0, i, start
    while j < len(lines):
        s = lines[j]
        while k < len(s):
            c = s[k]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    return False
            elif depth == 1 and s.startswith("dir=", k) and (k == 0 or not (s[k - 1].isalnum() or s[k - 1] == "_")):
                return True
            k += 1
        j, k = j + 1, 0
    return False


def scan(root):
    hits = []
    for f in files(root):
        p = os.path.join(root, f)
        try:
            lines = open(p, encoding="utf-8", errors="replace").read().split("\n")
        except (FileNotFoundError, IsADirectoryError):
            continue
        py = f.endswith(".py")
        for i, line in enumerate(lines):
            code = line.lstrip()
            if code.startswith("#") or (i == 0 and code.startswith("#!")):
                continue
            if not py:
                # a trailing shell comment is not code
                code = re.sub(r"(^|\s)#[^'\"]*$", "", line)
            if TMP_PATH.search(LOCKS.sub("", code)):
                hits.append((f, i + 1, "tmp-path", line.strip()))
            if not py and MKTEMP.search(code):
                hits.append((f, i + 1, "mktemp-bare", line.strip()))
            if not py and MKTEMP_TMPDIR.search(code):
                hits.append((f, i + 1, "mktemp-tmpdir", line.strip()))
            if py:
                for m in PY_TEMP.finditer(code):
                    if not py_call_has_dir(lines, i, line.index(m.group(0)) + len(m.group(0)) - 1):
                        hits.append((f, i + 1, "py-tempfile", line.strip()))
    return hits


def main():
    a = sys.argv[1:]
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if a[:1] == ["--root"] and len(a) == 2:
        root = a[1]
    elif a:
        print(__doc__.strip().split("\n")[-1], file=sys.stderr)
        return 2
    hits = scan(root)
    for f, ln, form, text in hits:
        print(f"{f}:{ln} [{form}] {text[:200]}")
    if hits:
        print(f"no_tmp_scan: FAIL {len(hits)} finding(s); put temp files under \"${{ZF_SCRATCH:-$HOME/.cache/zenith-float}}/JOB\" or \"$(zb_scratch)\"")
        return 1
    print("no_tmp_scan: ok (no /tmp work files, no bare or TMPDIR mktemp, no tempfile without dir=)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
