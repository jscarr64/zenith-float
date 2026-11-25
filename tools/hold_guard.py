#!/usr/bin/env python3
"""hold_guard: no HOLD file changes without a recorded Jeff ruling (P11 / HOLD.md; Paul 2026-11-25 9:31 AM CT, after the
H05 breach of 4e1d8e9, WORK2.md RULE BREACH 2026-11-25 9:32 AM CT).

tools/hold_manifest.txt pins every held file: one line `KEY  PATH  RULING`, KEY one of
  <sha256>   the file is held with exactly this content (HOLD.md is pinned the same way);
  live       the file was taken live (its HOLD header removed);
  deleted    the file was deleted;
and RULING the Jeff ruling that allows that state; it must name `Jeff YYYY-MM-DD`.
FAILS (exit 1) when
  - a tracked file whose first line is `// HOLD Hnn:` or `# HOLD Hnn:`, or HOLD.md, has no manifest line;
  - a pinned file's content differs from its sha256 (an edit), is missing, or lost its HOLD header;
  - a `live` file still has a HOLD header or is missing; a `deleted` file exists;
  - a line is malformed, repeats a path, or its RULING does not name `Jeff YYYY-MM-DD`.
So a held file only changes together with a manifest line that cites Jeff's ruling for the new state.
usage: tools/hold_guard.sh [--root DIR]   (exit 0 clean, 1 findings, 2 usage)
"""
import hashlib
import re
import subprocess
import sys

MANIFEST = "tools/hold_manifest.txt"
HOLD = re.compile(rb"^(//|#) HOLD H\d+:")
RULING = re.compile(r"\bJeff \d{4}-\d{2}-\d{2}\b")
KEY = re.compile(r"^([0-9a-f]{64}|live|deleted)$")


def main(argv):
    root = "."
    if argv[:1] == ["--root"] and len(argv) == 2:
        root = argv[1]
    elif argv:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    out = subprocess.run(["git", "-C", root, "ls-files", "-z"], capture_output=True, check=True)
    tracked = [f for f in out.stdout.decode().split("\0") if f]
    content = {}
    for f in tracked:
        try:
            content[f] = open(f"{root}/{f}", "rb").read()
        except (FileNotFoundError, IsADirectoryError):
            continue
    held = {f for f, b in content.items() if HOLD.match(b)}
    if "HOLD.md" in content:
        held.add("HOLD.md")
    bad = []
    pinned = {}
    try:
        lines = open(f"{root}/{MANIFEST}", encoding="utf-8").read().splitlines()
    except FileNotFoundError:
        lines = []
        if held:
            bad.append(f"{MANIFEST}: missing, {len(held)} held file(s) unpinned")
    for n, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split(None, 2)
        if len(parts) != 3 or not KEY.match(parts[0]):
            bad.append(f"{MANIFEST}:{n}: malformed line: {line[:120]}")
            continue
        key, path, ruling = parts
        if not RULING.search(ruling):
            bad.append(f"{MANIFEST}:{n}: {path}: the ruling does not name `Jeff YYYY-MM-DD`: {ruling[:120]}")
        if path in pinned:
            bad.append(f"{MANIFEST}:{n}: {path} listed twice")
        pinned[path] = (key, n)
        body = content.get(path)
        if key == "deleted":
            if body is not None:
                bad.append(f"{path}: recorded deleted ({MANIFEST}:{n}) but present")
        elif key == "live":
            if body is None:
                bad.append(f"{path}: recorded live ({MANIFEST}:{n}) but missing")
            elif HOLD.match(body):
                bad.append(f"{path}: recorded live ({MANIFEST}:{n}) but still has a HOLD header")
        elif body is None:
            bad.append(f"{path}: held ({MANIFEST}:{n}) but missing (a deletion needs a `deleted` line with Jeff's ruling)")
        else:
            if path != "HOLD.md" and not HOLD.match(body):
                bad.append(f"{path}: held ({MANIFEST}:{n}) but its HOLD header is gone (taking it live needs a `live` line)")
            got = hashlib.sha256(body).hexdigest()
            if got != key:
                bad.append(f"{path}: CHANGED (sha256 {got[:12]}.., pinned {key[:12]}.. at {MANIFEST}:{n}); a held "
                           "file changes only with a manifest line citing Jeff's ruling")
    for f in sorted(held - set(pinned)):
        bad.append(f"{f}: held file not in {MANIFEST}")
    for b in bad:
        print(f"hold_guard: {b}")
    if bad:
        print(f"hold_guard: FAIL {len(bad)} finding(s)")
        return 1
    print(f"hold_guard: ok ({len(pinned)} pinned)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
