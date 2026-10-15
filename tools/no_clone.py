#!/usr/bin/env python3
"""Zero-clone gate for zenith-float library code.

Fails (exit 1) on any of these outside #[cfg(test)] items and tests/benches/examples:
  .clone()  .to_vec()  .cloned()  Clone::clone  .clone_from(  impl Clone for  #[derive(.. Clone ..)] without Copy
Exceptions live in tools/clone_allowlist.txt as `path:line-independent-key` lines
(`path<TAB>fn<TAB>reason<TAB>approval`), each needing Jeff's approval. The file is empty by policy.
"""
import os, re, sys, collections

ROOT = os.environ.get("NO_CLONE_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIRS = ["zenith-float-num/src", "zenith-float-macro/src", "src", "zenith-float-compare/src"]
TOKENS = [
    ("clone()", re.compile(r"\.clone\(\)")),
    ("to_vec()", re.compile(r"\.to_vec\(\)")),
    ("cloned()", re.compile(r"\.cloned\(\)")),
    ("Clone::clone", re.compile(r"\bClone::clone\b")),
    ("clone_from", re.compile(r"\.clone_from\(")),
    ("impl Clone", re.compile(r"\bimpl(?:<[^>]*>)?\s+Clone\s+for\b")),
]
DERIVE = re.compile(r"#\[derive\(([^)]*)\)\]")


def test_mask(lines):
    skip = [False] * len(lines)
    i = 0
    while i < len(lines):
        if re.match(r"\s*#\[cfg\(test\)\]", lines[i]):
            j, depth, seen = i, 0, False
            while j < len(lines):
                depth += lines[j].count("{") - lines[j].count("}")
                seen = seen or "{" in lines[j]
                skip[j] = True
                if (seen and depth <= 0) or (not seen and lines[j].rstrip().endswith(";")):
                    break
                j += 1
            i = j + 1
        else:
            i += 1
    return skip


def load_allow():
    path = os.path.join(ROOT, "tools", "clone_allowlist.txt")
    out = set()
    if os.path.exists(path):
        for l in open(path):
            l = l.strip()
            if l and not l.startswith("#"):
                parts = l.split("\t")
                out.add((parts[0], parts[1]))
    return out


def main():
    allow = load_allow()
    hits = []
    for d in DIRS:
        for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
            dn[:] = [x for x in dn if x not in ("target", ".git")]
            for f in sorted(fn):
                if not f.endswith(".rs"):
                    continue
                p = os.path.join(dp, f)
                rel = os.path.relpath(p, ROOT)
                if rel.endswith("ops/tests.rs"):
                    continue
                lines = open(p, errors="replace").read().split("\n")
                skip = test_mask(lines)
                cur = "?"
                for i, l in enumerate(lines):
                    m = re.search(r"\bfn\s+(\w+)", l)
                    if m:
                        cur = m.group(1)
                    if skip[i]:
                        continue
                    code = l.split("//")[0]
                    for name, rx in TOKENS:
                        for _ in rx.finditer(code):
                            hits.append((rel, i + 1, cur, name))
                    dm = DERIVE.search(code)
                    if dm:
                        items = [x.strip() for x in dm.group(1).split(",")]
                        if "Clone" in items and "Copy" not in items:
                            hits.append((rel, i + 1, cur, "derive(Clone)"))
    bad = [h for h in hits if (h[0], h[2]) not in allow]
    by = collections.Counter(h[3] for h in bad)
    per = collections.Counter(h[0] for h in bad)
    if "--list" in sys.argv:
        for h in bad:
            print("%s:%d [%s] %s" % h)
    print("no_clone: %d violation(s) %s" % (len(bad), dict(by)))
    if "--files" in sys.argv:
        for f, n in per.most_common():
            print("%5d %s" % (n, f))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
