#!/usr/bin/env python3
"""Zero-clone gate for zenith-float code (library AND test code: nothing is skipped).

Fails (exit 1) on any of these in the four library source dirs and in the test dirs (TEST_DIRS); `#[cfg(test)]` items,
ops/tests.rs, tests/, benches/ and the nostd-tests crate are counted too. tools/ and target dirs are never scanned:
  .clone()  .to_vec()  .cloned()  Clone::clone  .clone_from(  impl Clone for  #[derive(.. Clone ..)] without Copy
Exceptions live in tools/clone_allowlist.txt, one entry per (file, fn):
  path<TAB>fn<TAB>reason<TAB>approval<TAB>doc<TAB>sites
`sites` is the exact number of `.clone()` sites allowed in that fn, `doc` is a repo file that must contain the
section "Exception: constant caches" (tools/NO_CLONE_RULES.md). Only `.clone()` can be excepted. An entry that matches
no `.clone()` site is STALE and fails; sites above `sites` fail as ordinary violations; the total number of excepted
sites may never exceed tools/exception_count_baseline.txt (enforced by rules_scan.py). Entries need Jeff's approval.

Second counter (`--helpers`): uses of the hidden copy helpers (see HELPERS) in non-definition positions, plus
Rc::clone / Arc::clone in zenith-float-compare. Ratchet: tools/helper_copy_count_baseline.txt (rules_scan.py).
"""
import os, re, sys, collections

ROOT = os.environ.get("NO_CLONE_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIRS = ["zenith-float-num/src", "zenith-float-macro/src", "src", "zenith-float-compare/src"]
# Test code outside src/ is counted too (the integration tests, benches and the no_std test crate). DIRS stays the
# library set: rules_scan.py scopes the add_ref/sub_ref/mul_ref/div_ref ratchets to it.
TEST_DIRS = ["zenith-float-num/tests", "zenith-float-num/benches", "zenith-float-num/examples", "zenith-float-num/fuzz",
             "tests", "benches", "examples", "nostd-tests"]
ALL_DIRS = DIRS + TEST_DIRS
TOKENS = [
    ("clone()", re.compile(r"\.clone\(\)")),
    ("to_vec()", re.compile(r"\.to_vec\(\)")),
    ("cloned()", re.compile(r"\.cloned\(\)")),
    ("Clone::clone", re.compile(r"\bClone::clone\b")),
    ("clone_from", re.compile(r"\.clone_from\(")),
    ("impl Clone", re.compile(r"\bimpl(?:<[^>]*>)?\s+Clone\s+for\b")),
]
DERIVE = re.compile(r"#\[derive\(([^)]*)\)\]")

# --- helper-copy counter -------------------------------------------------------------------------------------
# Hidden copies: each of these is a `.clone()` (plus set_inexact(false) or a precision change) the clone scan cannot see.
HELPERS = ("exact_copy_ext", "exact_copy", "exact_operand", "copy_p", "svd_copy_prec", "lift_num")
HELPER_USE = re.compile(r"\b(%s)\b" % "|".join(HELPERS))
# definition positions, removed before counting: `fn NAME`, `macro_rules! NAME`, `use ...NAME...;`
HELPER_DEF = re.compile(r"\bfn\s+(?:%s)\b|\bmacro_rules!\s*(?:%s)\b" % ("|".join(HELPERS), "|".join(HELPERS)))
USE_LINE = re.compile(r"^\s*(?:pub(?:\([^)]*\))?\s+)?use\s")
HANDLE_COPY = re.compile(r"\b(Rc|Arc)::clone\b")  # handle shares in the compare crate, invisible to TOKENS
HANDLE_DIR = "zenith-float-compare/src"
STR = re.compile(r'"(?:\\.|[^"\\])*"')


def strip_code(line):
    """drop string literals, /* */ spans on one line and // comments"""
    s = STR.sub('""', line)
    s = re.sub(r"/\*.*?\*/", "", s)
    return s.split("//")[0]


def test_mask(lines):
    """Lines inside #[cfg(test)] items. NOT used by this scan any more (test code is counted); kept for
    rules_scan.py, whose diff scan still exempts added test lines from the non-clone rules (unwrap, panic, f64...)."""
    skip = [False] * len(lines)
    i = 0
    while i < len(lines):
        if re.match(r"\s*#\[cfg\((test\)|all\(test\b)", lines[i]):
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


def rs_files(dirs=None):
    """every .rs file of the scan (library dirs + test dirs; tools/ and target dirs are never walked)"""
    for d in (ALL_DIRS if dirs is None else dirs):
        for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
            dn[:] = sorted(x for x in dn if x not in ("target", ".git"))
            for f in sorted(fn):
                if f.endswith(".rs"):
                    p = os.path.join(dp, f)
                    yield os.path.relpath(p, ROOT), p


def doc_ok(doc):
    p = os.path.join(ROOT, doc)
    try:
        return "Exception: constant caches" in open(p, errors="replace").read()
    except OSError:
        return False


def load_allow():
    """-> (entries {(path, fn): dict}, errors [str])"""
    path = os.path.join(ROOT, "tools", "clone_allowlist.txt")
    out, errs = {}, []
    if os.path.exists(path):
        for n, l in enumerate(open(path), 1):
            l = l.rstrip("\n")
            if not l.strip() or l.lstrip().startswith("#"):
                continue
            parts = l.split("\t")
            where = "tools/clone_allowlist.txt:%d" % n
            if len(parts) < 6 or not all(x.strip() for x in parts[:6]):
                errs.append("%s: malformed entry (need path, fn, reason, approval, doc, sites)" % where)
                continue
            try:
                sites = int(parts[5])
                assert sites >= 1
            except (ValueError, AssertionError):
                errs.append("%s: sites must be an integer >= 1" % where)
                continue
            if not doc_ok(parts[4].strip()):
                errs.append("%s: doc %s is missing or has no 'Exception: constant caches' section" % (where, parts[4]))
                continue
            key = (parts[0].strip(), parts[1].strip())
            if key in out:
                errs.append("%s: duplicate entry for %s %s" % (where, key[0], key[1]))
                continue
            out[key] = {"reason": parts[2], "approval": parts[3], "doc": parts[4], "sites": sites, "line": n}
    return out, errs


def collect():
    """every banned site in the library and test dirs, test code included: (rel, line, fn, token)"""
    hits = []
    for rel, p in rs_files():
        lines = open(p, errors="replace").read().split("\n")
        cur = "?"
        for i, l in enumerate(lines):
            m = re.search(r"\bfn\s+(\w+)", l)
            if m:
                cur = m.group(1)
            code = l.split("//")[0]
            for name, rx in TOKENS:
                for _ in rx.finditer(code):
                    hits.append((rel, i + 1, cur, name))
            dm = DERIVE.search(code)
            if dm:
                items = [x.strip() for x in dm.group(1).split(",")]
                if "Clone" in items and "Copy" not in items:
                    hits.append((rel, i + 1, cur, "derive(Clone)"))
    return hits


def scan():
    """-> dict(hits, bad, allowed, stale, errors, entries). allowed = excepted `.clone()` sites."""
    entries, errors = load_allow()
    hits = collect()
    used = collections.Counter()
    bad, allowed = [], []
    for h in hits:
        key = (h[0], h[2])
        if h[3] == "clone()" and key in entries and used[key] < entries[key]["sites"]:
            used[key] += 1
            allowed.append(h)
        else:
            bad.append(h)
    seen = {(h[0], h[2]) for h in hits if h[3] == "clone()"}
    stale = [k for k in entries if k not in seen]
    return {"hits": hits, "bad": bad, "allowed": allowed, "stale": stale, "errors": errors, "entries": entries}


def helper_sites():
    """call/use sites of the hidden copy helpers, definitions excluded -> [(rel, line, name)].
    Counted: every identifier occurrence of a HELPERS name in code (comments and string literals ignored) after
    removing `fn NAME` / `macro_rules! NAME` definitions and whole `use` lines. Plus Rc::clone / Arc::clone in
    zenith-float-compare/src. Test code (cfg(test), tests.rs) is included."""
    out = []
    for rel, p in rs_files():
        lines = open(p, errors="replace").read().split("\n")
        in_compare = rel.startswith(HANDLE_DIR + "/")
        for i, l in enumerate(lines):
            code = strip_code(l)
            if USE_LINE.match(code):
                continue
            code = HELPER_DEF.sub("", code)
            for m in HELPER_USE.finditer(code):
                out.append((rel, i + 1, m.group(1)))
            if in_compare:
                for m in HANDLE_COPY.finditer(code):
                    out.append((rel, i + 1, m.group(1) + "::clone"))
    return out


def helper_copy_count():
    return len(helper_sites())


def main():
    r = scan()
    bad, allowed = r["bad"], r["allowed"]
    by = collections.Counter(h[3] for h in bad)
    per = collections.Counter(h[0] for h in bad)
    if "--helpers" in sys.argv:
        sites = helper_sites()
        names = collections.Counter(s[2] for s in sites)
        if "--list" in sys.argv:
            for s in sites:
                print("%s:%d %s" % s)
        print("helper_copies: %d %s" % (len(sites), dict(sorted(names.items()))))
        return 0
    if "--list" in sys.argv:
        for h in bad:
            print("%s:%d [%s] %s" % h)
        for h in allowed:
            print("%s:%d [%s] %s (excepted)" % h)
    for e in r["errors"]:
        print("no_clone: ALLOWLIST ERROR: " + e)
    for k in r["stale"]:
        print("no_clone: STALE allowlist entry %s %s: no .clone() site matches it any more; remove it" % k)
    print("no_clone: %d violation(s) %s; %d excepted site(s); %d stale entr(y/ies)" % (len(bad), dict(by), len(allowed), len(r["stale"])))
    if "--files" in sys.argv:
        for f, n in per.most_common():
            print("%5d %s" % (n, f))
    return 1 if (bad or r["stale"] or r["errors"]) else 0


if __name__ == "__main__":
    sys.exit(main())
