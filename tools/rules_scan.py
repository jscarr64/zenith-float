#!/usr/bin/env python3
"""rules_scan: scan the ADDED lines of a diff for banned constructs in library code.

Run on EVERY model-written diff before review, and in tools/../full_gate (step `rules-scan`).

Usage:
  tools/rules_scan.sh                      # origin/main..HEAD, plus the clone-count baseline check
  tools/rules_scan.sh A..B                 # any range (a bare rev means REV..HEAD)
  tools/rules_scan.sh --staged             # git diff --cached
  tools/rules_scan.sh --working            # git diff HEAD (staged + unstaged)
  tools/rules_scan.sh --diff-file F        # a saved unified diff (used by the self-test)
  tools/rules_scan.sh --count-only         # only the count baseline checks (clone, helper-copy, ref ratchets) and the exception checks
  tools/rules_scan.sh --update-baseline    # lower the baselines to the current counts (never raises one)

LIBRARY PATHS (clone, to_vec, cloned, to_owned, unsafe, allow, unwrap, expect, panic): .rs files under
  zenith-float-num/src/  zenith-float-macro/src/  zenith-float-compare/src/  src/
  excluding any path with a tools/ component, and any added line inside a #[cfg(test)] item of the file at the
  diff's tip. (tests/ benches/ examples/ components and files named tests.rs are NOT excluded any more: test code is
  counted. The counters below also walk the test dirs, see no_clone.TEST_DIRS.)
CALCULATION PATHS (f32, f64): the library paths above except zenith-float-compare/src/ (benchmark adapters).

Flagged on added lines only (comments and string literals are ignored):
  .clone(  .to_vec(  .cloned(  .to_owned(  unsafe  #[allow / #![allow  unwrap(  expect(  panic!  f32  f64

ALLOWLIST: tools/rules_allowlist.txt, one `path:hash<TAB>reason` per line; the key is printed as allow-key=
for every flagged site (hash = sha1 of the added line with whitespace collapsed). The file starts EMPTY.
ONLY JEFF MAY EXTEND IT (.github/CODEOWNERS requires his review). Models never add entries.

CLONE COUNT: total `.clone(` occurrences in ALL scanned code, test code included (the four library dirs and
no_clone.TEST_DIRS; tools/ and target dirs excluded) must never go up versus tools/clone_count_baseline.txt. The
baseline is only ever lowered, with --update-baseline.

HELPER-COPY COUNT: total uses of the hidden copy helpers (no_clone.HELPERS: exact_copy_ext, exact_copy, exact_operand,
copy_p, svd_copy_prec, lift_num) in non-definition positions, plus Rc::clone / Arc::clone in zenith-float-compare/src
(see no_clone.helper_sites for exactly how they are counted). Same ratchet against tools/helper_copy_count_baseline.txt.

EXCEPTIONS: tools/clone_allowlist.txt (see tools/no_clone.py and tools/NO_CLONE_RULES.md). Checked here on every
count check: a stale or malformed entry fails, and the number of excepted `.clone()` sites may never exceed
tools/exception_count_baseline.txt (lowered with --update-baseline, never raised).

REF-SUM COUNT: total uses of `add_ref(` / `sub_ref(` (the borrowed-borrowed sums, sanctioned only where both
operands are still live) in non-test library code (REF_SCOPE = "library": scoped on purpose, test uses are not
counted; the report line also prints the with-tests figure), definitions (`fn add_ref`) excluded. Same ratchet as the
clone count against tools/ref_sum_count_baseline.txt: it may only go down; --update-baseline lowers it.

MUL-REF COUNT: the same ratchet, separately, for `mul_ref(` / `div_ref(` (borrowed-borrowed products and
quotients, sanctioned only where both operands are still live) against tools/mul_ref_count_baseline.txt.
"""
import hashlib, os, re, subprocess, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("RULES_SCAN_ROOT") or os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ["NO_CLONE_ROOT"] = ROOT
import no_clone  # test_mask, DIRS

LIB = re.compile(r"^(zenith-float-num/src|zenith-float-macro/src|zenith-float-compare/src|src)/.*\.rs$")
CALC = re.compile(r"^(zenith-float-num/src|zenith-float-macro/src|src)/.*\.rs$")
EXCL_PARTS = {"tools"}
# Scope of the add_ref/sub_ref/mul_ref/div_ref ratchets: the library only (the pre-test-counting meaning).
REF_SCOPE = "library"
REF_EXCL_PARTS = {"tests", "benches", "examples", "tools"}
PATTERNS = [
    ("clone", re.compile(r"\.clone\(")),
    ("to_vec", re.compile(r"\.to_vec\(")),
    ("cloned", re.compile(r"\.cloned\(")),
    ("to_owned", re.compile(r"\.to_owned\(")),
    ("unsafe", re.compile(r"\bunsafe\b")),
    ("allow", re.compile(r"#!?\[allow\b|^\s*#!?\[.*\ballow\(")),
    ("unwrap", re.compile(r"\bunwrap\(")),
    ("expect", re.compile(r"\bexpect\(")),
    ("panic", re.compile(r"\bpanic!")),
    ("f32", re.compile(r"\bf32\b|\d_?f32\b")),
    ("f64", re.compile(r"\bf64\b|\d_?f64\b")),
]
STR = re.compile(r'"(?:\\.|[^"\\])*"')


def is_lib(path):
    """library-dir .rs files, test files included (only tools/ is excluded)"""
    parts = path.split("/")
    if not LIB.match(path):
        return False
    return not (EXCL_PARTS & set(parts[:-1]))


def is_lib_only(path):
    """the old library-only scope (no tests/ benches/ examples/ components, no tests.rs): used by the ref ratchets"""
    parts = path.split("/")
    if not LIB.match(path) or parts[-1] == "tests.rs":
        return False
    return not (REF_EXCL_PARTS & set(parts[:-1]))


def strip(line):
    s = STR.sub('""', line)
    s = re.sub(r"/\*.*?\*/", "", s)
    return s.split("//")[0]


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True)


def tip_text(path, mode, rev):
    if mode == "file":
        return None
    if mode == "working":
        p = os.path.join(ROOT, path)
        return open(p, errors="replace").read() if os.path.exists(p) else None
    spec = ":" + path if mode == "staged" else "%s:%s" % (rev, path)
    r = git("show", spec)
    return r.stdout if r.returncode == 0 else None


def parse_diff(text):
    """yield (path, newline_no, content) for each added line"""
    path, n = None, 0
    for raw in text.split("\n"):
        if raw.startswith("+++ "):
            path = raw[4:].strip()
            path = path[2:] if path.startswith("b/") else (None if path == "/dev/null" else path)
        elif raw.startswith("--- ") or raw.startswith("diff --git"):
            continue
        elif raw.startswith("@@"):
            m = re.search(r"\+(\d+)", raw)
            n = int(m.group(1)) if m else 0
        elif path and raw.startswith("+"):
            yield path, n, raw[1:]
            n += 1
        elif path and raw.startswith(" "):
            n += 1


def load_allow():
    out = set()
    p = os.path.join(ROOT, "tools", "rules_allowlist.txt")
    if os.path.exists(p):
        for l in open(p):
            l = l.rstrip("\n")
            if l.strip() and not l.lstrip().startswith("#"):
                out.add(l.split("\t")[0].strip())
    return out


def scan_diff(text, mode, rev):
    allow, hits, masks = load_allow(), [], {}
    for path, n, content in parse_diff(text):
        if not is_lib(path):
            continue
        if path not in masks:
            t = tip_text(path, mode, rev)
            masks[path] = no_clone.test_mask(t.split("\n")) if t else None
        m = masks[path]
        if m is not None and 0 < n <= len(m) and m[n - 1]:
            continue
        code = strip(content)
        key = "%s:%s" % (path, hashlib.sha1(" ".join(content.split()).encode()).hexdigest()[:12])
        for label, rx in PATTERNS:
            if label in ("f32", "f64") and not CALC.match(path):
                continue
            if rx.search(code) and key not in allow:
                hits.append((path, n, label, " ".join(content.split()), key))
    return hits


def clone_count():
    """`.clone(` occurrences in every scanned file (library + test dirs), test code included"""
    total = 0
    for rel, p in no_clone.rs_files():
        for l in open(p, errors="replace").read().split("\n"):
            total += len(re.findall(r"\.clone\(", strip(l)))
    return total


def helper_copy_count():
    return no_clone.helper_copy_count()


REF_KINDS = (
    # (label, function names, baseline file)
    ("add_ref/sub_ref", ("add_ref", "sub_ref"), "ref_sum_count_baseline.txt"),
    ("mul_ref/div_ref", ("mul_ref", "div_ref"), "mul_ref_count_baseline.txt"),
)


def ref_sum_count(names=("add_ref", "sub_ref"), with_tests=False):
    """uses of the named ref functions, definitions (`fn NAME`), comments and strings excluded.
    Default scope is the library only (REF_SCOPE): #[cfg(test)] items, tests.rs and tests/benches/examples components
    are skipped. with_tests=True counts every scanned file instead (library + test dirs); it is for the report only."""
    use = re.compile(r"\b(?:%s)\(" % "|".join(names))
    deff = re.compile(r"\bfn\s+(?:%s)\b" % "|".join(names))
    total = 0
    for rel, p in no_clone.rs_files(None if with_tests else no_clone.DIRS):
        if not with_tests and not is_lib_only(rel):
            continue
        lines = open(p, errors="replace").read().split("\n")
        mask = [False] * len(lines) if with_tests else no_clone.test_mask(lines)
        for i, l in enumerate(lines):
            code = strip(l)
            if not mask[i] and not deff.search(code):
                total += len(use.findall(code))
    return total


def baseline_path():
    return os.path.join(ROOT, "tools", "clone_count_baseline.txt")


def ref_baseline_path(name="ref_sum_count_baseline.txt"):
    return os.path.join(ROOT, "tools", name)


def read_baseline():
    try:
        return int(open(baseline_path()).read().split()[0])
    except (OSError, ValueError, IndexError):
        return None


def count_check(update=False):
    cur, base = clone_count(), read_baseline()
    if base is None:
        print("rules_scan: missing or unreadable tools/clone_count_baseline.txt")
        return 1
    if cur > base:
        print("rules_scan: clone count %d is ABOVE the baseline %d (it may only go down)" % (cur, base))
        return 1
    if cur < base:
        if update:
            open(baseline_path(), "w").write("%d\n" % cur)
            print("rules_scan: baseline lowered %d -> %d" % (base, cur))
        else:
            print("rules_scan: clone count %d is below baseline %d: lower the baseline (rules_scan.sh --update-baseline)" % (cur, base))
        return 0
    print("rules_scan: clone count %d == baseline" % cur)
    return 0


def ref_count_check(update=False):
    rc = 0
    for label, names, fname in REF_KINDS:
        cur = ref_sum_count(names)
        path = ref_baseline_path(fname)
        try:
            base = int(open(path).read().split()[0])
        except (OSError, ValueError, IndexError):
            print("rules_scan: missing or unreadable tools/%s" % fname)
            rc = 1
            continue
        if cur > base:
            print("rules_scan: %s use count %d is ABOVE the baseline %d (it may only go down)" % (label, cur, base))
            rc = 1
        elif cur < base:
            if update:
                open(path, "w").write("%d\n" % cur)
                print("rules_scan: %s baseline lowered %d -> %d" % (label, base, cur))
            else:
                print("rules_scan: %s use count %d is below baseline %d: lower the baseline (rules_scan.sh --update-baseline)" % (label, cur, base))
        else:
            print("rules_scan: %s use count %d == baseline" % (label, cur))
        wt = ref_sum_count(names, with_tests=True)
        if wt != cur:
            print("rules_scan: note: %s uses incl. test code = %d (ratchet scope: %s only, %d test uses not counted)" % (label, wt, REF_SCOPE, wt - cur))
    return rc


def ratchet(label, fname, cur, update, what="count", extra=""):
    """cur vs tools/<fname>: above -> fail; below -> hint (or lower with update); equal -> ok. Returns rc."""
    path = os.path.join(ROOT, "tools", fname)
    try:
        base = int(open(path).read().split()[0])
    except (OSError, ValueError, IndexError):
        print("rules_scan: missing or unreadable tools/%s" % fname)
        return 1
    if cur > base:
        print("rules_scan: %s %s %d is ABOVE the baseline %d (it may only go down)%s" % (label, what, cur, base, extra))
        return 1
    if cur < base:
        if update:
            open(path, "w").write("%d\n" % cur)
            print("rules_scan: %s baseline lowered %d -> %d" % (label, base, cur))
        else:
            print("rules_scan: %s %s %d is below baseline %d: lower the baseline (rules_scan.sh --update-baseline)" % (label, what, cur, base))
        return 0
    print("rules_scan: %s %s %d == baseline%s" % (label, what, cur, extra))
    return 0


def helper_count_check(update=False):
    return ratchet("helper-copy", "helper_copy_count_baseline.txt", helper_copy_count(), update)


def exception_check(update=False):
    """clone allowlist: malformed or stale entries fail; the excepted-site count may only go down"""
    r = no_clone.scan()
    rc = 0
    for e in r["errors"]:
        print("rules_scan: clone allowlist error: " + e)
        rc = 1
    for k in r["stale"]:
        print("rules_scan: STALE clone allowlist entry %s %s: no .clone() site matches it any more; remove it" % k)
        rc = 1
    for h in r["allowed"]:
        print("rules_scan: excepted site %s:%d [%s]" % (h[0], h[1], h[2]))
    return max(rc, ratchet("clone exceptions", "exception_count_baseline.txt", len(r["allowed"]), update, "site count"))


def both_count_checks(update=False):
    rc = max(count_check(update), ref_count_check(update), helper_count_check(update), exception_check(update))
    if rc == 0:
        print("rules_scan: all count baselines hold")
    return rc


def main(argv):
    if "--count-only" in argv:
        return both_count_checks()
    if "--update-baseline" in argv:
        return both_count_checks(update=True)
    mode, rev, with_count = "range", "HEAD", False
    if "--diff-file" in argv:
        mode = "file"
        text = open(argv[argv.index("--diff-file") + 1], errors="replace").read()
    elif "--staged" in argv:
        mode = "staged"
        text = git("diff", "--cached", "--no-color", "--no-ext-diff", "-U0").stdout
    elif "--working" in argv:
        mode = "working"
        text = git("diff", "HEAD", "--no-color", "--no-ext-diff", "-U0").stdout
    else:
        args = [a for a in argv if not a.startswith("--")]
        rng = args[0] if args else "origin/main..HEAD"
        if ".." not in rng:
            rng += "..HEAD"
        rev = rng.split("..")[-1].lstrip(".") or "HEAD"
        r = git("diff", rng, "--no-color", "--no-ext-diff", "-U0")
        if r.returncode != 0:
            print("rules_scan: git diff %s failed: %s" % (rng, r.stderr.strip()))
            return 2
        text, with_count = r.stdout, True
    hits = scan_diff(text, mode, rev)
    for path, n, label, content, key in hits:
        print("%s:%d [%s] %s  allow-key=%s" % (path, n, label, content[:120], key))
    by = collections.Counter(h[2] for h in hits)
    print("rules_scan: %d flagged %s" % (len(hits), dict(by)))
    rc = 1 if hits else 0
    if with_count and not os.environ.get("RULES_SCAN_SKIP_COUNT"):
        rc = max(rc, both_count_checks())
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
