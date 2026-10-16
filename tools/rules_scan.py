#!/usr/bin/env python3
"""rules_scan: scan the ADDED lines of a diff for banned constructs in library code.

Run on EVERY model-written diff before review, and in tools/../full_gate (step `rules-scan`).

Usage:
  tools/rules_scan.sh                      # origin/main..HEAD, plus the clone-count baseline check
  tools/rules_scan.sh A..B                 # any range (a bare rev means REV..HEAD)
  tools/rules_scan.sh --staged             # git diff --cached
  tools/rules_scan.sh --working            # git diff HEAD (staged + unstaged)
  tools/rules_scan.sh --diff-file F        # a saved unified diff (used by the self-test)
  tools/rules_scan.sh --list-test-panics   # every unwrap/expect/panic!/... site in test code (file:line [kind]), then the count
  tools/rules_scan.sh --count-only         # only the count baseline checks (clone, helper-copy, ref, neg_via_mul ratchets) and the exception checks
  tools/rules_scan.sh --update-baseline    # lower the baselines to the current counts (never raises one)

LIBRARY PATHS (clone, to_vec, cloned, to_owned, unsafe, allow, unwrap, expect, panic): .rs files under
  zenith-float-num/src/  zenith-float-macro/src/  zenith-float-compare/src/  src/
  excluding any path with a tools/ component. Test code (an added line inside a #[cfg(test)] item of the file at the
  diff's tip, files named tests.rs, tests/ benches/ examples/ components) is no longer excluded: the clone rules
  (clone, to_vec, cloned) apply to it. It stays exempt from the other rules (unsafe, allow, unwrap, expect, panic,
  f32, f64), whose meaning is unchanged. The counters below also walk the test dirs, see no_clone.TEST_DIRS.
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
count check: a stale or malformed entry fails, and the number of excepted `.clone()` sites, plus 1 for the
approved Ziv working-copy helper `ziv_pass_operand` (the one fn of that name, in zenith-float-num/src/num/ziv_operand.rs, with
a CODEOWNERS entry; any second helper of that shape fails), may never exceed tools/exception_count_baseline.txt (now 1;
lowered with --update-baseline, never raised).

REF-SUM COUNT: total uses of `add_ref(` / `sub_ref(` (the borrowed-borrowed sums, sanctioned only where both
operands are still live) in non-test library code (REF_SCOPE = "library": scoped on purpose, test uses are not
counted; the report line also prints the with-tests figure), definitions (`fn add_ref`) excluded. Same ratchet as the
clone count against tools/ref_sum_count_baseline.txt: it may only go down; --update-baseline lowers it.

MUL-REF COUNT: the same ratchet, separately, for `mul_ref(` / `div_ref(` (borrowed-borrowed products and
quotients, sanctioned only where both operands are still live) against tools/mul_ref_count_baseline.txt.

NEG-VIA-MUL COUNT: total uses of `neg_via_mul(` (negation through a fresh `-1` and the general multiplication
kernel, for an operand that has to stay borrowed) in ALL scanned code: the four library dirs and no_clone.TEST_DIRS
(src, tests, benches, examples, nostd-tests), test code included, because a test use is still a use that shows the
idiom spreading. NOT counted: the definitions (`fn neg_via_mul`, ExactNum in ext.rs and ExactComplex in complex.rs, the
approved exception), comments, doc comments and string literals, and names that merely start with the word (a test
called `neg_via_mul_equals_neg_...` is not a use). The count only exists to stop growth: it may only go down, against
tools/neg_via_mul_count_baseline.txt (CODEOWNERS: Jeff), lowered with --update-baseline. Prefer `neg` on an owned
operand; a use that survives must be a borrowed operand that is still live.
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
    ("unwrap", re.compile(r"\bunwrap(?:_err)?\(")),
    ("expect", re.compile(r"\bexpect(?:_err)?\(")),
    ("panic", re.compile(r"\bpanic!")),
    ("unreachable", re.compile(r"\bunreachable!")),
    ("todo", re.compile(r"\btodo!")),
    ("unimplemented", re.compile(r"\bunimplemented!")),
    ("f32", re.compile(r"\bf32\b|\d_?f32\b")),
    ("f64", re.compile(r"\bf64\b|\d_?f64\b")),
]
CLONE_RULES = ("clone", "to_vec", "cloned")  # applied to test code too
# Jeff: no crash in TEST code either. These apply to test code as well (counted by the test-panic ratchet, and flagged
# on added test lines by the --staged / --working / --diff-file scans). assert!/assert_eq!/debug_assert! are fine.
PANIC_RULES = ("unwrap", "expect", "panic", "unreachable", "todo", "unimplemented")
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


def is_test_file(rel):
    """whole-file test code: the test dirs (no_clone.TEST_DIRS), any tests.rs, any tests/ benches/ examples/ component"""
    parts = rel.split("/")
    if any(rel == d or rel.startswith(d + "/") for d in no_clone.TEST_DIRS):
        return True
    if rel.startswith("zenith-float-compare/src/"):  # benchmark adapters: a crash there is a crash in a test run
        return True
    return parts[-1] == "tests.rs" or bool(REF_EXCL_PARTS & set(parts[:-1]))


def is_test_dir_file(path):
    """.rs files in no_clone.TEST_DIRS (outside the library dirs): scanned for the clone and PANIC rules only"""
    return path.endswith(".rs") and "tools" not in path.split("/")[:-1] and \
        any(path.startswith(d + "/") for d in no_clone.TEST_DIRS)


def test_panic_sites(with_lib=False):
    """-> [(rel, line, label)] of the PANIC_RULES constructs in TEST code (cfg(test) items and whole test files, the
    compare crate counts as test code) of every scanned file (library + test dirs). with_lib=True also returns the same
    constructs in the non-test library code (which must be 0: the diff scan flags them as library violations)."""
    out = []
    for rel, p in no_clone.rs_files():
        lines = open(p, errors="replace").read().split("\n")
        whole = is_test_file(rel)
        mask = no_clone.test_mask(lines)
        for i, l in enumerate(lines):
            if not (whole or mask[i]) and not with_lib:
                continue
            code = strip(l)
            for label, rx in PATTERNS:
                if label in PANIC_RULES:
                    for _ in rx.finditer(code):
                        out.append((rel, i + 1, label))
    return out


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
        if not (is_lib(path) or is_test_dir_file(path)):
            continue
        if path not in masks:
            t = tip_text(path, mode, rev)
            masks[path] = no_clone.test_mask(t.split("\n")) if t else None
        m = masks[path]
        # test code (#[cfg(test)] items, tests.rs, tests/ benches/ examples/ files) is exempt from the NON-clone rules
        # (unwrap, expect, panic!, unsafe, allow, f32/f64) as before, but the clone rules apply to it: it is counted.
        in_test = (m is not None and 0 < n <= len(m) and m[n - 1]) or not is_lib_only(path)
        # PANIC_RULES (unwrap, expect, panic!, unreachable!, todo!, unimplemented!) apply to test code too (Jeff). In a
        # RANGE scan the legacy test sites make the added lines of old commits fail, so there the test-panic RATCHET
        # (tools/test_panic_count_baseline.txt) is the enforcement; --staged/--working/--diff-file flag them per line.
        strict_panic = mode != "range"
        code = strip(content)
        key = "%s:%s" % (path, hashlib.sha1(" ".join(content.split()).encode()).hexdigest()[:12])
        for label, rx in PATTERNS:
            if label in ("f32", "f64") and not CALC.match(path):
                continue
            if in_test and label not in CLONE_RULES and not (strict_panic and label in PANIC_RULES):
                continue
            if label not in CLONE_RULES + PANIC_RULES and not is_lib(path):
                continue  # test-dir files outside the library dirs: clone and PANIC rules only
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


def neg_via_mul_count():
    """uses of `neg_via_mul(` in all scanned code (library + test dirs); definitions, comments and strings excluded."""
    use = re.compile(r"\bneg_via_mul\(")
    deff = re.compile(r"\bfn\s+neg_via_mul\b")
    total = 0
    for rel, p in no_clone.rs_files(None):
        for l in open(p, errors="replace").read().split("\n"):
            code = strip(l)
            if not deff.search(code):
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
    rc, ziv = ziv_helper_check()
    for e in r["errors"]:
        print("rules_scan: clone allowlist error: " + e)
        rc = 1
    for k in r["stale"]:
        print("rules_scan: STALE clone allowlist entry %s %s: no .clone() site matches it any more; remove it" % k)
        rc = 1
    for h in r["allowed"]:
        print("rules_scan: excepted site %s:%d [%s]" % (h[0], h[1], h[2]))
    return max(rc, ratchet("clone exceptions", "exception_count_baseline.txt", len(r["allowed"]) + ziv, update, "site count"))


# The one approved Ziv working-copy helper (Jeff, via Paul: one helper, one file, CODEOWNERS-protected, a non-clone-shaped name,
# one allocation plus one limb copy per retry pass at that pass's precision). Its definition counts 1 in the exception total
# (tools/exception_count_baseline.txt, 0 -> 1, only down). A second helper of that kind, a second definition, a definition
# in another file or a missing CODEOWNERS entry fails regardless of any baseline.
ZIV_NAME = no_clone.ZIV_HELPER
ZIV_FILE = no_clone.ZIV_HELPER_FILE
_ZIV_WORDS = r"(?:pass|operand|copy|clone|dup|widen|working|extend)"
ZIV_DEF = re.compile(r"\bfn\s+(\w*ziv\w*%s\w*|\w*%s\w*ziv\w*)\b" % (_ZIV_WORDS, _ZIV_WORDS))


def ziv_helper_defs():
    """[(rel, line, name)] of every fn that looks like a Ziv working-copy helper (a name with `ziv` and pass/operand/copy/
    clone/dup/widen/working/extend), in all scanned code, tests included."""
    out = []
    for rel, p in no_clone.rs_files(None):
        for i, l in enumerate(open(p, errors="replace").read().split("\n")):
            for m in ZIV_DEF.finditer(strip(l)):
                out.append((rel, i + 1, m.group(1)))
    return out


def ziv_helper_check():
    """-> (rc, count): count is the helper's contribution to the exception total (1 when the approved definition exists)"""
    defs = ziv_helper_defs()
    rc = 0
    good = [d for d in defs if d[2] == ZIV_NAME and d[0] == ZIV_FILE]
    seen_good = False
    for d in defs:
        if d in good and not seen_good:
            seen_good = True
            continue
        print("rules_scan: second Ziv working-copy helper %s:%d fn %s (only one `fn %s`, in %s, is approved)" % (d + (ZIV_NAME, ZIV_FILE)))
        rc = 1
    if any(d[2] == ZIV_NAME and d[0] != ZIV_FILE for d in defs) or (defs and not good and any(d[2] == ZIV_NAME for d in defs)):
        print("rules_scan: %s must live in %s" % (ZIV_NAME, ZIV_FILE))
        rc = 1
    if good:
        try:
            owners = open(os.path.join(ROOT, ".github", "CODEOWNERS"), errors="replace").read().split("\n")
        except OSError:
            owners = []
        if not any(l.split()[:1] == ["/" + ZIV_FILE] and len(l.split()) > 1 for l in owners if not l.lstrip().startswith("#")):
            print("rules_scan: %s has no CODEOWNERS entry (/%s @owner)" % (ZIV_FILE, ZIV_FILE))
            rc = 1
        print("rules_scan: Ziv working-copy helper %s: 1 definition in %s" % (ZIV_NAME, ZIV_FILE))
    return rc, len(defs)


def test_panic_check(update=False):
    """unwrap/expect/panic!/unreachable!/todo!/unimplemented! in test code: may only go down"""
    sites = test_panic_sites()
    by = collections.Counter(x[2] for x in sites)
    return ratchet("test-panic", "test_panic_count_baseline.txt", len(sites), update, "site count", " %s" % dict(sorted(by.items())))


def neg_via_mul_check(update=False):
    return ratchet("neg_via_mul", "neg_via_mul_count_baseline.txt", neg_via_mul_count(), update, "use count")


def both_count_checks(update=False):
    rc = max(count_check(update), ref_count_check(update), helper_count_check(update), exception_check(update),
             test_panic_check(update), neg_via_mul_check(update))
    if rc == 0:
        print("rules_scan: all count baselines hold")
    return rc


def main(argv):
    if "--list-test-panics" in argv:
        sites = test_panic_sites()
        for rel, n, label in sites:
            print("%s:%d [%s]" % (rel, n, label))
        print("rules_scan: %d test-panic sites %s" % (len(sites), dict(sorted(collections.Counter(x[2] for x in sites).items()))))
        return 0
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
