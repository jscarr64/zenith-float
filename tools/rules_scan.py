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

BORROW RULE: ExactNum::{add_ref, sub_ref, mul_ref, div_ref} (zenith-float-num/src/ext.rs) are the sanctioned
borrowing ops: both operands read in place, only the result allocated. The rule fails closed, by mechanism, not by
spelling. REACH: every non-test fn and macro_rules! of the whole zenith-float-num/src reachable from them, through
every call form (method, UFCS / `<T as Tr>::f`, paths, turbofish and generic calls, fns used as values, operators into
op-trait impls, macro invocations). A call is resolved by impl type, field / local types and arity; one that cannot be
resolved (unknown or generic receiver, trait or unknown path, a type without its own such fn) follows EVERY crate fn
of that name; only inherent fns of std types (Vec::new) are not followed. COPY: any mechanism that duplicates data:
the clone / copy family as methods, UFCS paths or fn values (clone, clone_from, clone_into, clone_from_slice,
copy_from_slice, swap_with_slice, cloned, copied, to_owned, into_owned, to_vec*, extend_from_slice/within); slice ->
owned conversions (`Vec|Box|Rc|...::from`, `from_iter`, `.as_slice().into()`, `x[a..b].into()`); a closure or pattern
that binds an item by value or derefs it (`|&x|`, `|x| *x`, typed or braced, `for &x in`, `let Some(&x)`); `extend`
fed by anything but an owned iterator; `vec![e; n]` with a non-literal e and `[e; n]` with operand data; helper-copy
calls. Comments, strings and char literals are stripped first. A site passes only as a NAMED exception
(BORROW_EXCEPTIONS: file, fn, token, max count, reason); the scan prints the site count and the exception count and
fails on any unlisted site, naming the op and the reach chain: `BORROW file:line in fn f: `tok` (add_ref -> g -> f)`.
A root that is not found fails; --borrow-reach lists the reached fns. (This replaces the former add_ref/sub_ref and mul_ref/div_ref use-count ratchets: a use of a borrowing op copies
nothing, so it needs no count.) Skipped when zenith-float-num/src/ext.rs does not exist.

NEG-VIA-MUL COUNT: total uses of `neg_via_mul(` (negation through a fresh `-1` and the general multiplication
kernel, for an operand that has to stay borrowed) in ALL scanned code: the four library dirs and no_clone.TEST_DIRS
(src, tests, benches, examples, nostd-tests), test code included, because a test use is still a use that shows the
idiom spreading. NOT counted: the definitions (`fn neg_via_mul`, ExactNum in ext.rs and ExactComplex in complex.rs, the
approved exception), comments, doc comments and string literals, and names that merely start with the word (a test
called `neg_via_mul_equals_neg_...` is not a use). The count only exists to stop growth: it may only go down, against
tools/neg_via_mul_count_baseline.txt (CODEOWNERS: Jeff), lowered with --update-baseline. Prefer `neg` on an owned
operand; a use that survives must be a borrowed operand that is still live.

HELPER-COPY RULING (Jeff, 2026-10-16, tools/NO_CLONE_RULES.md "Ruling: helper_copy"): (1) COMBINED: clone count +
helper-copy count may never exceed tools/clone_helper_combined_baseline.txt (only down). (2) PER COMMIT (range mode):
in every commit of the range after tools/helper_copy_rule_base.txt (every commit when that file is absent), the
helper-copy count may rise only if the clone count falls by at least as much in that same commit. (3) PINNED SITES:
every use of `ziv_pass_operand` outside its own file must sit in a (file, fn) listed in tools/ziv_shape_a_sites.txt
(`path<TAB>fn<TAB>reason`), and that list must be covered by a .github/CODEOWNERS line owned by @jscarr64. (1) and
the CODEOWNERS check apply whenever tools/ziv_shape_a_sites.txt exists; (3) applies to every use.
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


# ---- helper-copy ruling (Jeff, 2026-10-16) ----
PINNED_FILE = "tools/ziv_shape_a_sites.txt"
COMBINED_FILE = "clone_helper_combined_baseline.txt"
RULE_BASE_FILE = "tools/helper_copy_rule_base.txt"
FN_DEF = re.compile(r"\bfn\s+(\w+)")


def pinned_sites():
    """{(path, fn)} from tools/ziv_shape_a_sites.txt, or None when the file is absent"""
    try:
        text = open(os.path.join(ROOT, PINNED_FILE), errors="replace").read()
    except OSError:
        return None
    out = set()
    for l in text.split("\n"):
        if not l.strip() or l.lstrip().startswith("#"):
            continue
        f = l.split("\t")
        if len(f) >= 2:
            out.add((f[0].strip(), f[1].strip()))
    return out


def enclosing_fn(rel, line):
    lines = open(os.path.join(ROOT, rel), errors="replace").read().split("\n")
    for l in reversed(lines[:line]):
        m = FN_DEF.search(strip(l))
        if m:
            return m.group(1)
    return None


def codeowners_cover(path, owner="@jscarr64"):
    try:
        owners = open(os.path.join(ROOT, ".github", "CODEOWNERS"), errors="replace").read().split("\n")
    except OSError:
        return False
    return any(l.split()[:1] == ["/" + path] and owner in l.split()[1:] for l in owners if not l.lstrip().startswith("#"))


def ruling_check(update=False):
    rc = 0
    pinned = pinned_sites()
    for rel, line, name in no_clone.helper_sites():
        if name != ZIV_NAME:
            continue
        fn = enclosing_fn(rel, line)
        if pinned is None or (rel, fn) not in pinned:
            print("rules_scan: ziv_pass_operand use at %s:%d (fn %s) is not on the pinned shape-A list %s" % (rel, line, fn, PINNED_FILE))
            rc = 1
    if pinned is not None:
        if not codeowners_cover(PINNED_FILE):
            print("rules_scan: %s is not covered by CODEOWNERS for @jscarr64 (/%s @jscarr64)" % (PINNED_FILE, PINNED_FILE))
            rc = 1
        rc = max(rc, ratchet("clone+helper_copy combined", COMBINED_FILE, clone_count() + helper_copy_count(), update))
    return rc


def counts_at(rev, tmp):
    """(clone, helper_copy) of the tree at rev, counted by this scanner"""
    d = os.path.join(tmp, rev)
    os.makedirs(d, exist_ok=True)
    a = subprocess.run(["git", "-C", ROOT, "archive", rev], capture_output=True)
    if a.returncode != 0:
        return None
    t = subprocess.run(["tar", "-x", "-C", d], input=a.stdout, capture_output=True)
    if t.returncode != 0:
        return None
    env = dict(os.environ, RULES_SCAN_ROOT=d, NO_CLONE_ROOT=d)
    r = subprocess.run([sys.executable, os.path.abspath(__file__), "--print-counts"], env=env, capture_output=True, text=True)
    f = r.stdout.split()
    if r.returncode != 0 or len(f) != 2:
        return None
    return int(f[0]), int(f[1])


EXCEPTIONS_FILE = "tools/helper_copy_commit_exceptions.txt"


def helper_copy_exceptions():
    """{full 40-char SHA: allowed helper_copy rise} from `<sha> <count> <reason>` lines (Rusty rulings only)"""
    try:
        lines = open(os.path.join(ROOT, EXCEPTIONS_FILE)).read().splitlines()
    except OSError:
        return {}
    out = {}
    for l in lines:
        f = l.split()
        if l.startswith("#") or len(f) < 3 or len(f[0]) != 40 or not f[1].isdigit():
            continue
        out[f[0]] = int(f[1])
    return out


def per_commit_check(rng):
    import tempfile
    commits = [c for c in git("rev-list", "--reverse", "--first-parent", rng).stdout.split()]
    try:
        base = open(os.path.join(ROOT, RULE_BASE_FILE)).read().split()[0]
        after = set(git("rev-list", base + ".." + rng.split("..")[-1]).stdout.split())
        commits = [c for c in commits if c in after]
    except (OSError, IndexError):
        pass
    rc = 0
    with tempfile.TemporaryDirectory() as tmp:
        cache = {}
        for c in commits:
            parent = git("rev-parse", "--verify", "-q", c + "^").stdout.strip()
            if not parent:
                continue
            for r in (parent, c):
                if r not in cache:
                    cache[r] = counts_at(r, tmp)
            a, b = cache[parent], cache[c]
            if a is None or b is None:
                print("rules_scan: could not count commit %s" % c[:10])
                rc = 1
                continue
            dh, dc = b[1] - a[1], a[0] - b[0]
            allowed = helper_copy_exceptions().get(c)
            if dh > 0 and dc < dh and allowed is not None:
                if dh > allowed:
                    print("rules_scan: commit %s raises helper_copy by %d, above its documented exception count %d in %s" % (c[:10], dh, allowed, EXCEPTIONS_FILE))
                    rc = 1
                else:
                    print("rules_scan: commit %s raises helper_copy by %d (documented exception, allowed %d, in %s)" % (c[:10], dh, allowed, EXCEPTIONS_FILE))
                continue
            if dh > 0 and dc < dh:
                print("rules_scan: commit %s raises helper_copy by %d but lowers clones by only %d (helper_copy may rise only with an equal-or-larger clone drop)" % (c[:10], dh, dc))
                rc = 1
    if rc == 0:
        print("rules_scan: per-commit helper_copy rule holds on %d commit(s)" % len(commits))
    return rc


# ---- borrowing ops: no copy anywhere in their call graph ----
# ExactNum::{add_ref, sub_ref, mul_ref, div_ref} read both operands in place and allocate only the result. This rule
# walks their static call graph (type-aware resolution, see the BORROW RULE in the header) and fails on a copy token in
# any reached body, printing the reach chain.
BORROW_ROOTS = (("zenith-float-num/src/ext.rs", ("add_ref", "sub_ref", "mul_ref", "div_ref")),)
# the whole crate src (non-test library code): a helper anywhere in it can be reached
BORROW_SCOPE = ("zenith-float-num/src/",)
# Copy mechanisms (see the BORROW RULE). Spelling-independent pieces are regexes; extend(...) and vec![e; n] / [e; n]
# need their argument parsed (borrow_copies).
BORROW_COPY_METHODS = ("clone", "clone_from", "clone_into", "clone_from_slice", "copy_from_slice", "swap_with_slice", "cloned", "copied",
                       "to_owned", "into_owned", "extend_from_slice", "extend_from_within", r"to_vec\w*")
BORROW_LEAD = r"(?<![|\w)\]])\|(?!\|)"  # a closure's opening pipe (not `||`, not a bitor after an operand)
BORROW_BANNED = re.compile(
    r"(?:\.|::)\s*(?:%s)\b(?:\s*\()?" % "|".join(BORROW_COPY_METHODS)  # method calls, UFCS paths and fn values, any spacing
    + r"|\b(?:Vec|Box|Rc|Arc|String|VecDeque|SmallVec|WordBuf|From|FromIterator)\s*(?:::\s*<[^>]*>\s*)?::\s*(?:from|from_iter)\b"
    + r"|\b(?:as_slice|as_mut_slice|as_ref|borrow|deref)\s*\(\s*\)\s*\.\s*(?:try_)?into\s*\("  # slice -> owned conversion
    + r"|\[[^\[\]]*\.\.[^\[\]]*\]\s*\.\s*(?:try_)?into\s*\("  # x[a..b].into()
    + r"|" + BORROW_LEAD + r"\s*\(?\s*&"  # |&x| / |(&a, &b)|: a closure that binds an item by value
    + r"|" + BORROW_LEAD + r"[^|;{}]*\|\s*(?:->\s*[^{|;]*?)?\{?\s*\*"  # |x| *x, |x: &T| -> T { *x }: a deref closure
    + r"|\bfor\s+\(?\s*&|\blet\s+(?:Some|Ok)?\s*\(?\s*&\s*(?:mut\s+)?[a-z_(]"  # for &x in / let &x = / if let Some(&x)
    + r"|\b(?:%s)\b" % "|".join(no_clone.HELPERS))
BORROW_LITERAL = re.compile(r"-?\d[\w.]*|(?:\w+\s*::\s*)*[A-Z][A-Z0-9_]*|true|false|None|(?:\w+\s*::\s*)*(?:new|default|zero)\s*\(\s*\)")
BORROW_OWNED_ITER = re.compile(r"\binto_iter\s*\(|\bdrain\s*\(|\brepeat\w*\s*\(|\.\s*map\s*\(|\bfrom_fn\s*\(|\bonce\s*\(|\bsuccessors\s*\(|^\s*\(?\s*[\w.]+\s*\.\.")


def _close(code, i):
    """index of the bracket closing the one at code[i]"""
    depth = 0
    for j in range(i, len(code)):
        if code[j] in "([{":
            depth += 1
        elif code[j] in ")]}":
            depth -= 1
            if depth == 0:
                return j
    return len(code)


# Named exceptions: (file, fn, token, max sites, reason). A copy site of the borrow call graph passes only when its
# (file, fn, token) is listed here and that fn holds no more such sites than listed; anything else fails. Kept to the
# copies an algorithm needs: the result's own words, a working buffer the kernel mutates, scalar word reads.
M = "zenith-float-num/src/mantissa/"
SCALAR = "scalar Word read into a local (one u64 value in a register); no buffer is duplicated"
BORROW_EXCEPTIONS = (
    ("zenith-float-num/src/common/buf.rs", "resize_len", ".extend_from_slice(", 1,
     "inline->heap spill of the result's own WordBuf: inline words cannot change storage by ownership"),
    ("zenith-float-num/src/common/int.rs", "copy_from", ".copy_from_slice(", 1,
     "SliceWithSign::copy_from: FFT accumulator / mul.rs intermediates written into a kernel's &mut working slice"),
    (M + "base.rs", "abs_add", ".copy_from_slice(", 2,
     "addend below the precision: the result is the other operand's words plus a sticky bit"),
    (M + "div.rs", "div_basic", ".clone_from_slice(", 2,
     "schoolbook division with d == 1: the remainder buffer the kernel subtracts in place starts as the operands"),
    (M + "div.rs", "div_basic", ".copy_from_slice(", 1, "writes the computed remainder into the result buffer"),
    (M + "div.rs", "div_recursive", ".copy_from_slice(", 7,
     "Burnikel-Ziegler: the remainder (= dividend when it is short) and quotient/remainder blocks assembled in "
     "working buffers the recursion mutates in place"),
    (M + "div.rs", "div_unbalanced", ".copy_from_slice(", 5,
     "unbalanced division: working remainder starts as the dividend and is reduced in place; quotient blocks "
     "assembled into the result"),
    (M + "fft.rs", "fft_decompose", ".copy_from_slice(", 3,
     "FFT input: operand words split into zero-padded parts the transform overwrites in place"),
    (M + "fft.rs", "fft_mul", ".copy_from_slice(", 1, "writes a reduced product part back into the transform buffer"),
    (M + "toom3.rs", "toom3", ".copy_from_slice(", 3, "writes the evaluated parts s0/s4 into the result buffer d3"),
    (M + "base.rs", "find_bit_len", "for&", 1, SCALAR),
    (M + "base.rs", "leading_zero_bits", "for&", 1, SCALAR),
    (M + "base.rs", "round_mantissa", "for&", 1, SCALAR),
    (M + "base.rs", "word_at", "letSome(&w", 1, SCALAR),
    (M + "util.rs", "next", "letSome(&h", 2, SCALAR + " (NormalizedView)"),
    (M + "util.rs", "next", "letSome(&l", 1, SCALAR + " (NormalizedView)"),
    (M + "util.rs", "next", "letSome(&v", 2, SCALAR + " (NormalizedView)"),
    ("zenith-float-num/src/ieee_soft/array.rs", "get", ".copied(", 1,
     SCALAR + "; IEEE soft-float bits, reached only by the name fallback of `.get(`"),
    ("zenith-float-num/src/ieee_soft/array.rs", "map", "|&", 1,
     SCALAR + "; IEEE soft-float bits, reached only by the name fallback of `.map(`"),
)


def borrow_copies(body):
    """[(offset, token)] of every copy mechanism in a stripped body"""
    out = [(m.start(), re.sub(r"\s+", "", m.group(0))) for m in BORROW_BANNED.finditer(body)]
    for m in re.finditer(r"\.\s*extend\s*\(", body):  # Extend<&T>: only an owned iterator may feed it
        a = m.end() - 1
        arg = body[a + 1:_close(body, a)]
        if arg.lstrip().startswith("&") or not BORROW_OWNED_ITER.search(arg):
            out.append((m.start(), ".extend(<borrowed>)"))
    for m in re.finditer(r"(\bvec\s*!\s*)?\[", body):  # vec![e; n] and [e; n] repeat operand data unless e is a literal
        a = m.end() - 1
        inner, depth, k = body[a + 1:_close(body, a)], 0, None
        for q, ch in enumerate(inner):
            depth += ch in "([{"; depth -= ch in ")]}"
            if ch == ";" and depth == 0:
                k = q; break
        if k is None:
            continue
        e = inner[:k].strip()
        if m.group(1) and not BORROW_LITERAL.fullmatch(e) or not m.group(1) and re.search(r"[.\[]", e):
            out.append((m.start(), "vec![<operand>;n]" if m.group(1) else "[<operand>;n]"))
    return sorted(out)
def strip_rust(src):
    """src with comments (line, nested block, doc) and string / char literal contents blanked; newlines and
    columns kept, so offsets map to the same lines. Lifetimes ('a) stay."""
    out, i, n = [], 0, len(src)
    blank = lambda s: "".join("\n" if c == "\n" else " " for c in s)
    while i < n:
        c = src[i]
        if src.startswith("//", i):
            j = src.find("\n", i)
            j = n if j < 0 else j
            out.append(blank(src[i:j])); i = j
        elif src.startswith("/*", i):
            depth, j = 1, i + 2
            while j < n and depth:
                if src.startswith("/*", j):
                    depth, j = depth + 1, j + 2
                elif src.startswith("*/", j):
                    depth, j = depth - 1, j + 2
                else:
                    j += 1
            out.append(blank(src[i:j])); i = j
        elif re.match(r'b?r#*"', src[i:i + 40]) and (i == 0 or not (src[i - 1].isalnum() or src[i - 1] == "_")):
            m = re.match(r'(b?r)(#*)"', src[i:])
            end = src.find('"' + m.group(2), i + m.end())
            j = n if end < 0 else end + 1 + len(m.group(2))
            out.append('""' + blank(src[i + 2:j])); i = j
        elif c == '"' or (c == "b" and src.startswith('b"', i) and (i == 0 or not (src[i - 1].isalnum() or src[i - 1] == "_"))):
            j = i + (2 if c == "b" else 1)
            while j < n and src[j] != '"':
                j += 2 if src[j] == "\\" else 1
            j = min(j + 1, n)
            out.append('""' + blank(src[i + 2:j])); i = j
        elif c == "'":
            m = re.match(r"'(?:\\(?:u\{[0-9a-fA-F]+\}|x[0-9a-fA-F]{2}|.)|[^\\'\n])'", src[i:])
            if m:
                out.append("' '" + blank(src[i + 3:i + m.end()])); i += m.end()
            else:
                out.append(c); i += 1
        else:
            out.append(c); i += 1
    return "".join(out)


BORROW_TYPE = r"&?\s*(?:'[a-z_]\w*\s+)?(?:mut\s+)?(?:(?:crate|super|self|[a-z_]\w*)\s*::\s*)*([A-Z]\w*|\[|[a-z]\w*)"
# std types whose methods are not crate code (Vec derefs to a slice: crate impls on slices are `[]`)
BORROW_EXTERN = {"Vec": "[]", "String": "str"}
BORROW_FIELD = re.compile(r"\b([a-z_]\w*)\s*:(?!:)\s*(?:pub(?:\([^)]*\))?\s+)?" + BORROW_TYPE)
BORROW_PARAM = re.compile(r"\b([a-z_]\w*)\s*:(?!:)\s*" + BORROW_TYPE)
BORROW_LET_T = re.compile(r"\blet\s+(?:mut\s+)?([a-z_]\w*)\s*(?::\s*" + BORROW_TYPE + r"|=\s*&?\s*(?:mut\s+)?([A-Z]\w*)\s*(?:::|\{|\())")
BORROW_LOCAL = re.compile(r"\b(?:let|for)\s+(?:mut\s+)?([a-z_]\w*)|\|\s*&?\s*(?:mut\s+)?([a-z_]\w*)\s*[|,:]|\b([a-z_]\w*)\s*:(?!:)")
BORROW_PAT = re.compile(r"\blet\s+([^=;]*)=|(?:[{,]|=>)\s*([^{},;=]*?)=>")
BORROW_ITEM = re.compile(r"^[ \t]*(?:pub(?:\([^)]*\))?\s+)?(?:unsafe\s+)?((?:impl|trait)\b)", re.M)
BORROW_FN = re.compile(r"\bfn\s+([a-z_][a-z0-9_]*)|\bmacro_rules!\s*([a-z_][a-z0-9_]*)")
BORROW_REF = re.compile(r"(<[^<>]*(?:<[^<>]*(?:<[^<>]*>[^<>]*)*>[^<>]*)*>\s*::\s*|(?:\b(?:[A-Za-z_]\w*)\s*(?:::\s*<[^<>]*(?:<[^<>]*>[^<>]*)*>)?\s*::\s*)+)?"
                        r"(\.\s*)?\b([a-z_][a-z0-9_]*)\b(\s*!)?")
# operator -> trait method names (dispatched only when an operand is a local / param of a crate type with that impl)
BORROW_OPS = (("+", ("add", "add_assign")), ("-", ("sub", "sub_assign", "neg")), ("*", ("mul", "mul_assign")),
              ("/", ("div", "div_assign")), ("%", ("rem", "rem_assign")))


def _impl_type(header):
    """`impl<T> Tr<X> for &'a Ty<T>` -> 'Ty' (the self type), or the trait name for `trait Tr`."""
    h = header.strip()
    if h.startswith("trait"):
        m = re.match(r"trait\s+([A-Z]\w*)", h)
        return m.group(1) if m else None
    h = h[4:].lstrip()
    if h.startswith("<"):
        depth = 0
        for k, ch in enumerate(h):
            depth += (ch == "<") - (ch == ">")
            if depth == 0:
                h = h[k + 1:]
                break
    m = re.search(r"\bfor\s+(.*)", h, re.S)
    if m:
        h = m.group(1)
    m = re.search(BORROW_TYPE, h.split(" where ")[0].strip())
    return _norm_type(m.group(1)) if m else None


def _norm_type(t):
    """a type head from BORROW_TYPE -> the key impls are indexed under ('[]' for slices, None when unknown)"""
    if t == "[":
        return "[]"
    if t in ("impl", "dyn", "mut", "fn"):
        return None
    return BORROW_EXTERN.get(t, t)


def _span_end(code, j):
    """index of the `}` matching the first `{` at or after j (or of a `;` that ends a body-less item first)"""
    depth, started = 0, False
    while j < len(code):
        ch = code[j]
        if ch == "{":
            depth, started = depth + 1, True
        elif ch == "}":
            depth -= 1
        if started and depth <= 0:
            return j
        if not started and ch == ";":
            return j
        j += 1
    return len(code) - 1


fields = {}  # struct name -> {field: type key}, filled by borrow_fns


def borrow_fns():
    """[Fn] for every non-test fn and macro_rules! of BORROW_SCOPE: dict(name, rel, line, body, bline, ty, free,
    file_stem). Bodies come from brace matching over comment/string-stripped text; ty is the enclosing impl's self type
    (or trait), free = not inside an impl/trait block. Macros are named `name!`."""
    out = []
    for rel, p in no_clone.rs_files(["zenith-float-num/src"]):
        if not rel.startswith(BORROW_SCOPE) or not is_lib_only(rel):
            continue
        code = strip_rust(open(p, errors="replace").read())
        lines = code.split("\n")
        mask = no_clone.test_mask(lines)
        code = "\n".join("" if mask[i] else l for i, l in enumerate(lines))
        impls = []
        for m in BORROW_ITEM.finditer(code):
            b = code.find("{", m.start(1))
            if b < 0 or ";" in code[m.start(1):b]:
                continue
            impls.append((m.start(1), _span_end(code, b), _impl_type(code[m.start(1):b])))
        for m in re.finditer(r"\bstruct\s+([A-Z]\w*)[^{;(]*\{", code):
            for g in BORROW_FIELD.findall(code[m.end():_span_end(code, m.end() - 1)]):
                fields.setdefault(m.group(1), {})[g[0]] = _norm_type(g[1])
        stem = re.sub(r"(/mod)?\.rs$", "", rel).split("/")
        for m in BORROW_FN.finditer(code):
            name = m.group(1) or (m.group(2) + "!")
            j = _span_end(code, m.end())
            encl = [t for s, e, t in impls if s < m.start() < e]
            ar, slf = _sig(code[m.end():j + 1]) if m.group(1) else (None, False)
            out.append(dict(arity=ar, self=slf, name=name, rel=rel, line=code.count("\n", 0, m.start()) + 1, body=code[m.end():j + 1],
                            bline=code.count("\n", 0, m.end()) + 1, ty=encl[-1] if encl else None,
                            free=not encl, mods=set(stem)))
    return out


def _known_types(f):
    """local / param name -> crate-visible type name, from the signature and `let` lines; self -> impl type"""
    kt = {}
    for g in BORROW_PARAM.findall(f["body"]):
        kt[g[0]] = _norm_type(g[1])
    for g in BORROW_LET_T.findall(f["body"]):
        kt[g[0]] = _norm_type(g[1] or g[2])
    if f["ty"]:
        kt["self"] = f["ty"]
    for k, v in list(kt.items()):
        if v == "Self":
            kt[k] = f["ty"]
    return kt


def _top_args(text, i):
    """number of top-level comma-separated items in the parenthesised list opening at text[i] == '(' (None when a
    closure makes the split ambiguous)"""
    depth, n, nonempty, j = 0, 0, False, i
    while j < len(text):
        ch = text[j]
        if ch in "([{" or (ch == "<" and depth > 0 and not text[j - 1] in " =<"):
            depth += 1
        elif ch in ")]}" or (ch == ">" and depth > 1 and text[j - 1] not in "-="):
            depth -= 1
            if depth == 0:
                return n + 1 if nonempty else 0
        elif ch == "|" and depth == 1:
            return None
        elif ch == "," and depth == 1:
            if not text[j + 1:].lstrip().startswith(")"):
                n += 1
        elif depth == 1 and not ch.isspace():
            nonempty = True
        j += 1
    return None


def _sig(body):
    """(param count, takes self) of a fn from its text after the name, or (None, False)"""
    j, depth = 0, 0
    while j < len(body) and body[j] != "(":
        if body[j] == "{":
            return None, False
        j += 1
    if j >= len(body):
        return None, False
    k = body.find(")", j)
    takes_self = bool(re.match(r"\(\s*(?:&\s*(?:'[a-z_]\w*\s+)?)?(?:mut\s+)?self\b", body[j:]))
    return _top_args(body, j), takes_self


def _recv_type(before, kt):
    """type key of the receiver chain `a.b.c.` that ends `before`, through locals and struct fields; None if unknown"""
    r = re.search(r"((?:[A-Za-z_]\w*\s*\.\s*)*[A-Za-z_]\w*)\s*\.\s*$", before)
    if not r or re.search(r"[.)\]?:>]\s*$", before[:r.start()]):
        return None
    chain = [x.strip() for x in r.group(1).split(".")]
    ty = kt.get(chain[0])
    for fld in chain[1:]:
        ty = fields.get(ty, {}).get(fld) if ty else None
    return ty


# std *types* (not traits: `Add::add`, `From::from`, `Clone::clone` dispatch into crate impls and are followed)
BORROW_STD_TYPES = {"Vec", "Box", "String", "Option", "Result", "Some", "Ok", "Err", "Rc", "Arc", "Mutex", "RwLock",
                    "RefCell", "Cell", "VecDeque", "BTreeMap", "HashMap", "Ordering", "Wrapping", "PhantomData"}


def _close_angle(code, i):
    """index of the `>` closing the `<` at code[i] (`->` is not a bracket)"""
    depth = 0
    for j in range(i, len(code)):
        if code[j] == "<":
            depth += 1
        elif code[j] == ">" and code[j - 1] != "-":
            depth -= 1
            if depth == 0:
                return j
    return len(code) - 1


def _callees(f, by_name, types):
    """names-and-types resolution of the fns a body can reach (an over-approximation within each rule):
    `name!` -> macros; `T::name` / `<T as Tr>::name` -> fns of impl T (Self -> own impl); `mod::name` -> free fns of
    files under that module; `.name(` on self / a typed local -> that type's fns, on anything else -> every fn of that
    name in an impl; bare `name(` or `name` as a value -> free fns (and own-impl fns); operators on a typed operand ->
    that type's op-trait fns."""
    body, kt = f["body"], _known_types(f)
    local = {x for g in BORROW_LOCAL.findall(body) for x in g if x}
    for pat in BORROW_PAT.findall(body):  # let / if let / match-arm patterns bind every lowercase ident in them
        local.update(re.findall(r"\b[a-z_]\w*\b", " ".join(pat)))
    res = []
    for m in BORROW_REF.finditer(body):
        path, dot, name, bang = m.group(1), m.group(2), m.group(3), m.group(4)
        if bang:
            res += [g for g in by_name.get(name + "!", [])]
            continue
        cands = by_name.get(name, [])
        if not cands:
            continue
        end = m.end()
        tf = re.match(r"\s*::\s*<", body[end:])
        if tf:  # turbofish: name::<T>( is a (generic) call, name::<T>:: a path
            end = _close_angle(body, end + tf.end() - 1) + 1
        rest = body[end:end + 3].lstrip()
        if rest.startswith("::"):
            continue  # a path segment (module or type), not a fn
        call = rest.startswith("(")
        if call:  # Rust has no overloading or default args: a callee must take exactly this many args
            na = _top_args(body, body.index("(", end))
            if na is not None:
                cands = [g for g in cands if g["arity"] is None
                         or g["arity"] - (1 if dot and g["self"] else 0) == na
                         or (not dot and g["self"] and g["arity"] == na)]
        if path:
            bare = path
            while re.search(r"<[^<>]*>", bare):  # drop turbofish generics: Vec::<T>::new is Vec::new
                bare = re.sub(r"(?:::\s*)?<[^<>]*>", "", bare)
            segs = re.findall(r"[A-Za-z_]\w*", bare)
            qual = re.match(r"<\s*&?\s*(?:mut\s+)?([A-Z]\w*)", path)  # <T as Tr>::f / <T>::f
            last = qual.group(1) if qual else segs[-1]
            if last == "Self":
                last = f["ty"]
            if segs and segs[0] in ("core", "std", "alloc"):
                continue  # outside the crate: no crate body to reach
            if last in BORROW_STD_TYPES:
                continue  # an inherent fn of a std type (Vec::new): std code calls no crate fn back
            if last and last[0].isupper():  # Type::f; a trait, generic param or default method -> every method f
                res += [g for g in cands if g["ty"] == last] or [g for g in cands if not g["free"]]
            else:  # mod::f; an unknown module (a `use` alias) -> every free f
                res += ([g for g in cands if g["free"] and (last in g["mods"] or last in ("crate", "super", "self"))]
                        or [g for g in cands if g["free"]])
        elif dot:
            if not call:
                continue  # a field
            ty = _recv_type(body[:m.start() + len(dot)], kt)
            typed = [g for g in cands if g["ty"] == ty] if ty else []
            # unknown receiver, a generic / trait-bound one, or a type without its own such method (trait impls for
            # slices or primitives, Deref, trait objects): every method of that name
            res += typed or [g for g in cands if not g["free"]]
        else:
            if not call and (name in local or re.match(r"\s*(?:[:{!]|=[^=>])", body[end:])
                             or re.search(r"\b(?:fn|let|mut|struct|enum|mod|use|as)\s*$|[\'.]\s*$", body[:m.start()])):
                continue  # a fn value is any other use of a free fn's name: binding, argument, field init, return
            res += [g for g in cands if g["free"] or (g["ty"] and g["ty"] == f["ty"] and not call)]
    for op, names in BORROW_OPS:
        rx = re.compile(r"\b([a-z_]\w*)\s*%s=?(?![>=%s])|(?<![-=<>!&|*/%%+])%s=?\s*&?\s*(?:mut\s+)?([a-z_]\w*)"
                        % (re.escape(op), re.escape(op), re.escape(op)))
        tys = {kt[x] for g in rx.findall(body) for x in g if x in kt}
        for n in names:
            res += [g for g in by_name.get(n, []) if g["ty"] in tys]
    return res


def borrow_closure():
    """-> (reached [(rel, line, name)], hits [(rel, line, name, token, path)]); path = root -> ... -> name (BFS)"""
    fns = borrow_fns()
    by_name = collections.defaultdict(list)
    for g in fns:
        by_name[g["name"]].append(g)
    types = {g["ty"] for g in fns if g["ty"]}
    root_rel = dict((n, rel) for rel, names in BORROW_ROOTS for n in names)
    todo = collections.deque((g, (g["name"],)) for rel, names in BORROW_ROOTS for n in names
                             for g in by_name.get(n, []) if g["rel"] == rel)
    seen, reached, hits = set(), [], []
    while todo:
        f, path = todo.popleft()
        key = (f["rel"], f["line"], f["name"])
        if key in seen:
            continue
        seen.add(key)
        reached.append(key)
        for off, tok in borrow_copies(f["body"]):
            hits.append((f["rel"], f["bline"] + f["body"].count("\n", 0, off), f["name"], tok, path))
        for g in _callees(f, by_name, types):
            if (g["rel"], g["line"], g["name"]) not in seen:
                todo.append((g, path + (g["name"],)))
    return reached, hits


def borrow_check(update=False):
    if not os.path.exists(os.path.join(ROOT, BORROW_ROOTS[0][0])):
        return 0
    reached, hits = borrow_closure()
    roots = [n for _, names in BORROW_ROOTS for n in names]
    missing = [n for n in roots if not any(r[2] == n for r in reached)]
    for n in missing:
        print("rules_scan: borrowing op %s not found (BORROW_ROOTS is stale)" % n)
    allowed = {(r, f, t): n for r, f, t, n, _ in BORROW_EXCEPTIONS}
    groups = collections.Counter((rel, name, tok) for rel, _, name, tok, _ in hits)
    bad = [h for h in hits if groups[(h[0], h[2], h[3])] > allowed.get((h[0], h[2], h[3]), 0)]
    for rel, line, name, tok, path in bad:
        print("rules_scan: BORROW %s:%d in fn %s: `%s` (%s)" % (rel, line, name, tok, " -> ".join(path)))
    for k, n in sorted(allowed.items()):
        if groups.get(k, 0) < n:
            print("rules_scan: note: borrow exception %s %s `%s` allows %d, found %d: lower it" % (k + (n, groups.get(k, 0))))
    print("rules_scan: borrowing ops reach %d fns, %d copy sites (%d named exceptions, %d unlisted)"
          % (len(reached), len(hits), len(hits) - len(bad), len(bad)))
    return 1 if bad or missing else 0


def both_count_checks(update=False):
    rc = max(count_check(update), helper_count_check(update), exception_check(update),
             test_panic_check(update), neg_via_mul_check(update), ruling_check(update), borrow_check(update))
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
    if "--print-counts" in argv:
        print("%d %d" % (clone_count(), helper_copy_count()))
        return 0
    if "--borrow-reach" in argv:
        for rel, line, name in sorted(borrow_closure()[0]):
            print("%s:%d %s" % (rel, line, name))
        return borrow_check()
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
        rc = max(rc, both_count_checks(), per_commit_check(rng))
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
