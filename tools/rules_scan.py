#!/usr/bin/env python3
"""rules_scan: scan the ADDED lines of a diff for banned constructs in library code.

Run on EVERY model-written diff before review, and in tools/../full_gate (step `rules-scan`).

Usage:
  tools/rules_scan.sh                      # origin/main..HEAD, plus the clone-count baseline check
  tools/rules_scan.sh A..B                 # any range (a bare rev means REV..HEAD)
  tools/rules_scan.sh --staged             # git diff --cached
  tools/rules_scan.sh --working            # git diff HEAD (staged + unstaged)
  tools/rules_scan.sh --diff-file F        # a saved unified diff (used by the self-test)
  tools/rules_scan.sh --count-only         # only the clone-count and add_ref/sub_ref-count baseline checks
  tools/rules_scan.sh --update-baseline    # lower tools/clone_count_baseline.txt to the current count (never raises it)

LIBRARY PATHS (clone, to_vec, cloned, to_owned, unsafe, allow, unwrap, expect, panic): .rs files under
  zenith-float-num/src/  zenith-float-macro/src/  zenith-float-compare/src/  src/
  excluding any path with a tests/ benches/ examples/ tools/ component, files named tests.rs, and
  any added line inside a #[cfg(test)] item of the file at the diff's tip.
CALCULATION PATHS (f32, f64): the library paths above except zenith-float-compare/src/ (benchmark adapters).

Flagged on added lines only (comments and string literals are ignored):
  .clone(  .to_vec(  .cloned(  .to_owned(  unsafe  #[allow / #![allow  unwrap(  expect(  panic!  f32  f64

ALLOWLIST: tools/rules_allowlist.txt, one `path:hash<TAB>reason` per line; the key is printed as allow-key=
for every flagged site (hash = sha1 of the added line with whitespace collapsed). The file starts EMPTY.
ONLY JEFF MAY EXTEND IT (.github/CODEOWNERS requires his review). Models never add entries.

CLONE COUNT: total `.clone(` occurrences in non-test library code (the four library dirs) must never go
up versus tools/clone_count_baseline.txt. The baseline is only ever lowered, with --update-baseline.

REF-SUM COUNT: total uses of `add_ref(` / `sub_ref(` (the borrowed-borrowed sums, sanctioned only where both
operands are still live) in non-test library code, definitions (`fn add_ref`) excluded. Same ratchet as the
clone count against tools/ref_sum_count_baseline.txt: it may only go down; --update-baseline lowers it.
"""
import hashlib, os, re, subprocess, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("RULES_SCAN_ROOT") or os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ["NO_CLONE_ROOT"] = ROOT
import no_clone  # test_mask, DIRS

LIB = re.compile(r"^(zenith-float-num/src|zenith-float-macro/src|zenith-float-compare/src|src)/.*\.rs$")
CALC = re.compile(r"^(zenith-float-num/src|zenith-float-macro/src|src)/.*\.rs$")
EXCL_PARTS = {"tests", "benches", "examples", "tools"}
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
    parts = path.split("/")
    if not LIB.match(path) or parts[-1] == "tests.rs":
        return False
    return not (EXCL_PARTS & set(parts[:-1]))


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
    total = 0
    for d in no_clone.DIRS:
        for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
            dn[:] = [x for x in dn if x not in ("target", ".git")]
            for f in fn:
                if not f.endswith(".rs"):
                    continue
                rel = os.path.relpath(os.path.join(dp, f), ROOT)
                if not is_lib(rel):
                    continue
                lines = open(os.path.join(dp, f), errors="replace").read().split("\n")
                mask = no_clone.test_mask(lines)
                for i, l in enumerate(lines):
                    if not mask[i]:
                        total += len(re.findall(r"\.clone\(", strip(l)))
    return total


REF_SUM = re.compile(r"\b(?:add_ref|sub_ref)\(")


def ref_sum_count():
    total = 0
    for d in no_clone.DIRS:
        for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
            dn[:] = [x for x in dn if x not in ("target", ".git")]
            for f in fn:
                if not f.endswith(".rs"):
                    continue
                rel = os.path.relpath(os.path.join(dp, f), ROOT)
                if not is_lib(rel):
                    continue
                lines = open(os.path.join(dp, f), errors="replace").read().split("\n")
                mask = no_clone.test_mask(lines)
                for i, l in enumerate(lines):
                    code = strip(l)
                    if not mask[i] and not re.search(r"\bfn\s+(add_ref|sub_ref)\b", code):
                        total += len(REF_SUM.findall(code))
    return total


def baseline_path():
    return os.path.join(ROOT, "tools", "clone_count_baseline.txt")


def ref_baseline_path():
    return os.path.join(ROOT, "tools", "ref_sum_count_baseline.txt")


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
    cur = ref_sum_count()
    try:
        base = int(open(ref_baseline_path()).read().split()[0])
    except (OSError, ValueError, IndexError):
        print("rules_scan: missing or unreadable tools/ref_sum_count_baseline.txt")
        return 1
    if cur > base:
        print("rules_scan: add_ref/sub_ref use count %d is ABOVE the baseline %d (it may only go down)" % (cur, base))
        return 1
    if cur < base:
        if update:
            open(ref_baseline_path(), "w").write("%d\n" % cur)
            print("rules_scan: ref-sum baseline lowered %d -> %d" % (base, cur))
        else:
            print("rules_scan: add_ref/sub_ref use count %d is below baseline %d: lower the baseline (rules_scan.sh --update-baseline)" % (cur, base))
        return 0
    print("rules_scan: add_ref/sub_ref use count %d == baseline" % cur)
    return 0


def both_count_checks(update=False):
    return max(count_check(update), ref_count_check(update))


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
