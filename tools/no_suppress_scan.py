#!/usr/bin/env python3
"""no_suppress_scan: zero lint suppressions and zero unsafe code in zenith-float (Jeff 2026-11-24 1:08 PM CT: "we should
have 0 [allowed] dead code or unsafe code so the [allowed] shjould not be anywhere inside of zenith float whatsoever";
1:17 PM CT go). Scans the WHOLE TREE (every git-tracked file at the working tree), not a diff. There is no allowlist and
no reason file: the only passing count is 0.

Usage: tools/no_suppress_scan.sh [--root DIR]      (exit 0 clean, 1 findings, 2 scanner error)

FAILS on (file:line [form] text), every git-tracked file, all crates, tests/, benches/, examples/, fuzz/, nostd-tests/:
  [allow] [expect]      `#[allow(..)]`, `#![allow(..)]`, `#[expect(..)]`, `#![expect(..)]` in any .rs (comments, strings
                        and char literals stripped first, so prose and fixtures in strings do not count)
  [cfg_attr-allow]      `#[cfg_attr(<pred>, allow(..))]` / `expect(..)`, inner or outer
  [unsafe]              the token `unsafe` in code, anywhere except the known-pending inventory below
  [cargo-lint]          any `[lints.*]` / `[workspace.lints.*]` entry of any Cargo.toml at level "allow" or "warn"
                        (string or `{ level = .. }` form)
  [clippy-config]       any clippy.toml / .clippy.toml key that relaxes a lint by a number: `*-threshold`, `max-*`,
                        `*-limit` (e.g. too-many-arguments-threshold, max-fn-params-bools, type-complexity-threshold)
  [lint-flag]           `-A <lint>`, `-A<lint>`, `--allow`, `--cap-lints` passed to rustc/cargo/clippy (RUSTFLAGS,
                        RUSTDOCFLAGS, cargo or clippy command lines) in scripts/, tools/, .github/, .githooks/, .cargo/
                        and top-level shell/make/just files; `.cargo/config*` rustflags arrays holding "-A"/"--cap-lints"
  [no-forbid]           a crate root (each package's lib target and bin targets) without `#![forbid(unsafe_code)]`,
                        unless its Cargo.toml has `[lints.rust] unsafe_code = "forbid"` (or `[lints] workspace = true`
                        with that in `[workspace.lints.rust]`)

KNOWN-PENDING UNSAFE (Paul call 2026-11-24, recorded in WORK2; Jeff 1:17 PM CT: "unsafe tests as a separate step after"):
tools/unsafe_pending.tsv lists `path<TAB>count<TAB>note` for the test-only unsafe that a later step removes (the MPFR FFI
tests under zenith-float-num/tests/mpfr/ and nostd-tests/src/bin/nostd-host.rs). The scan prints every pending file with
its count and the total. A pending file must hold EXACTLY its pinned count: more fails (it cannot grow), fewer fails too
(lower the pin in the same commit, so it cannot grow back). A library file (a crate's src/ outside nostd-tests) can never
be pending, and an unsafe line in any file not on the list fails. A pending bin root is exempt from [no-forbid] only
while it is listed. Only lines can leave this file; nothing is ever added to it.

Test: tools/no_suppress_selftest.sh (fixture repos: each planted form fails by name; comments, strings, raw strings,
forbid/deny attributes and identifiers like unsafe_code pass).
"""
import os, re, sys, subprocess, tomllib

PENDING_FILE = "tools/unsafe_pending.tsv"
SELF = {"tools/no_suppress_scan.py", "tools/no_suppress_selftest.sh", "tools/no_suppress_scan.sh"}
LIB_DIRS = ("src/", "zenith-float-num/src/", "zenith-float-bound/src/", "zenith-float-macro/src/",
            "zenith-float-compare/src/", "zenith-float-testkit/src/")


def strip_rust(src):
    """comments, string / raw string / byte string and char literals -> spaces (newlines kept)."""
    out, i, n = [], 0, len(src)
    blank = lambda s: re.sub(r"[^\n]", " ", s)
    while i < n:
        c = src[i]
        if src.startswith("//", i):
            j = src.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i)); i = j; continue
        if src.startswith("/*", i):
            d, j = 0, i
            while j < n:
                if src.startswith("/*", j):
                    d += 1; j += 2
                elif src.startswith("*/", j):
                    d -= 1; j += 2
                else:
                    j += 1
                if d == 0:
                    break
            out.append(blank(src[i:j])); i = j; continue
        prev_ident = i > 0 and (src[i - 1].isalnum() or src[i - 1] == "_")
        m = re.match(r'(b?r)(#*)"', src[i:i + 70]) if not prev_ident else None
        if m:
            end = '"' + m.group(2)
            j = src.find(end, i + len(m.group(0)))
            j = n if j < 0 else j + len(end)
            out.append(blank(src[i:j])); i = j; continue
        if c == '"' or (c == "b" and src.startswith('b"', i) and not prev_ident):
            j = i + (2 if c == "b" else 1)
            while j < n and src[j] != '"':
                j += 2 if src[j] == "\\" else 1
            j += 1
            out.append(blank(src[i:j])); i = j; continue
        if c == "'":
            m = re.match(r"'(?:\\u\{[0-9a-fA-F]+\}|\\x[0-9a-fA-F]{2}|\\.|[^\\'\n])'", src[i:i + 14])
            if m:
                out.append(" " * len(m.group(0))); i += len(m.group(0)); continue
        out.append(c); i += 1
    return "".join(out)


def line_of(s, pos):
    return s.count("\n", 0, pos) + 1


def balanced(s, j):
    """s[j-1] is '(' -> index just past its matching ')'"""
    d = 1
    while j < len(s) and d:
        d += {"(": 1, ")": -1}.get(s[j], 0); j += 1
    return j


ATTR = re.compile(r"#\s*(!)?\s*\[\s*(allow|expect|cfg_attr)\s*\(")
UNSAFE = re.compile(r"\bunsafe\b")


def scan_rs(rel, raw, find):
    s = strip_rust(raw)
    lines = raw.split("\n")
    for m in ATTR.finditer(s):
        end = balanced(s, m.end())
        body = s[m.end():end - 1]
        kind = m.group(2)
        ln = line_of(s, m.start())
        if kind == "cfg_attr":
            if not re.search(r"\b(allow|expect)\s*\(", body):
                continue
            kind = "cfg_attr-allow"
        find(rel, ln, kind, lines[ln - 1].strip())
    return [line_of(s, m.start()) for m in UNSAFE.finditer(s)], s


def is_lib(rel):
    return rel.endswith(".rs") and rel.startswith(LIB_DIRS) and "/tools/" not in "/" + rel


def lint_level(v):
    if isinstance(v, str):
        return v
    if isinstance(v, dict):
        return v.get("level")
    return None


FLAG_FILE = re.compile(r"(^|/)(scripts|tools|\.github|\.githooks|\.cargo)/|^[^/]*\.(sh|bash|mk)$|^(Makefile|justfile|Justfile)$")
FLAG_CTX = re.compile(r"RUSTFLAGS|RUSTDOCFLAGS|\bcargo\b|\bclippy\b|\brustc\b|rustflags")
FLAG = re.compile(r"(?<![\w-])(-A\s*[\"']?[a-z_]+(?:::[a-z_]+)?\b|--allow\b|--cap-lints\b)")


def main(argv):
    root = None
    if "--root" in argv:
        root = argv[argv.index("--root") + 1]
    else:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        files = subprocess.run(["git", "-C", root, "ls-files", "-z"], check=True, capture_output=True,
                               text=True).stdout.split("\0")
    except (subprocess.CalledProcessError, OSError) as e:
        print(f"no_suppress_scan: ERROR git ls-files in {root}: {e}", file=sys.stderr); return 2
    files = sorted(f for f in files if f and os.path.isfile(os.path.join(root, f)))
    findings = []
    find = lambda f, ln, kind, text: findings.append(f"{f}:{ln} [{kind}] {text[:160]}")
    pending = {}
    pp = os.path.join(root, PENDING_FILE)
    if os.path.exists(pp):
        for k, l in enumerate(open(pp, encoding="utf-8"), 1):
            l = l.rstrip("\n")
            if not l.strip() or l.lstrip().startswith("#"):
                continue
            t = l.split("\t")
            if len(t) < 3 or not t[1].isdigit() or not t[2].strip():
                findings.append(f"{PENDING_FILE}:{k} [pending] malformed (want path<TAB>count<TAB>note): {l[:120]}")
                continue
            pending[t[0]] = int(t[1])
    for f in pending:
        if is_lib(f):
            findings.append(f"{PENDING_FILE}: {f} [pending] is a library file: library unsafe is never pending")
    unsafe_by_file, stripped = {}, {}
    pkgs = []
    for f in files:
        p = os.path.join(root, f)
        base = os.path.basename(f)
        if f.endswith(".rs"):
            raw = open(p, encoding="utf-8", errors="replace").read()
            us, s = scan_rs(f, raw, find)
            stripped[f] = s
            if us:
                unsafe_by_file[f] = us
        elif base == "Cargo.toml":
            try:
                doc = tomllib.loads(open(p, encoding="utf-8").read())
            except tomllib.TOMLDecodeError as e:
                findings.append(f"{f}:0 [cargo-lint] unparsable Cargo.toml: {e}"); continue
            txt = open(p, encoding="utf-8").read().split("\n")
            for top, tab in (("lints", doc.get("lints")), ("workspace.lints", (doc.get("workspace") or {}).get("lints"))):
                if not isinstance(tab, dict):
                    continue
                for tool, entries in tab.items():
                    if not isinstance(entries, dict):
                        continue
                    for name, v in entries.items():
                        lv = lint_level(v)
                        if lv in ("allow", "warn"):
                            ln = next((i + 1 for i, t in enumerate(txt) if re.match(r"\s*" + re.escape(name) + r"\s*=", t)), 0)
                            find(f, ln, "cargo-lint", f"{top}.{tool}.{name} = {lv}")
            if "package" in doc:
                pkgs.append((os.path.dirname(f), doc))
        elif base in ("clippy.toml", ".clippy.toml"):
            for k, l in enumerate(open(p, encoding="utf-8").read().split("\n"), 1):
                m = re.match(r"\s*([A-Za-z0-9_-]+)\s*=", l)
                if m and re.search(r"(threshold$|^max-|-limit$)", m.group(1)):
                    find(f, k, "clippy-config", l.strip())
        if f in SELF or f.endswith(".rs"):
            continue
        if FLAG_FILE.search(f):
            try:
                content = open(p, encoding="utf-8").read()
            except UnicodeDecodeError:
                continue
            cfg = "/.cargo/" in "/" + f
            for k, l in enumerate(content.split("\n"), 1):
                code = l.split("#", 1)[0] if f.endswith((".sh", ".bash", ".toml", ".yml", ".yaml")) or "/.githooks/" in "/" + f else l
                if cfg and re.search(r"[\"'](-A|--cap-lints|--allow)", code):
                    find(f, k, "lint-flag", l.strip()); continue
                if FLAG_CTX.search(code) and FLAG.search(code):
                    find(f, k, "lint-flag", l.strip())
    # unsafe: library never, elsewhere only the pinned pending inventory, exact counts
    total_pending = 0
    for f, lns in sorted(unsafe_by_file.items()):
        if f in pending and not is_lib(f):
            if len(lns) != pending[f]:
                findings.append(f"{f}: {len(lns)} unsafe, pinned {pending[f]} [pending] "
                                f"({'it grew' if len(lns) > pending[f] else 'lower the pin in tools/unsafe_pending.tsv'})")
            total_pending += len(lns)
            continue
        for ln in lns:
            find(f, ln, "unsafe", open(os.path.join(root, f), encoding="utf-8", errors="replace").read().split("\n")[ln - 1].strip())
    for f in pending:
        if f not in unsafe_by_file and not is_lib(f):
            findings.append(f"{f}: 0 unsafe, pinned {pending[f]} [pending] (remove the line from {PENDING_FILE})")
    # crate roots must forbid unsafe_code
    ws = {}
    rt = os.path.join(root, "Cargo.toml")
    if os.path.exists(rt):
        ws = (tomllib.loads(open(rt, encoding="utf-8").read()).get("workspace") or {})
    ws_forbid = lint_level(((ws.get("lints") or {}).get("rust") or {}).get("unsafe_code")) == "forbid"
    for d, doc in pkgs:
        lints = doc.get("lints") or {}
        forbid = lint_level((lints.get("rust") or {}).get("unsafe_code")) == "forbid" or \
            (lints.get("workspace") is True and ws_forbid)
        j = lambda x: os.path.normpath(os.path.join(d, x)) if d else os.path.normpath(x)
        roots = []
        lib = doc.get("lib") or {}
        if lib.get("path"):
            roots.append(j(lib["path"]))
        elif j("src/lib.rs") in stripped:
            roots.append(j("src/lib.rs"))
        for b in doc.get("bin") or []:
            if b.get("path"):
                roots.append(j(b["path"]))
            elif b.get("name") and j(f"src/bin/{b['name']}.rs") in stripped:
                roots.append(j(f"src/bin/{b['name']}.rs"))
        if j("src/main.rs") in stripped:
            roots.append(j("src/main.rs"))
        roots += [f for f in stripped if f.startswith(j("src/bin") + "/") and f.count("/") == j("src/bin").count("/") + 1]
        for r in sorted(set(roots)):
            if r not in stripped:
                continue
            if forbid or re.search(r"#!\s*\[\s*forbid\s*\(\s*([^)]*\b)?unsafe_code\b", stripped[r]):
                continue
            if r in pending:
                continue
            findings.append(f"{r} [no-forbid] crate root without #![forbid(unsafe_code)] (and no Cargo [lints.rust] unsafe_code = \"forbid\")")
    for x in findings:
        print(x)
    if pending:
        print(f"no_suppress_scan: known-pending unsafe (separate step; exact pins in {PENDING_FILE}): "
              f"{len(pending)} files, {total_pending} sites")
        for f in sorted(pending):
            print(f"  pending {f}\t{len(unsafe_by_file.get(f, []))}")
    if findings:
        print(f"no_suppress_scan: FAIL {len(findings)} finding(s); fix the code, nothing is allow-listed")
        return 1
    print(f"no_suppress_scan: ok ({len(files)} files; 0 allow/expect, 0 lint relaxations, 0 unsafe outside the pending inventory)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
