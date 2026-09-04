#!/usr/bin/env python3
"""Shared parsers for zenith-float prepublish checks (leaves, Precision rustdoc, hygiene)."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

IDENT = re.compile(r"`([a-zA-Z_][a-zA-Z0-9_]*)`")
LEAF_ARM = re.compile(r'"([a-zA-Z_][a-zA-Z0-9_]*)"')
UNWRAP = re.compile(r"\bunwrap\s*\(")
TODO = re.compile(r"\b(TODO|FIXME|HACK)\b")
PROOF = re.compile(r"unwrap", re.IGNORECASE)


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def section_body(text: str, heading_substr: str) -> str:
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith("## ") and heading_substr.lower() in line.lower():
            start = i
            break
    if start is None:
        raise SystemExit(f"missing heading containing {heading_substr!r}")
    end = len(lines)
    for j in range(start + 1, len(lines)):
        if lines[j].startswith("## "):
            end = j
            break
    return "\n".join(lines[start:end])


def backtick_idents_first_column(md_table: str) -> list[str]:
    names: list[str] = []
    for line in md_table.splitlines():
        if not line.startswith("|") or line.startswith("| ---") or line.startswith("|--"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if not cells:
            continue
        if cells[0].lower() in {"method", "item", "api"}:
            continue
        for span in re.findall(r"`([^`]+)`", cells[0]):
            m = re.match(r"[A-Za-z_][A-Za-z0-9_]*", span)
            if m:
                names.append(m.group(0))
    return names


def complete_list_after(text: str, marker: str) -> list[str]:
    idx = text.find(marker)
    if idx < 0:
        raise SystemExit(f"missing marker {marker!r}")
    rest = text[idx + len(marker) :]
    tick = rest.find("`")
    if tick < 0:
        raise SystemExit(f"no backtick list after {marker!r}")
    end = rest.find("\n\n", tick)
    block = rest[tick:] if end < 0 else rest[tick:end]
    return IDENT.findall(block)


def match_arm_names(src: str) -> list[str]:
    """Leaf names from `\"foo\" =>` / `\"a\" | \"b\" =>` in a match."""
    names: list[str] = []
    for line in src.splitlines():
        if "=>" not in line:
            continue
        left = line.split("=>", 1)[0]
        if '"' not in left:
            continue
        names.extend(LEAF_ARM.findall(left))
    return names


def rust_files(root: Path) -> list[Path]:
    out: list[Path] = []
    for sub in (
        root / "zenith-float-num" / "src",
        root / "zenith-float-macro" / "src",
        root / "src",
    ):
        if sub.is_dir():
            out.extend(sorted(sub.rglob("*.rs")))
    return [p for p in out if p.name != "tests.rs"]


def production_span(text: str) -> str:
    cut = text.find("#[cfg(test)]")
    return text if cut < 0 else text[:cut]


def has_precision_for(name: str, files: list[Path]) -> bool:
    fn = re.compile(rf"\b(?:pub\s+)?fn\s+{re.escape(name)}\s*\(")
    wrap = re.compile(rf"^\s*{re.escape(name)}\s*,\s*$")
    for path in files:
        lines = path.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines):
            hit = bool(fn.search(line)) or bool(wrap.match(line))
            if not hit:
                continue
            window = "\n".join(lines[max(0, i - 40) : i + 1])
            if "# Precision" in window:
                return True
    return False


def check_precision(root: Path) -> int:
    cap = (root / "doc" / "ZENITH_FLOAT_CAPABILITIES.md").read_text(encoding="utf-8")
    names = []
    names.extend(backtick_idents_first_column(section_body(cap, "Special functions")))
    # Unique, skip expr-only aliases that are not method names.
    skip = {"yes"}
    ordered = []
    seen = set()
    for n in names:
        if n in skip or n in seen:
            continue
        seen.add(n)
        ordered.append(n)
    src = rust_files(root)
    missing = [n for n in ordered if not has_precision_for(n, src)]
    if missing:
        print("error: missing # Precision rustdoc for:", ", ".join(missing), file=sys.stderr)
        return 1
    print(f"check_precision_docs: {len(ordered)} specials have # Precision")
    return 0


def extract_call_match(src: str) -> str:
    """The match that lists expr!/cexpr! function names (first `\"recip\" =>`)."""
    idx = src.find('"recip" =>')
    if idx < 0:
        raise SystemExit("no expr/cexpr call match starting at recip")
    return src[idx:]


def check_leaves(root: Path) -> int:
    lib = (root / "doc" / "LIBRARY.md").read_text(encoding="utf-8")
    eidx = lib.find("### `expr!(expression, context)`")
    cidx = lib.find("### `cexpr!(expression, context)`")
    if eidx < 0 or cidx < 0:
        print("error: LIBRARY.md missing expr!/cexpr! headings", file=sys.stderr)
        return 1
    expr_md = complete_list_after(lib[eidx:cidx], "**Function leaves (complete list):**")
    cexpr_md = complete_list_after(lib[cidx:], "**Function leaves (complete list):**")
    expr_rs = match_arm_names(extract_call_match((root / "zenith-float-macro" / "src" / "lib.rs").read_text()))
    cexpr_rs = match_arm_names(extract_call_match((root / "zenith-float-macro" / "src" / "cplx.rs").read_text()))
    failed = 0
    for label, doc, code in (
        ("expr!", expr_md, expr_rs),
        ("cexpr!", cexpr_md, cexpr_rs),
    ):
        ds, cs = set(doc), set(code)
        if ds != cs:
            failed = 1
            extra_doc = sorted(ds - cs)
            extra_code = sorted(cs - ds)
            if extra_doc:
                print(f"error: {label} in LIBRARY.md not in macro: {extra_doc}", file=sys.stderr)
            if extra_code:
                print(f"error: {label} in macro not in LIBRARY.md: {extra_code}", file=sys.stderr)
        else:
            print(f"check_leaves: {label} {len(ds)} leaves match")
    return failed


def nearby_proof(lines: list[str], i: int) -> bool:
    # Same-line comment is the proof (existing kernel style).
    if "//" in lines[i]:
        return True
    lo, hi = max(0, i - 2), min(len(lines), i + 3)
    for j in range(lo, hi):
        s = lines[j]
        if "//" not in s and "/*" not in s:
            continue
        comment = s.split("//", 1)[-1] if "//" in s else s
        if PROOF.search(comment):
            return True
    return False


def check_hygiene(root: Path) -> int:
    failed = 0
    for path in rust_files(root):
        rel = path.relative_to(root)
        prod = production_span(path.read_text(encoding="utf-8"))
        lines = prod.splitlines()
        for i, line in enumerate(lines):
            if TODO.search(line):
                print(f"error: {rel}:{i+1}: TODO/FIXME/HACK", file=sys.stderr)
                failed = 1
            if UNWRAP.search(line) and not nearby_proof(lines, i):
                print(f"error: {rel}:{i+1}: unwrap() without proof comment", file=sys.stderr)
                failed = 1
    if failed == 0:
        print("check_hygiene: no TODO/FIXME/HACK; production unwrap() has a proof comment")
    return failed


def check_version(root: Path) -> int:
    cargo = (root / "Cargo.toml").read_text(encoding="utf-8")
    cap = (root / "doc" / "ZENITH_FLOAT_CAPABILITIES.md").read_text(encoding="utf-8")
    m = re.search(r"(?m)^version\s*=\s*\"([^\"]+)\"", cargo)
    if not m:
        print("error: no version in Cargo.toml", file=sys.stderr)
        return 1
    ver = m.group(1)
    vm = re.search(r"\*\*Version:\*\*\s*(\S+)", cap)
    if not vm or vm.group(1) != ver:
        print(f"error: CAPABILITIES version {vm.group(1) if vm else '?'} != Cargo.toml {ver}", file=sys.stderr)
        return 1
    print(f"check_version: {ver}")
    return 0


def check_verify(root: Path) -> int:
    doc = root / "doc"
    for path in sorted(doc.glob("*.md")):
        text = path.read_text(encoding="utf-8")
        if "<!-- verify -->" in text:
            print(f"error: leftover <!-- verify --> in {path.name}", file=sys.stderr)
            return 1
    print("check_verify: no <!-- verify --> markers")
    return 0


def main() -> int:
    os.chdir(repo_root())
    root = repo_root()
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    dispatch = {
        "precision": check_precision,
        "leaves": check_leaves,
        "hygiene": check_hygiene,
        "version": check_version,
        "verify": check_verify,
    }
    if cmd not in dispatch:
        print("usage: prepublish_lib.py {precision|leaves|hygiene|version|verify}", file=sys.stderr)
        return 2
    return dispatch[cmd](root)


if __name__ == "__main__":
    sys.exit(main())
