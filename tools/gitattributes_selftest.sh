#!/usr/bin/env bash
# Self-test for .gitattributes (ZF-TD41): goldens must be byte-exact (-text)
# and ordinary text files must check out LF. Fails on a tree with no
# .gitattributes, which is the red state for the .gitattributes change.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
fail() { echo "GITATTRIBUTES SELFTEST FAIL: $*"; exit 1; }
[ -f "$here/.gitattributes" ] || fail "no .gitattributes at repo root"
attr() { git -C "$here" check-attr -a -- "$1"; }
golden="$(attr zenith-float-num/tests/oracle/goldens/arith.jsonl)"
echo "$golden" | grep -q "text: unset" || fail "arith.jsonl check-attr: $golden"
text="$(attr README.md)"
echo "$text" | grep -q "eol: lf" || fail "README.md check-attr: $text"
echo "gitattributes selftest: ok"
