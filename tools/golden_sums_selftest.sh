#!/usr/bin/env bash
# Self-test for tools/golden_sums_check.sh on a scratch copy (mktemp -d under $TMPDIR): it must fail closed when
# SHA256SUMS is missing, empty, does not list a required golden, or a hash mismatches; and pass on a correct pin.
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d)"; trap 'rm -rf "$t"' EXIT
fail() { echo "GOLDEN SUMS SELFTEST FAIL: $*"; exit 1; }
chk() { "$here/golden_sums_check.sh" "$t" borrow_ops.jsonl; }
printf '{"id":"a"}\n{"id":"b"}\n' > "$t/borrow_ops.jsonl"; printf 'x\n' > "$t/other.jsonl"
sum() { sha256sum "$t/$1" | cut -d' ' -f1; }
n=0
rm -f "$t/SHA256SUMS"; chk >/dev/null 2>&1 && fail "missing SHA256SUMS passed"; n=$((n + 1))
: > "$t/SHA256SUMS"; chk >/dev/null 2>&1 && fail "empty SHA256SUMS passed"; n=$((n + 1))
printf '%s  other.jsonl\n' "$(sum other.jsonl)" > "$t/SHA256SUMS"; chk >/dev/null 2>&1 && fail "SHA256SUMS without borrow_ops.jsonl passed"; n=$((n + 1))
printf '%s  borrow_ops.jsonl\n' "$(sum borrow_ops.jsonl)" > "$t/SHA256SUMS"; printf '{"id":"c"}\n' >> "$t/borrow_ops.jsonl"
chk >/dev/null 2>&1 && fail "drifted borrow_ops.jsonl passed"; n=$((n + 1))
printf '%s  borrow_ops.jsonl\n' "$(sum borrow_ops.jsonl)" > "$t/SHA256SUMS"
out="$(chk 2>&1)" || fail "correct pin failed: $out"; n=$((n + 1))
[ "$n" = 5 ] || fail "ran $n cases, want 5"
echo "golden_sums selftest: ok ($n cases)"
