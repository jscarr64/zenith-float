#!/usr/bin/env bash
# Proves tools/mutants_apply_check.sh can fail: a fixture repo with one patch that applies and one stale patch (its context is gone)
# must FAIL naming only the stale id; with the stale row removed it must pass. Then runs the real check on this checkout.
# Why: 11 arith patches went stale at c52b43b and only the 3h mutants step noticed (INVALID). usage: tools/mutants_apply_check_selftest.sh
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap mutants_apply_check_selftest "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
here="$(cd "$(dirname "$0")" && pwd)"; . "$here/scratch.sh"
fail() { echo "mutants-apply-check selftest FAIL: $*"; exit 1; }
t=$(mktemp -d "$(zb_scratch)/mutants-apply-check-selftest.XXXXXX") || fail "mktemp failed"; runlog_on_exit "rm -rf \"$t\""
r=$t/repo; mkdir -p "$r/src" "$r/tools/mutants/manifest.d" "$r/tools/mutants/s" || fail "mkdir"
git -C "$r" init -q && printf 'a\nb\nc\n' > "$r/src/x.rs" && git -C "$r" add . && git -C "$r" -c user.name=t -c user.email=t@t commit -qm x || fail "fixture repo"
printf 'a\nB\nc\n' > "$r/src/x.rs"; git -C "$r" diff > "$r/tools/mutants/s/001-OK.patch"
printf 'a\nb\nc\n' > "$r/src/x.rs"
printf -- '--- a/src/x.rs\n+++ b/src/x.rs\n@@ -1,3 +1,3 @@\n a\n-gone\n+B\n c\n' > "$r/tools/mutants/s/002-STALE.patch"
printf '# c\nsrc/x.rs:f\t001-OK\ts/001-OK.patch\tharness\n' > "$r/tools/mutants/manifest.txt"
printf 'src/x.rs:f\t002-STALE\ts/002-STALE.patch\tharness\n' > "$r/tools/mutants/manifest.d/s.txt"
out=$("$here/mutants_apply_check.sh" "$r" 2>&1); rc=$?
[ "$rc" -ne 0 ] || fail "stale patch not flagged (rc 0): $out"
grep -qx 'STALE 002-STALE: s/002-STALE.patch does not apply' <<<"$out" || fail "stale id not named: $out"
! grep -q '001-OK' <<<"$out" || fail "good patch flagged: $out"
: > "$r/tools/mutants/manifest.d/s.txt"
out=$("$here/mutants_apply_check.sh" "$r" 2>&1) || fail "all-good fixture failed: $out"
git -C "$r" diff --quiet || fail "check modified the fixture tree"
"$here/mutants_apply_check.sh" "$here/.." || fail "a manifest patch does not apply to this checkout (above)"
echo "mutants-apply-check selftest: PASS"
