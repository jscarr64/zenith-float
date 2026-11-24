#!/usr/bin/env bash
# Selftest for .githooks/pre-push (Jeff 2026-11-24 1:17 PM CT: "pre-push hook"). The hook runs, for every pushed
# commit, the CI rules-scan command (tools/rules_scan.sh "${RATCHET_BASE:-c74e834^}..<sha>", .github/workflows/
# rules-scan.yml) and tools/no_suppress_scan.py on that commit's tree. Fixture: a --shared clone of this repo's HEAD
# with one base commit that drops clippy.toml threshold keys (so the base is no_suppress-clean); RATCHET_BASE = that
# base, so each case scans only its own fixture commit.
# usage: tools/pre_push_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
hook="$root/.githooks/pre-push"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
w=$(mktemp -d "${TMPDIR:-/tmp}/pre-push-selftest.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
if [ ! -x "$hook" ]; then
  echo "FAIL: $hook missing or not executable"; echo "pre_push_selftest: FAIL"; exit 1
fi
c="$w/clone"
git clone -q --shared "$root" "$c" || exit 2
g() { git -C "$c" -c user.name=t -c user.email=t@t "$@"; }
orig=$(g rev-parse HEAD)
sed -i -E '/^[a-z-]*(threshold|max-[a-z-]*|-limit) *=/d' "$c/clippy.toml"
g commit -qam "fixture base: no clippy thresholds"; base=$(g rev-parse HEAD)
mkc() { # mkc FILE CONTENT -> sha of a commit on base adding FILE
  g checkout -q --detach "$base"; mkdir -p "$c/$(dirname "$1")"; printf '%b' "$2" > "$c/$1"
  g add -A; g commit -qm "case $1"; g rev-parse HEAD
}
zero=0000000000000000000000000000000000000000
# case WANT LABEL SHA [OUTPUT-REGEX]
case_() {
  local want=$1 label=$2 sha=$3 rx=${4:-}; n=$((n + 1))
  local out rc
  out=$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$sha" "$zero" \
        | RATCHET_BASE="$base" bash "$hook" origin "$w/remote" 2>&1); rc=$?
  if [ "$want" = pass ] && [ $rc -ne 0 ]; then bad "$label: expected pass, rc=$rc"; echo "$out" | tail -5
  elif [ "$want" = fail ] && [ $rc -eq 0 ]; then bad "$label: expected fail, rc=0"
  elif [ -n "$rx" ] && ! grep -Eq "$rx" <<<"$out"; then bad "$label: output lacks /$rx/"; echo "$out" | tail -5
  else echo "ok $n $label (rc=$rc)"; fi
}
clean=$(mkc zenith-float-num/src/zz_fixture.rs '/// fixture\npub fn zz(x: u32) -> u32 {\n    x + 1\n}\n')
clone=$(mkc zenith-float-num/src/zz_fixture.rs 'pub fn zz(x: &Vec<u8>) -> Vec<u8> {\n    x.clone()\n}\n')
tallow=$(mkc tests/zz_fixture.rs '#[allow(dead_code)]\nfn zz() {}\n')
lallow=$(mkc zenith-float-num/src/zz_fixture.rs '#[allow(clippy::too_many_arguments)]\npub fn zz() {}\n')
unsafe_=$(mkc zenith-float-num/src/zz_fixture.rs 'pub fn zz(p: *const u8) -> u8 {\n    unsafe { *p }\n}\n')
g checkout -q --detach "$tallow"   # the working tree is NOT the pushed commit: the hook must scan the pushed tree
case_ pass "clean commit passes both scans" "$clean" 'pre-push: ok'
case_ fail "added .clone() in library code fails rules-scan" "$clone" 'FAIL rules-scan'
case_ fail "#[allow] in tests/ passes rules-scan but fails no-suppress" "$tallow" 'ok rules-scan.*|FAIL no-suppress'
case_ fail "#[allow] in library code fails" "$lallow" 'FAIL (rules-scan|no-suppress)'
case_ fail "unsafe in library code fails no-suppress" "$unsafe_" 'FAIL no-suppress'
case_ fail "a relaxed clippy.toml threshold fails no-suppress (this repo's HEAD)" "$orig" '\[clippy-config\]'
n=$((n + 1)); out=$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$zero" "$clean" | RATCHET_BASE="$base" bash "$hook" origin "$w/remote" 2>&1)
if [ $? -eq 0 ]; then echo "ok $n branch delete (nothing to scan) passes"; else bad "branch delete: $out"; fi
grep -q 'FAIL no-suppress' <<<"$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$tallow" "$zero" | RATCHET_BASE="$base" bash "$hook" origin x 2>&1)" \
  && grep -q 'ok rules-scan' <<<"$(cd "$c" && printf 'refs/heads/x %s refs/heads/x %s\n' "$tallow" "$zero" | RATCHET_BASE="$base" bash "$hook" origin x 2>&1)" \
  && echo "ok $((n + 1)) tests/ #[allow]: rules-scan ok, no-suppress FAIL (why the hook runs both)" || bad "tests/ allow split"
[ $fail -eq 0 ] && echo "pre_push_selftest: ok ($n cases + split)" || echo "pre_push_selftest: FAIL"
exit $fail
