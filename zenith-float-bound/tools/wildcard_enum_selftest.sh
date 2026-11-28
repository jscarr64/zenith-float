#!/usr/bin/env bash
# S10 check (TEST_PLAN P5, section 21): clippy::wildcard_enum_match_arm must be DENIED for every target (src, tests), so a new rounding
# mode or enum variant can never fall through a `_ =>` or a `other =>` arm. The deny lives in Cargo.toml [lints.clippy]; `cargo clippy
# --all-targets -- -D warnings` (the gate's clippy step) then enforces it on the real tree. `match_wildcard_for_single_variants = "deny"` closes the
# case wildcard_enum_match_arm exempts (a `_` that stands for exactly one variant). This script proves the enforcement can fail:
#   1. Cargo.toml [lints.clippy] contains `wildcard_enum_match_arm = "deny"`;
#   2. a scratch crate (under /tmp, built from the repo's own [lints.*] blocks) with a BAD sample (`_ =>` and a binding arm over an enum)
#      FAILS clippy with the lint named;
#   3. the same crate with the lint line stripped PASSES the bad sample (so step 2 is really caused by the line);
#   4. a GOOD sample (every variant named) passes with the repo's lints;
#   5. Cargo.toml also has `match_wildcard_for_single_variants = "deny"`; a SINGLE-variant sample (`A => .., _ => ..` over a two-variant enum) FAILS
#      with the repo's lints and PASSES with that line stripped.
# No repo file is touched. Offline, no dependencies. usage: tools/wildcard_enum_selftest.sh
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap wildcard_enum_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
export LC_ALL=C
root=${WILDCARD_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
fail() { echo "wildcard-enum selftest FAIL: $*"; exit 1; }
[ -f "$root/Cargo.toml" ] || fail "no Cargo.toml in $root"
t=$(mktemp -d "$(zb_scratch)"/wildcard-selftest.XXXXXX); runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"

lints() { # the [lints.*] tables of the repo's Cargo.toml; arg: keep | strip (no wildcard_enum_match_arm line) | stripsingle (no match_wildcard_for_single_variants line)
  awk '/^\[lints\./ { on = 1 } /^\[/ && !/^\[lints\./ { on = 0 } on' "$root/Cargo.toml" | {
    case "${1:-keep}" in strip) grep -v '^wildcard_enum_match_arm' ;; stripsingle) grep -v '^match_wildcard_for_single_variants' ;; *) cat ;; esac
  }
}
ck; grep -q '^wildcard_enum_match_arm = "deny"$' "$root/Cargo.toml" || fail "Cargo.toml [lints.clippy] has no wildcard_enum_match_arm = \"deny\" line"
ck; grep -q '^match_wildcard_for_single_variants = "deny"$' "$root/Cargo.toml" || fail "Cargo.toml [lints.clippy] has no match_wildcard_for_single_variants = \"deny\" line"
ck; lints | grep -c '^\[lints.clippy\]' > /dev/null || fail "Cargo.toml has no [lints.clippy] table"

mkcrate() { # mkcrate <dir> <strip|keep> <bad|good>
  mkdir -p "$1/src"
  cp "$root/rust-toolchain.toml" "$1/" 2>/dev/null
  { printf '[package]\nname = "wildcard_sample"\nversion = "0.0.0"\nedition = "2021"\n\n[dependencies]\n\n'; lints "$2"; } > "$1/Cargo.toml"
  if [ "$3" = single ]; then
    cat > "$1/src/lib.rs" <<'RS'
pub enum Two { A, B }
pub fn single(m: &Two) -> u8 {
    match m {
        Two::A => 1,
        _ => 0,
    }
}
RS
  elif [ "$3" = bad ]; then
    cat > "$1/src/lib.rs" <<'RS'
pub enum Mode { A, B, C }
pub fn wild(m: &Mode) -> u8 {
    match m {
        Mode::A => 1,
        _ => 0,
    }
}
pub fn binding(m: Mode) -> Mode {
    match m {
        Mode::A => Mode::B,
        other => other,
    }
}
RS
  else
    cat > "$1/src/lib.rs" <<'RS'
pub enum Mode { A, B, C }
pub fn named(m: &Mode) -> u8 {
    match m {
        Mode::A => 1,
        Mode::B | Mode::C => 0,
    }
}
RS
  fi
}
run() { (cd "$1" && CARGO_TARGET_DIR="$1/target" cargo clippy --offline --all-targets -- -D warnings 2>&1); }

mkcrate "$t/bad" keep bad;  out=$(run "$t/bad");  rc=$?
ck; [ "$rc" != 0 ] || fail "the bad sample passed clippy with the repo's lints: $out"
ck; echo "$out" | grep -c 'wildcard match will also match any future added variants\|wildcard_enum_match_arm' > /dev/null || fail "the bad sample failed, but not on the wildcard lint: $out"
mkcrate "$t/stripped" strip bad; out=$(run "$t/stripped"); rc=$?
ck; [ "$rc" = 0 ] || fail "the bad sample fails even without the wildcard lint line (so the lint line proves nothing): $out"
mkcrate "$t/good" keep good; out=$(run "$t/good"); rc=$?
ck; [ "$rc" = 0 ] || fail "the good sample was rejected: $out"
mkcrate "$t/single" keep single; out=$(run "$t/single"); rc=$?
ck; [ "$rc" != 0 ] || fail "the single-variant sample passed clippy with the repo's lints: $out"
ck; echo "$out" | grep -c 'wildcard matches only a single variant' > /dev/null || fail "the single-variant sample failed, but not on the single-variant lint: $out"
mkcrate "$t/single2" stripsingle single; out=$(run "$t/single2"); rc=$?
ck; [ "$rc" = 0 ] || fail "the single-variant sample fails even without the single-variant lint line (so that line proves nothing): $out"
mkcrate "$t/wildonly" strip single; out=$(run "$t/wildonly"); rc=$?
ck; [ "$rc" != 0 ] || fail "the single-variant sample passes with only the single-variant lint stripped from the other side: the two lints are not independent: $out"
st_done wildcard-enum 11   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "wildcard-enum selftest: ok (both deny lines present; bad and single-variant samples rejected; stripped configs accept them; good sample accepted)"
