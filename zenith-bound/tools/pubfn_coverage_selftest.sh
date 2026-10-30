#!/usr/bin/env bash
# Self-test for tools/pubfn_coverage.sh: builds fixture trees under /tmp at run time and proves the check
# (1) lists exactly the expected public items, (2) PASSES on a clean fixture, (3) FAILS on each failure mode
# (UNCOVERED, STALE, NOTEST), (4) honours the exemption list for UNCOVERED only, reads a present list (6) and handles an empty one with no ERR (7), and FAILS STALE-EXEMPT on an exempt entry that is not an item (8). No repo files are touched.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap pubfn_coverage_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d "$(zb_scratch)"/pubfn-selftest.XXXXXX)"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
chk() { PUBFN_ROOT="$1" PUBFN_EXEMPT="${2:-$1/tools/pubfn_coverage_exempt.txt}" "$here/pubfn_coverage.sh" "${@:3}"; }

mkfix() { # mkfix <dir>: a clean, fully covered fixture crate
  local d="$1"; mkdir -p "$d/src/n" "$d/tests" "$d/tools"
  cat > "$d/src/lib.rs" <<'RS'
//! fixture
pub mod m;
pub mod n;
pub fn top() {}
RS
  cat > "$d/src/m.rs" <<'RS'
pub struct Foo { x: u8 }
impl Foo {
    pub fn new() -> Foo { Foo { x: 0 } }
    pub const fn c() -> u8 { 1 }
    fn private() {}
    pub(crate) fn hidden() {}
}
impl Ord for Foo {
    fn cmp(&self, o: &Foo) -> core::cmp::Ordering { self.x.cmp(&o.x) }
}
pub fn free() -> u8 { let s = "pub fn in_a_string() {}"; 0 }
// pub fn in_a_comment() {}
pub enum E {
    A,
    B(u8),
    C { z: u8 },
}
RS
  printf 'pub mod inner {\n    pub fn deep() {}\n}\n' > "$d/src/n/mod.rs"
  cat > "$d/tests/t.rs" <<'RS'
//! COVERS: crate::top, Foo::new, Foo::c, Ord for Foo, m::free, E::A, E::B, E::C, n::inner::deep
#[test]
fn covers_everything() -> Result<(), String> {
    top(); Foo::new(); Foo::c(); let _ = Ord::cmp; free(); let _ = (E::A, E::B(1), E::C { z: 1 }); deep();
    Ok(())
}
RS
  : > "$d/tools/pubfn_coverage_exempt.txt"
}

want='Foo::c
Foo::new
E::A
E::B
E::C
Ord for Foo
crate::top
m::free
n::inner::deep'
want="$(printf '%s\n' "$want" | LC_ALL=C sort)"

# 1. item extraction
c="$t/clean"; mkfix "$c"
ck; got="$(chk "$c" "" --list)" || fail "--list failed"
ck; [ "$got" = "$want" ] || fail "item list differs. got:
$got
want:
$want"

# 2. clean fixture passes
ck; out="$(chk "$c")" || fail "clean fixture failed: $out"
ck; echo "$out" | grep -c "ok (9 public items)" > /dev/null || fail "clean fixture: unexpected summary: $out"

# 3a. UNCOVERED: a new public fn that no COVERS line names
u="$t/uncovered"; mkfix "$u"; printf 'impl Foo {\n    pub fn extra() {}\n}\n' >> "$u/src/m.rs"
ck; out="$(chk "$u")" && fail "uncovered public fn passed"
ck; echo "$out" | grep -c "UNCOVERED: Foo::extra" > /dev/null || fail "uncovered: wrong report: $out"
# 4. exemption rescues UNCOVERED
printf 'Foo::extra\tfixture exemption\n' > "$u/tools/pubfn_coverage_exempt.txt"
ck; chk "$u" >/dev/null || fail "exempt item still failed"
# 3b. STALE: a COVERS name that no longer exists (an exemption must NOT rescue it)
s="$t/stale"; mkfix "$s"; sed -i 's/^\(\/\/! COVERS:.*\)$/\1, Foo::gone/' "$s/tests/t.rs"
printf 'Foo::gone\tshould not rescue\n' > "$s/tools/pubfn_coverage_exempt.txt"
ck; out="$(chk "$s")" && fail "stale COVERS name passed"
ck; echo "$out" | grep -c "STALE: tests/t.rs names Foo::gone" > /dev/null || fail "stale: wrong report: $out"
# 3c. NOTEST: declared, exists, but no #[test] in the declaring file mentions it (the mention sits in a non-test helper)
n="$t/notest"; mkfix "$n"
cat > "$n/tests/t.rs" <<'RS'
//! COVERS: crate::top, Foo::new, Foo::c, Ord for Foo, m::free, E::A, E::B, E::C, n::inner::deep
fn helper() { top(); Foo::new(); Foo::c(); free(); deep(); let _ = (E::A, E::B(1), E::C { z: 1 }); let _ = Ord::cmp; }
#[test]
fn covers_nothing_by_name() -> Result<(), String> { helper(); Ok(()) }
RS
ck; out="$(chk "$n")" && fail "name-less test passed"
ck; echo "$out" | grep -c "NOTEST: tests/t.rs declares Foo::new" > /dev/null || fail "notest: wrong report: $out"
# 3c'. NOTEST: a file with COVERS and no #[test] at all
z="$t/notestfile"; mkfix "$z"; printf '//! COVERS: crate::top, Foo::new, Foo::c, Ord for Foo, m::free, E::A, E::B, E::C, n::inner::deep\n' > "$z/tests/t.rs"
ck; out="$(chk "$z")" && fail "file with no #[test] passed"
ck; echo "$out" | grep -c "NOTEST" > /dev/null || fail "no-test file: wrong report: $out"

# 5. word match: a test that mentions only a longer identifier (new_thing) must not satisfy Foo::new
w="$t/word"; mkfix "$w"; sed -i 's/Foo::new();/Foo::new_thing();/' "$w/tests/t.rs"
ck; out="$(chk "$w")" && fail "substring-only mention passed (grep -w semantics broken)"
ck; echo "$out" | grep -c "NOTEST: tests/t.rs declares Foo::new" > /dev/null || fail "word: wrong report: $out"
# 6. exemptions present and read: a comment-headed list with two entries rescues exactly those two UNCOVERED items;
#    dropping the list (an empty exempt) must make the same tree FAIL, so a silently skipped load cannot pass
x="$t/exread"; mkfix "$x"; printf 'impl Foo {\n    pub fn ex1() {}\n    pub fn ex2() {}\n}\n' >> "$x/src/m.rs"
printf '# header comment\n\nFoo::ex1\tfixture\n  # indented comment\nFoo::ex2\tfixture\n' > "$x/tools/pubfn_coverage_exempt.txt"
ck; out="$(chk "$x")" || fail "two read exemptions still failed: $out"
ck; echo "$out" | grep -c "exempt list: 2 items" > /dev/null || fail "exempt read: count not reported: $out"
: > "$t/empty_exempt"
ck; out="$(chk "$x" "$t/empty_exempt")" && fail "dropped exemptions passed: $out"
ck; echo "$out" | grep -c "UNCOVERED: Foo::ex2" > /dev/null || fail "dropped exemptions: wrong report: $out"
# 7. empty exempt list (comments and blanks only, as shipped) is explicit and logs no ERR. The run joins this
#    self-test's own log as a nested sub-run so an ERR trap line from pubfn_coverage.sh would land there.
e="$t/emptyex"; mkfix "$e"; printf '# only a comment\n\n' > "$e/tools/pubfn_coverage_exempt.txt"
before=$(wc -l < "$POPLOG_FILE")
ck; out="$(export POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE; chk "$e")" || fail "empty exempt list failed: $out"
ck; echo "$out" | grep -c "exempt list is empty" > /dev/null || fail "empty exempt: not reported explicitly: $out"
ck; if tail -n "+$((before + 1))" "$POPLOG_FILE" | grep "ERR rc=" | grep -q "pubfn_coverage.sh"; then fail "empty exempt: ERR line logged: $(tail -n "+$((before + 1))" "$POPLOG_FILE" | grep "ERR rc=")"; fi
# 8. STALE-EXEMPT: an exempt entry that names no public item must FAIL (a clean, fully covered tree otherwise)
q="$t/staleexempt"; mkfix "$q"; printf 'Foo::ghost\tno such item\n' > "$q/tools/pubfn_coverage_exempt.txt"
ck; out="$(chk "$q")" && fail "stale exempt entry passed: $out"
ck; echo "$out" | grep -c "FAIL STALE-EXEMPT: Foo::ghost is in the exempt list but is not a public item" > /dev/null || fail "stale exempt: wrong report: $out"
st_done pubfn_coverage 24   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "pubfn_coverage selftest: ok"
