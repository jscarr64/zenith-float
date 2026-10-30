#!/usr/bin/env bash
# Self-test for tools/rules_scan.sh (zenith-bound): every banned pattern is flagged on ADDED lines in src/, tests/, benches/ and examples/ alike (there is no test-code exemption), removed/context lines, comments, strings and non-Rust paths are not flagged, the allowlist works, git modes work, and both count baselines only go down.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap rules_scan_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d)"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
scan() { "$here/rules_scan.sh" "$@"; }
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
mkdiff() { # mkdiff <path> <added line...>  -> fixture diff on stdout
  local p="$1"; shift
  printf 'diff --git a/%s b/%s\n--- a/%s\n+++ b/%s\n@@ -1,1 +1,%d @@\n context\n' "$p" "$p" "$p" "$p" "$(( $# + 1 ))"
  for l in "$@"; do printf '+%s\n' "$l"; done
}
LIB=src/foo.rs
export RULES_SCAN_ROOT="$t"; mkdir -p "$t/tools"; : > "$t/tools/rules_allowlist.txt"; echo 0 > "$t/tools/clone_count_baseline.txt"; echo 0 > "$t/tools/ref_sum_count_baseline.txt"
export RULES_SCAN_SKIP_COUNT=1   # fixture diffs have no tree to count

mkdiff $LIB 'let a = 1 + 2;' 'let s = a.max(b);' > "$t/clean.diff"
ck; scan --diff-file "$t/clean.diff" >/dev/null || fail "clean diff flagged"

while IFS='|' read -r name line; do
  mkdiff $LIB "$line" > "$t/d.diff"
  ck; out="$(scan --diff-file "$t/d.diff")" && fail "not flagged: $name"
  ck; echo "$out" | grep -c "\[$name\]" > /dev/null || fail "wrong pattern label for $name: $out"
done <<'CASES'
clone|let y = x.clone();
clone_from|y.clone_from(&x);
clone_path|let y = Clone::clone(&x);
to_vec|let y = s.to_vec();
cloned|let y = it.cloned();
copied|let y = it.copied();
to_owned|let y = s.to_owned();
derive_clone|#[derive(Clone)]
derive_clone|#[derive(Debug, Clone, Copy)]
impl_clone|impl Clone for Big {
impl_clone|impl Copy for Big {
unsafe|unsafe { core::hint::unreachable_unchecked() }
allow|#[allow(dead_code)]
allow|#![allow(clippy::all)]
unwrap|let v = r.unwrap();
expect|let v = r.expect("boom");
panic|panic!("no");
unreachable|unreachable!();
todo|todo!();
unimplemented|unimplemented!();
unwrap_variant|let v = r.unwrap_or_default();
unwrap_variant|let v = r.unwrap_unchecked();
unwrap_variant|let e = r.unwrap_err();
unwrap_variant|let e = r.unwrap_err_unchecked();
unwrap_variant|let e = r.expect_err("x");
assert|assert!(a == b);
assert|assert_eq!(a, b);
assert|assert_ne!(a, b);
assert|assert_matches!(a, Some(_));
assert|debug_assert!(a);
assert|debug_assert_eq!(a, b);
f32|let a: f32 = 1.0;
f64|let a = 1.0_f64; let b: f64 = 2.0;
hwfloat|let r = a.mul_add(b, c);
hwfloat|use core::arch::x86_64::_mm_add_pd;
neg_sub_method|let d = a.sub(&b, p, rm).neg();
neg_sub_method|let d = a.sub_ref(&b).neg();
neg_sub_call|let d = neg(a.sub(&b));
neg_sub_call|let d = Neg::neg(Sub::sub(a, b));
neg_sub_op|let d = -(b - a);
neg_sub_op|let d = -(b-a);
neg_sub_op|let d = x + -(b.exp - a.exp);
neg_sub_op|return -(b - a);
neg_sub_op|let d = -a.sub(&b);
neg_sub_op|let d = -(a.sub(&b));
CASES

# S9: the neg-of-sub forms are refused in test paths as well, and clean subtractions are not flagged
for p in tests/props.rs tests/model/mod.rs; do
  mkdiff $p 'let d = a.sub(&b).neg();' > "$t/p.diff"
  ck; out="$(scan --diff-file "$t/p.diff")" && fail "neg_sub_method not flagged in test path: $p"
  ck; echo "$out" | grep -c "\[neg_sub_method\]" > /dev/null || fail "wrong pattern label for neg_sub_method in $p: $out"
  mkdiff $p 'let d = -(b - a);' > "$t/p.diff"
  ck; out="$(scan --diff-file "$t/p.diff")" && fail "neg_sub_op not flagged in test path: $p"
  ck; echo "$out" | grep -c "\[neg_sub_op\]" > /dev/null || fail "wrong pattern label for neg_sub_op in $p: $out"
done
mkdiff $LIB 'let d = a.sub(&b, p, rm);' 'let d = x - (a - b);' 'let d = -1 - (a - b);' 'let d = a.sub_lower(&b);' 'let n = x.neg();' 'let d = a.sub(&b).abs();' 'let d = -a;' 'let d = -(a + b);' 'let d = a - -b;' 'fn f(x: u8) -> (u8, u8) { (x - 1, x) }' 'let d = (a - b).neg();' 'let d = neg(a).sub(b);' 'let s = "-(b - a) and a.sub(b).neg()";' '// -(b - a)' 'let d = x.neg().sub(y);' > "$t/negok.diff"
ck; scan --diff-file "$t/negok.diff" >/dev/null || fail "clean subtraction forms flagged: $(scan --diff-file "$t/negok.diff")"

# the zero-unwrap rule applies to test code and test helpers too
for p in tests/props.rs tests/model/mod.rs benches/b.rs examples/e.rs; do
  mkdiff $p 'let v = r.unwrap();' > "$t/p.diff"
  ck; out="$(scan --diff-file "$t/p.diff")" && fail "unwrap not flagged in test/bench/example path: $p"
  ck; echo "$out" | grep -c "\[unwrap\]" > /dev/null || fail "wrong pattern label for unwrap in $p: $out"
done

mkdiff tests/props.rs 'let y = x.clone();' > "$t/p.diff"
ck; out="$(scan --diff-file "$t/p.diff")" && fail "clone not flagged in test path: tests/props.rs"
ck; echo "$out" | grep -c "\[clone\]" > /dev/null || fail "wrong pattern label for clone in tests/props.rs: $out"

# the assert family and the unwrap variants are banned in test code too
for p in tests/props.rs tests/model/mod.rs benches/b.rs examples/e.rs; do
  mkdiff $p 'assert!(x);' > "$t/p.diff"
  ck; out="$(scan --diff-file "$t/p.diff")" && fail "assert not flagged in test/bench/example path: $p"
  ck; echo "$out" | grep -c "\[assert\]" > /dev/null || fail "wrong pattern label for assert in $p: $out"
  mkdiff $p 'let v = r.unwrap_or_default();' > "$t/p.diff"
  ck; out="$(scan --diff-file "$t/p.diff")" && fail "unwrap_or_default not flagged in test/bench/example path: $p"
  ck; echo "$out" | grep -c "\[unwrap_variant\]" > /dev/null || fail "wrong pattern label for unwrap_or_default in $p: $out"
done

# a derive spanning several lines is flagged when a following line of the block names Clone or Copy
for names in 'Clone,' 'Copy,' 'Debug, Clone,'; do
  mkdiff tests/props.rs '#[derive(' '    Debug,' "    $names" ')]' 'struct A;' > "$t/p.diff"
  ck; out="$(scan --diff-file "$t/p.diff")" && fail "multi-line derive with $names not flagged"
  ck; echo "$out" | grep -c "\[derive_clone\]" > /dev/null || fail "wrong pattern label for multi-line derive with $names: $out"
done
mkdiff tests/props.rs '#[derive(' '    Debug,' '    PartialEq,' ')]' 'struct A;' 'let c = Clone;' > "$t/p.diff"
ck; scan --diff-file "$t/p.diff" >/dev/null || fail "multi-line derive without Clone/Copy flagged: $(scan --diff-file "$t/p.diff")"

mkdiff tests/props.rs 'let a: f64 = 1.0;' > "$t/p.diff"
ck; out="$(scan --diff-file "$t/p.diff")" && fail "f64 not flagged in test path: tests/props.rs"
ck; echo "$out" | grep -c "\[f64\]" > /dev/null || fail "wrong pattern label for f64 in tests/props.rs: $out"

mkdiff tests/props.rs 'panic!("x");' > "$t/p.diff"
ck; out="$(scan --diff-file "$t/p.diff")" && fail "panic not flagged in test path: tests/props.rs"
ck; echo "$out" | grep -c "\[panic\]" > /dev/null || fail "wrong pattern label for panic in tests/props.rs: $out"

# not flagged: comment, string, removed line, context line, deny attr, unwrap_or, unwrap_or_else, expect_err_x, assertion helpers without the macro bang
mkdiff $LIB '// x.clone() unsafe unwrap() f64' '/// doc: panic! #[allow]' 'let s = "x.clone() unwrap() f32";' '#![forbid(unsafe_code)]' 'let v = r.unwrap_or(0);' 'let v = r.expect_err_x;' 'let v = r.unwrap_or_else(f);' 'let r = assert_eq_x(a, b);' 'let s = "assert!(x) unwrap_or_default()";' > "$t/ok.diff"
ck; scan --diff-file "$t/ok.diff" >/dev/null || fail "comment/string/benign lines flagged: $(scan --diff-file "$t/ok.diff")"
for p in README.md docs/CONTRACT.md scripts/gen.py tools/t.py; do
  mkdiff $p 'let y = x.clone(); unsafe {} r.unwrap(); let a: f64 = 1.0;' > "$t/p.diff"
  ck; scan --diff-file "$t/p.diff" >/dev/null || fail "non-lib path flagged: $p"
done
printf 'diff --git a/%s b/%s\n--- a/%s\n+++ b/%s\n@@ -1,2 +1,1 @@\n-let y = x.clone();\n context\n' $LIB $LIB $LIB $LIB > "$t/rm.diff"
ck; scan --diff-file "$t/rm.diff" >/dev/null || fail "removed line flagged"

# allowlist: file:hash<TAB>reason
mkdiff $LIB 'let y = x.clone();' > "$t/d.diff"
ck; key="$(scan --diff-file "$t/d.diff" | grep -o 'allow-key=[^ ]*' | head -1 | cut -d= -f2)"; [ -n "$key" ] || fail "no allow-key printed"
printf '%s\tproof + Jeff approval\n' "$key" > "$t/tools/rules_allowlist.txt"
ck; scan --diff-file "$t/d.diff" >/dev/null || fail "allowlisted site still flagged"
: > "$t/tools/rules_allowlist.txt"

# hidden-copy count: fn bodies under src/ that rebuild a value from the fields of another value of its type
hc="$t/hc"; mkdir -p "$hc/src" "$hc/tools"
hidden() { RULES_SCAN_ROOT="$hc" "$here/rules_scan.sh" --hidden-copies; }
flagged_hc() { # flagged_hc <label> <rust source>: exactly one site expected
  printf '%s\n' "$2" > "$hc/src/x.rs"
  ck; out="$(hidden)" && fail "hidden copy not flagged: $1"
  ck; echo "$out" | grep -c "hidden_copy_count 1\$" > /dev/null || fail "hidden copy count is not 1 for $1: $out"
}
clean_hc() { # clean_hc <label> <rust source>
  printf '%s\n' "$2" > "$hc/src/x.rs"
  ck; out="$(hidden)" || fail "false hidden-copy flag: $1: $out"
}
flagged_hc "replicate by name" 'impl M { fn replicate(&self) -> M { M::zero() } }'
flagged_hc "copy by name" 'fn copy_of(x: &M) -> M { M::zero() }'
flagged_hc "dup by name" 'fn dup(x: &M) -> M { M::zero() }'
flagged_hc "old replicate body, renamed (pattern-bound fields)" 'impl M { fn again(&self) -> M { match &self.inner { Inner::Infinity => M::infinity(), Inner::Finite { exp, mant } => M::from_canonical(*exp, *mant) } } }'
flagged_hc "field reads in a call" 'fn again(x: &M) -> M { M::from_canonical(x.exp, x.mant) }'
flagged_hc "accessor reads in a call" 'fn again(x: &M) -> Option<M> { M::finite(x.mant()?, x.exp()?).ok() }'
flagged_hc "struct literal from another value" 'fn again(x: &P) -> P { Self { a: x.a, b: x.b } }'
flagged_hc "struct literal from pattern-bound names" 'fn again(x: &P) -> P { let P { a, b } = x; P { a: *a, b: *b } }'
clean_hc "computed args" 'fn next(x: &M) -> M { M::from_canonical(x.exp, x.mant + 1) }'
clean_hc "one field" 'fn up(x: &M) -> M { M::from_canonical(x.exp, TOP) }'
clean_hc "plain params" 'fn successor(exp: i64, mant: u64) -> M { M::from_canonical(exp, mant) }'
clean_hc "pattern-bound but modified" 'fn f(x: &M) -> M { match &x.inner { Inner::Finite { exp, mant } => M::from_canonical(exp.max(0), *mant & 1), _ => M::zero() } }'
clean_hc "type definition and match pattern" 'struct P { a: u8, b: u8 }
fn f(x: &P) -> bool { matches!(x, P { a: 0, .. }) }
fn g(x: &Inner) -> Option<(i64, u64)> { match x { Inner::Finite { exp, mant } if *mant != 0 => Some((*exp, *mant)), _ => None } }'
clean_hc "comment and string" 'fn f() -> M { // M::from_canonical(x.exp, x.mant) and fn replicate
  let _s = "M::from_canonical(x.exp, x.mant)"; M::zero() }'
ck; : > "$hc/src/x.rs"; hidden >/dev/null || fail "empty tree flagged"
# the same rules apply under tests/ benches/ examples/ (no test exemption); one fixture file at a time
flagged_hc_at() { # flagged_hc_at <path> <label> <rust source>: exactly one site expected
  mkdir -p "$(dirname "$hc/$1")"; printf '%s\n' "$3" > "$hc/$1"
  ck; out="$(hidden)" && fail "hidden copy not flagged in $1: $2"
  ck; echo "$out" | grep -c "hidden_copy_count 1\$" > /dev/null || fail "hidden copy count is not 1 in $1 for $2: $out"
  ck; echo "$out" | grep -c "^$1:" > /dev/null || fail "hidden copy site is not reported under $1 for $2: $out"
  rm -f "$hc/$1"
}
clean_hc_at() { # clean_hc_at <path> <label> <rust source>
  mkdir -p "$(dirname "$hc/$1")"; printf '%s\n' "$3" > "$hc/$1"
  ck; out="$(hidden)" || fail "false hidden-copy flag in $1: $2: $out"
  rm -f "$hc/$1"
}
DUP='fn dup(x: &MagV) -> MagV { match x { MagV::Finite { exp, mant } => MagV::Finite { exp: *exp, mant: *mant }, MagV::Infinity => MagV::Infinity } }'
for d in tests/model/direct.rs tests/mag/support.rs benches/b.rs examples/e.rs; do
  flagged_hc_at "$d" "dup-style helper by name" "$DUP"
  flagged_hc_at "$d" "same helper renamed (pattern-bound fields)" 'fn again(x: &MagV) -> MagV { match x { MagV::Finite { exp, mant } => MagV::Finite { exp: *exp, mant: *mant }, MagV::Infinity => MagV::Infinity } }'
done
flagged_hc_at tests/model/direct.rs "field reads in a call" 'fn again(x: &MagV) -> MagV { MagV::new(x.exp, x.mant) }'
flagged_hc_at tests/model/direct.rs "same type behind Option/Result" 'fn again(x: &MagV) -> Option<MagV> { MagV::mk(x.exp, x.mant).ok() }'
flagged_hc_at tests/model/direct.rs "copy of a model value through Result with a pattern" 'fn again(m: &MagV) -> Result<MagV, String> { match m { MagV::Finite { exp, mant } => Ok(MagV::from_parts(*exp, *mant)), MagV::Infinity => Ok(MagV::Infinity) } }'
flagged_hc_at tests/mag/support.rs "a converter still counts when it is named dup" 'fn dup(m: &MagV) -> Result<Mag, String> { Mag::zero() }'
flagged_hc_at tests/model/direct.rs "a method on the type (self receiver) is not a converter" 'impl MagV { fn again(&self) -> Mag { match self { MagV::Finite { exp, mant } => Mag::from_parts(*exp, *mant), MagV::Infinity => Mag::infinity() } } }'
clean_hc_at tests/mag/support.rs "model to library converter (other types in and out)" 'pub fn from_magv(m: &MagV) -> Result<Mag, String> { match m { MagV::Finite { exp, mant } => Mag::finite(*mant, *exp).map_err(|e| format!("{:?}", e)), MagV::Infinity => Ok(Mag::infinity()) } }'
clean_hc_at tests/mag/support.rs "library to model converter (to_magv)" 'pub fn to_magv(m: &Mag) -> Result<MagV, String> { match m { Mag::Finite { exp, mant } => Ok(MagV::Finite { exp: *exp, mant: *mant }), Mag::Infinity => Ok(MagV::Infinity) } }'
clean_hc_at tests/mag/support.rs "model to dyadic converter (dy_of)" 'pub fn dy_of(m: &MagV) -> Result<Dy, String> { match m { MagV::Finite { exp, mant } => Ok(Dy::from_mag(*exp, *mant)), MagV::Infinity => Err(String::new()) } }'
clean_hc_at tests/model/dyadic.rs "model to dyadic converter" 'pub fn from_magv(m: &MagV) -> Option<Dy> { match m { MagV::Finite { exp, mant } => Some(Dy::from_mag(*exp, *mant)), MagV::Infinity => None } }'
clean_hc_at tests/model/direct.rs "from-parts rebuild with computed arguments" 'fn f(a: i128, m: u128) -> Option<MagV> { Some(MagV::Finite { exp: i64::try_from(a).ok()?, mant: u64::try_from(m).ok()? }) }'
clean_hc_at tests/model/direct.rs "comment and string" 'fn f() -> MagV { // MagV::new(x.exp, x.mant) and fn dup
  let _s = "fn dup(x: &MagV)"; zero() }'
ck; hidden >/dev/null || fail "tests fixtures left a flag behind: $(hidden)"

# Rusty's six converter-exemption evasions (green3 review A2 A3 A4 A5 A25 A31): a field-by-field copy of a Mag/MagV must be flagged
# even when the type is renamed, aliased or wrapped. Soft-fail so every plant is reported in one run.
SOFT=0
soft() { echo "SELFTEST FAIL: $*"; SOFT=$((SOFT + 1)); }
plant_hc() { # plant_hc <id> <path> <rust source>: exactly one hidden-copy site expected at <path>
  mkdir -p "$(dirname "$hc/$2")"; printf '%s\n' "$3" > "$hc/$2"
  out="$(hidden)"; rc=$?
  if [ $rc -eq 0 ] || ! echo "$out" | grep -c "hidden_copy_count 1\$" > /dev/null || ! echo "$out" | grep -c "^$2:" > /dev/null; then
    soft "plant $1 not flagged by the hidden-copy count: $(echo "$out" | tr '\n' ' ')"
  fi
  rm -f "$hc/$2"
}
INNER='match &m.inner { Inner::Finite { exp, mant } => Mag::from_canonical(*exp, *mant), Inner::Infinity => Mag::infinity() }'
plant_hc A2 src/mag/zz.rs "type M2 = Mag; pub fn twin(m: &M2) -> Mag { $INNER }"
plant_hc A3 src/mag/zz.rs "type M2 = Mag; pub fn twin(m: &Mag) -> M2 { $INNER }"
plant_hc A4 src/mag/zz.rs "use super::Mag as M; pub fn twin(m: &Mag) -> M { match &m.inner { Inner::Finite { exp, mant } => M::from_canonical(*exp, *mant), Inner::Infinity => M::infinity() } }"
plant_hc A5 src/mag/zz.rs "pub struct W(pub Mag); pub fn twin(w: &W) -> Mag { match &w.0.inner { Inner::Finite { exp, mant } => Mag::from_canonical(*exp, *mant), Inner::Infinity => Mag::infinity() } }"
plant_hc A25 tests/model/zz.rs "use super::dyadic::MagV; type V = MagV; pub fn again(x: &V) -> Option<MagV> { Some(MagV::Finite { exp: x.exp()?, mant: x.mant()? }) }"
plant_hc A31 src/mag/zz.rs "pub struct Pair(pub Mag, pub Mag); pub fn twin(m: &Mag) -> Pair { match &m.inner { Inner::Finite { exp, mant } => Pair(Mag::from_canonical(*exp, *mant), Mag::zero()), Inner::Infinity => Pair(Mag::infinity(), Mag::zero()) } }"
# no exemption at all under src/, even for a real converter; under tests/ only when every named type is Mag, MagV or Dy
plant_hc S1 src/mag/zz.rs 'pub fn from_magv(m: &MagV) -> Result<Mag, String> { match m { MagV::Finite { exp, mant } => Mag::finite(*mant, *exp).map_err(|e| format!("{:?}", e)), MagV::Infinity => Ok(Mag::infinity()) } }'
plant_hc T1 tests/mag/zz.rs 'pub fn from_magv(m: &MagV) -> Result<Mag, BoundError> { match m { MagV::Finite { exp, mant } => Mag::finite(*mant, *exp), MagV::Infinity => Ok(Mag::infinity()) } }'
plant_hc T2 tests/mag/zz.rs 'pub fn from_w(m: &W) -> Result<Mag, String> { match m { W::Finite { exp, mant } => Mag::finite(*mant, *exp).map_err(|e| format!("{:?}", e)), W::Infinity => Ok(Mag::infinity()) } }'
plant_hc T3 tests/mag/zz.rs 'pub fn conv(m: &MagV) -> Result<Pair, String> { match m { MagV::Finite { exp, mant } => Ok(Pair::new(*mant, *exp)), MagV::Infinity => Ok(Pair::inf()) } }'
# the aliasing and renaming of Mag/MagV/Dy is banned outright (rules scan label, no baseline, any hit fails)
while IFS='|' read -r name line; do
  mkdiff src/mag/zz.rs "$line" > "$t/d.diff"
  out="$(scan --diff-file "$t/d.diff")" && { soft "not flagged: $name: $line"; continue; }
  echo "$out" | grep -c "\[$name\]" > /dev/null || soft "wrong pattern label for $name: $out"
done <<'CASES'
type_alias|type M2 = Mag;
type_alias|pub type V = MagV;
type_alias|type D = Dy;
type_alias|pub(crate) type Pairs = Vec<(Mag, Mag)>;
type_alias|type R = Result<Mag, String>;
type_alias|type Output = Mag;
type_alias|type V<'a> = &'a MagV;
type_rename|use super::Mag as M;
type_rename|use crate::model::dyadic::{Dy as D, MagV};
type_rename|use crate::model::dyadic::{MagV, Dy as D};
type_rename|use zenith_bound::mag::Mag as Magnitude;
CASES
mkdiff src/mag/zz.rs 'type Pair = (u8, u8);' 'use std::fmt::Result as FmtResult;' 'use super::Mag;' 'let m = Mag::zero(); // type M2 = Mag; use x::Mag as M' 'let s = "type M2 = Mag; Mag as M";' 'type Word = u64;' 'impl Magnitude for X {}' > "$t/ok.diff"
scan --diff-file "$t/ok.diff" >/dev/null || soft "alias-free lines flagged: $(scan --diff-file "$t/ok.diff")"
# rustfmt-wrapped aliases and renames (Rusty's 7af7cc2 review G1): the labels see the joined added lines, so a wrap cannot hide them
wrapped() { # wrapped <label> <added line>...: the lines form ONE statement; the label must be reported
  local name="$1"; shift
  mkdiff src/mag/zz.rs "$@" > "$t/d.diff"
  out="$(scan --diff-file "$t/d.diff")" && { soft "wrapped form not flagged: $name: $*"; return; }
  echo "$out" | grep -c "\[$name\]" > /dev/null || soft "wrong pattern label for wrapped $name: $out"
}
wrapped type_alias 'type X =' '    Option<' '        Mag>;'
wrapped type_rename 'use super::{Mag' '    as M3};'
wrapped type_rename 'pub use crate::model::dyadic::{' '    MagV,' '    Dy' '        as D,' '};'
wrapped type_rename 'use super::Mag as' '    M;'
wrapped type_alias 'pub type Pairs<' '    T,' '> = Result<' '    Mag,' '    T,' '> where' '    T: Sized;'
wrapped type_alias 'pub type V<'"'"'a> =' '    &'"'"'a' '    MagV;'
wrapped type_alias 'type D' '    = Dy;'
mkdiff src/mag/zz.rs 'type Pair =' '    (u8,' '        u8);' 'use std::fmt::{' '    Result as FmtResult,' '};' 'type Word = u64;' 'let m = Mag::zero();' 'use super::{' '    Mag,' '    Inner,' '};' > "$t/ok.diff"
scan --diff-file "$t/ok.diff" >/dev/null || soft "wrapped alias-free lines flagged: $(scan --diff-file "$t/ok.diff")"
[ "$SOFT" -eq 0 ] || { echo "SELFTEST: $SOFT converter-exemption plant check(s) failed"; exit 1; }

# git modes in a throwaway repo
g="$t/repo"; mkdir -p "$g/src" "$g/tools"; cd "$g" || exit 1
git init -q . && git config user.email a@b && git config user.name t
cp "$t/tools/rules_allowlist.txt" tools/; echo 0 > tools/clone_count_baseline.txt; echo 0 > tools/ref_sum_count_baseline.txt; echo 0 > tools/hidden_copy_count_baseline.txt; echo 0 > tools/sanctioned_rebuild_count_baseline.txt
printf 'fn a() {}\n' > src/x.rs; git add -A; git commit -qm base; base=$(git rev-parse HEAD)
export RULES_SCAN_ROOT="$g"; unset RULES_SCAN_SKIP_COUNT
ck; scan "$base..HEAD" >/dev/null || fail "empty range flagged"
printf 'fn a() {}\nfn b(x: &S) { let y = x.clone(); }\n' > src/x.rs
ck; scan --working >/dev/null && fail "working diff clone not flagged"
ck; git add -A; scan --staged >/dev/null && fail "staged diff clone not flagged"
ck; scan --working >/dev/null && fail "working mode (git diff HEAD) must still see staged additions"
ck; git commit -qm bad; scan "$base..HEAD" >/dev/null && fail "range clone not flagged"
# a #[cfg(test)] module gets no exemption
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(x: &S) { let y = x.clone(); let _ = r.unwrap(); }\n}\n' > src/x.rs
ck; git add -A; git commit -qm t; scan "$base..HEAD" >/dev/null && fail "cfg(test) additions not flagged"

# clone-count baseline: may only go down
printf 'fn a(x: &S) { let y = x.clone(); }\n' > src/x.rs; git add -A; git commit -qm c
echo 1 > tools/clone_count_baseline.txt
ck; scan --count-only >/dev/null || fail "count == baseline should pass"
ck; echo 0 > tools/clone_count_baseline.txt; scan --count-only >/dev/null && fail "count above baseline passed"
printf 'fn a(x: &S) { let y = x.clone(); let z = it.copied(); }\n' > src/x.rs; echo 2 > tools/clone_count_baseline.txt
ck; { scan --count-only || true; } | grep -c "clone count 2 == baseline" > /dev/null || fail "a .copied( use must count with the clone family: $(scan --count-only)"
printf 'fn a(x: &S) { let y = x.clone(); }\n' > src/x.rs
ck; echo 5 > tools/clone_count_baseline.txt; { scan --count-only || true; } | grep -ci "lower" > /dev/null || fail "no hint to lower baseline"
ck; scan --update-baseline >/dev/null || fail "lowering baseline failed"; [ "$(tr -d '[:space:]' < tools/clone_count_baseline.txt)" = 1 ] || fail "baseline not lowered"
ck; echo 0 > tools/clone_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising baseline allowed"; [ "$(tr -d '[:space:]' < tools/clone_count_baseline.txt)" = 0 ] || fail "baseline raised"

# add_ref/sub_ref/mul_ref/div_ref use-count baseline: counts every use incl. tests (not definitions, comments, strings), may only go down
cat > src/x.rs <<'EOF'
impl N { pub fn add_ref(a: &N, b: &N) -> N { a.go(b) } pub fn sub_ref(a: &N, b: &N) -> N { a.go(b) } pub fn mul_ref(a: &N, b: &N) -> N { a.go(b) } pub fn div_ref(a: &N, b: &N) -> N { a.go(b) } }
fn u(a: &N, b: &N) { let s = N::add_ref(a, b); let d = N::sub_ref(a, b); let m = N::mul_ref(a, b); let q = N::div_ref(a, b); }
// N::add_ref(a, b) in a comment
fn s() { let x = "N::sub_ref(a, b)"; }
#[cfg(test)]
mod tests { fn t(a: &N, b: &N) { let s = N::add_ref(a, b); } }
EOF
git add -A; git commit -qm refs
echo 5 > tools/ref_sum_count_baseline.txt; echo 7 > tools/clone_count_baseline.txt
ck; { scan --count-only || true; } | grep -c "use count 5 == baseline" > /dev/null || fail "ref count should be 5 (definitions, comments, strings excluded; the test-module use counts): $(scan --count-only)"
ck; echo 4 > tools/ref_sum_count_baseline.txt; { scan --count-only || true; } | grep -c "ABOVE" > /dev/null || fail "ref count above baseline not reported"
ck; scan --count-only >/dev/null && fail "ref count above baseline passed"
ck; echo 9 > tools/ref_sum_count_baseline.txt; { scan --count-only || true; } | grep -ci "lower the baseline" > /dev/null || fail "no hint to lower ref baseline"
ck; scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/ref_sum_count_baseline.txt)" = 5 ] || fail "ref baseline not lowered"
ck; echo 0 > tools/ref_sum_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising ref baseline allowed"; [ "$(tr -d '[:space:]' < tools/ref_sum_count_baseline.txt)" = 0 ] || fail "ref baseline raised"
ck; rm -f tools/ref_sum_count_baseline.txt; scan --count-only >/dev/null && fail "missing ref baseline passed"
echo 0 > tools/ref_sum_count_baseline.txt

# hidden-copy baseline: counts fn bodies, may only go down
printf 'fn replicate(x: &S) -> S { S::new() }\nfn again(x: &S) -> S { S::new(x.a, x.b) }\n' > src/x.rs; git add -A; git commit -qm hc
echo 2 > tools/hidden_copy_count_baseline.txt
ck; { scan --count-only || true; } | grep -c "hidden-copy count 2 == baseline" > /dev/null || fail "hidden-copy count should be 2: $(scan --count-only)"
ck; echo 1 > tools/hidden_copy_count_baseline.txt; { scan --count-only || true; } | grep -c "hidden-copy count 2 is ABOVE" > /dev/null || fail "hidden-copy count above baseline not reported"
ck; scan --count-only >/dev/null && fail "hidden-copy count above baseline passed"
ck; echo 9 > tools/hidden_copy_count_baseline.txt; { scan --count-only || true; } | grep -ci "hidden-copy count 2 is below baseline" > /dev/null || fail "no hint to lower hidden-copy baseline"
ck; scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/hidden_copy_count_baseline.txt)" = 2 ] || fail "hidden-copy baseline not lowered"
ck; echo 0 > tools/hidden_copy_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising hidden-copy baseline allowed"; [ "$(tr -d '[:space:]' < tools/hidden_copy_count_baseline.txt)" = 0 ] || fail "hidden-copy baseline raised"
# the hidden-copy count includes tests/: a dup-style helper there raises the count above a 0 baseline
printf 'fn a() {}\n' > src/x.rs; mkdir -p tests; printf 'fn dup(x: &MagV) -> MagV { match x { MagV::Finite { exp, mant } => MagV::Finite { exp: *exp, mant: *mant }, MagV::Infinity => MagV::Infinity } }\n' > tests/m.rs; git add -A; git commit -qm hct
echo 0 > tools/hidden_copy_count_baseline.txt
ck; { scan --count-only || true; } | grep -c "hidden-copy count 1 is ABOVE" > /dev/null || fail "hidden copy in tests/ not counted: $(scan --count-only)"
ck; scan --count-only >/dev/null && fail "hidden copy in tests/ passed a 0 baseline"
ck; rm -f tools/hidden_copy_count_baseline.txt; scan --count-only >/dev/null && fail "missing hidden-copy baseline passed"

# sanctioned-rebuild baseline (O18-C): counts Mag::from_canonical(..) calls in src/ fed by a part of another value; may only go down
rb() { RULES_SCAN_ROOT="$g" "$here/rules_scan.sh" --sanctioned-rebuilds; }
rm -f tests/m.rs; echo 0 > tools/hidden_copy_count_baseline.txt; echo 0 > tools/clone_count_baseline.txt; echo 0 > tools/ref_sum_count_baseline.txt   # the other three counts are 0 and have baselines, so only the rebuild count can fail below
cat > src/x.rs <<'EOF'
impl Mag {
    fn from_canonical(exp: i64, mant: u64) -> Mag { Mag { exp, mant } }
    fn finite(exp: i64, mant: u64) -> Mag { Mag::from_canonical(exp, mant << 1) }
    fn min_pos() -> Mag { Mag::from_canonical(i64::MIN, TOP) }
    fn pack(e: i64, mant: u64) -> Mag { Self::from_canonical(e, mant) }
    fn next(e: i64, next: u64) -> Mag { Mag::from_canonical(e.checked_add(1)?, next) }
    // Mag::from_canonical(e, m.get()) in a comment
    fn s() -> Mag { let _x = "Mag::from_canonical(e, m.get())"; Mag::zero() }
    pub fn add_up(&self, other: &Mag) -> Mag {
        match (self.nonzero(), other.nonzero()) {
            (Some((e, m)), None) | (None, Some((e, m))) => Mag::from_canonical(e, m.get()),
            _ => Mag::zero(),
        }
    }
    pub fn sub_lower(&self, other: &Mag) -> Mag {
        match self.nonzero() {
            Some((e, m)) => Mag::from_canonical(
                e,
                m.get(),
            ),
            None => Mag::zero(),
        }
    }
}
EOF
git add -A; git commit -qm rb2
echo 2 > tools/sanctioned_rebuild_count_baseline.txt
ck; out="$(rb)"; echo "$out" | grep -c "sanctioned_rebuild_count 2\$" > /dev/null || fail "sanctioned rebuilds should be 2 (definition, computed and constant args, comment, string not counted; multi-line call counted): $out"
ck; { scan --count-only || true; } | grep -c "sanctioned-rebuild count 2 == baseline" > /dev/null || fail "sanctioned-rebuild count should equal baseline 2: $(scan --count-only)"
cat >> src/x.rs <<'EOF'
impl Mag { pub fn third(&self) -> Mag { match self.nonzero() { Some((e, m)) => Mag::from_canonical(e, m.get()), None => Mag::zero() } } }
EOF
git add -A; git commit -qm rb3
ck; { scan --count-only || true; } | grep -c "sanctioned-rebuild count 3 is ABOVE the baseline 2" > /dev/null || fail "a third rebuild not reported: $(scan --count-only)"
ck; scan --count-only >/dev/null && fail "a third rebuild passed the baseline of 2"
printf 'impl Mag { fn again(x: &Mag) -> Mag { Mag::from_canonical(x.exp, x.mant) } fn again2(x: &Mag) -> Mag { Mag::from_canonical(x.exp(), x.mant()) } }\n' > src/x.rs
ck; rb | grep -c "sanctioned_rebuild_count 2\$" > /dev/null || fail "field and accessor reads must count: $(rb)"
ck; echo 9 > tools/sanctioned_rebuild_count_baseline.txt; { scan --count-only || true; } | grep -ci "sanctioned-rebuild count 2 is below baseline 9" > /dev/null || fail "no hint to lower rebuild baseline"
ck; scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/sanctioned_rebuild_count_baseline.txt)" = 2 ] || fail "rebuild baseline not lowered"
ck; echo 0 > tools/sanctioned_rebuild_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising rebuild baseline allowed"; [ "$(tr -d '[:space:]' < tools/sanctioned_rebuild_count_baseline.txt)" = 0 ] || fail "rebuild baseline raised"
ck; rm -f tools/sanctioned_rebuild_count_baseline.txt; scan --count-only >/dev/null && fail "missing rebuild baseline passed"

# alias/rename ban on the whole tree (baseline-free, wrapped forms included): any hit fails --count-only and the gate's default run
echo 0 > tools/sanctioned_rebuild_count_baseline.txt; echo 0 > tools/hidden_copy_count_baseline.txt; rm -rf tests
printf 'use super::Mag;\ntype Word = u64;\n// type M2 = Mag; and Mag as M\nfn s() { let _x = "type M2 = Mag;"; }\n' > src/x.rs
ck; scan --count-only >/dev/null || fail "alias-free tree flagged: $(scan --count-only)"
printf 'type X =\n    Option<\n        Mag>;\n' > src/x.rs
ck; out="$(scan --count-only)" && fail "wrapped type alias in the tree passed --count-only"
ck; echo "$out" | grep -c "src/x.rs:1 \[type_alias\]" > /dev/null || fail "wrapped type alias not reported with its first line: $out"
printf 'use super::{Mag\n    as M3};\n' > src/x.rs
ck; out="$(scan --count-only)" && fail "wrapped rename in the tree passed --count-only"
ck; echo "$out" | grep -c "src/x.rs:1 \[type_rename\]" > /dev/null || fail "wrapped rename not reported: $out"
mkdir -p tests; printf 'pub use crate::model::dyadic::{\n    MagV,\n    Dy\n        as D,\n};\n' > tests/y.rs; printf 'fn a() {}\n' > src/x.rs
ck; out="$(scan --count-only)" && fail "wrapped pub use rename under tests/ passed --count-only"
ck; echo "$out" | grep -c "tests/y.rs:3 \[type_rename\]" > /dev/null || fail "rename under tests/ not reported at its line: $out"
rm -rf tests; printf 'fn a() {}\n' > src/x.rs
ck; scan --count-only >/dev/null || fail "tree not clean after the alias fixtures: $(scan --count-only)"

# API-freeze unit 4: no public signature in src/ names zf's i32 `Exponent` or `i32` (zb signatures use i64); baseline-free
printf 'pub fn e(&self) -> Option<i64> { None }\npub(crate) fn p(x: i32) {}\nfn q(x: Exponent) {}\n// pub fn c(x: i32)\n' > src/x.rs
ck; scan --count-only >/dev/null || fail "i64 / private i32 flagged as narrow exponent: $(scan --count-only)"
printf 'impl A {\n    pub fn exp(\n        &self,\n    ) -> Option<Exponent> { None }\n}\n' > src/x.rs
ck; out="$(scan --count-only)" && fail "pub fn returning Exponent passed"
ck; echo "$out" | grep -q "src/x.rs:2 \[pub_sig_narrow_exponent\]" || fail "wrapped Exponent signature not reported at its line: $out"
ck; printf 'pub fn s(e: i32) {}\n' > src/x.rs; scan --count-only >/dev/null && fail "pub fn taking i32 passed"
ck; printf 'pub struct S { pub e: i32 }\n' > src/x.rs; scan --count-only >/dev/null && fail "pub field of type i32 passed"
ck; printf 'pub const E: zenith_float_num::Exponent = 0;\npub type T = i32;\n' > src/x.rs; scan --count-only >/dev/null && fail "pub const/type naming Exponent/i32 passed"
printf 'fn a() {}\n' > src/x.rs
ck; scan --count-only >/dev/null || fail "tree not clean after the narrow-exponent fixtures: $(scan --count-only)"

# API-freeze unit 5: every prec in src/ uses the `Precision` alias, never a concrete integer width (params, fields, fn prec returns)
printf 'pub type Precision = usize;\npub struct B { prec: Precision }\npub fn f(prec: Precision) -> Precision { prec }\npub trait E { fn prec(&self) -> Precision; }\n// fn g(prec: u32)\nfn h(n: u32) {}\n' > src/x.rs
ck; scan --count-only >/dev/null || fail "Precision alias flagged: $(scan --count-only)"
for bad in 'pub fn f(prec: u32) {}' 'pub fn f(x: u8, mid_prec: usize) {}' 'struct B { prec: u64 }' 'pub trait E { fn prec(&self) -> u32; }' 'pub fn prec(\n    &self,\n) -> i64 { 0 }' 'fn f(prec: core::primitive::u128) {}'; do
  printf "$bad\n" > src/x.rs
  ck; out="$(scan --count-only)" && fail "concrete prec width passed: $bad"
  ck; echo "$out" | grep -q "src/x.rs:1 \[prec_concrete_width\]" || fail "concrete prec width not reported at line 1: $bad: $out"
done
printf 'fn a() {}\n' > src/x.rs
ck; scan --count-only >/dev/null || fail "tree not clean after the prec fixtures: $(scan --count-only)"
# no vacuous pass: --working and --staged outside a git repo (git diff fails) must exit nonzero, not report "0 flagged"
nogit=$(mktemp -d "$(zb_scratch)"/rules-nogit.XXXXXX)
for mode in --working --staged; do
  ck; out="$(RULES_SCAN_ROOT="$nogit" python3 "$here/rules_scan.py" $mode 2>&1)" && { rm -rf "$nogit"; fail "$mode passed although git diff failed: $out"; }
  ck; echo "$out" | grep -c "git diff .*failed" > /dev/null || { rm -rf "$nogit"; fail "$mode: no git failure message: $out"; }
done
rm -rf "$nogit"
st_done rules_scan 284   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "rules_scan selftest: ok"
