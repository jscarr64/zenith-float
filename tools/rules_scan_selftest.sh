#!/usr/bin/env bash
# Self-test for tools/rules_scan.sh: fixture diffs prove every banned pattern is flagged on ADDED lines,
# that removed/context lines, comments, strings, test code and non-lib paths are not, that the
# allowlist works, that staged/working/range modes work, that the clone-count, helper-copy and ref baselines only go down,
# and that the constant-cache exception allowlist works (unlisted clone fails, listed passes, stale entry fails, count only down).
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d)"; trap 'rm -rf "$t"' EXIT
scan() { "$here/rules_scan.sh" "$@"; }
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
mkdiff() { # mkdiff <path> <added line...>  -> fixture diff on stdout
  local p="$1"; shift
  printf 'diff --git a/%s b/%s\n--- a/%s\n+++ b/%s\n@@ -1,1 +1,%d @@\n context\n' "$p" "$p" "$p" "$p" "$(( $# + 1 ))"
  for l in "$@"; do printf '+%s\n' "$l"; done
}
LIB=zenith-float-num/src/foo.rs
export RULES_SCAN_ROOT="$t"; mkdir -p "$t/tools"; : > "$t/tools/rules_allowlist.txt"; echo 0 > "$t/tools/clone_count_baseline.txt"; echo 0 > "$t/tools/ref_sum_count_baseline.txt"; echo 0 > "$t/tools/mul_ref_count_baseline.txt"; echo 0 > "$t/tools/helper_copy_count_baseline.txt"; echo 0 > "$t/tools/exception_count_baseline.txt"
export RULES_SCAN_SKIP_COUNT=1   # fixture diffs have no tree to count

mkdiff $LIB 'let a = 1 + 2;' 'let s = a.max(b);' > "$t/clean.diff"
scan --diff-file "$t/clean.diff" >/dev/null || fail "clean diff flagged"

while IFS='|' read -r name line; do
  mkdiff $LIB "$line" > "$t/d.diff"
  out="$(scan --diff-file "$t/d.diff")" && fail "not flagged: $name"
  echo "$out" | grep -q "\[$name\]" || fail "wrong pattern label for $name: $out"
done <<'CASES'
clone|let y = x.clone();
to_vec|let y = s.to_vec();
cloned|let y = it.cloned();
to_owned|let y = s.to_owned();
unsafe|unsafe { core::hint::unreachable_unchecked() }
allow|#[allow(dead_code)]
allow|#![allow(clippy::all)]
unwrap|let v = r.unwrap();
expect|let v = r.expect("boom");
panic|panic!("no");
f32|let a: f32 = 1.0;
f64|let a = 1.0_f64; let b: f64 = 2.0;
CASES

# not flagged: comment, string, removed line, context line, deny attr, unwrap_or, unsafe_code, non-library paths (tools/, top-level tests/benches/examples, README); cfg(test) not visible in diff-file mode so paths only
mkdiff $LIB '// x.clone() unsafe unwrap() f64' '/// doc: panic! #[allow]' 'let s = "x.clone() unwrap() f32";' '#![forbid(unsafe_code)]' 'let v = r.unwrap_or(0);' 'let v = r.expect_err_x;' > "$t/ok.diff"
scan --diff-file "$t/ok.diff" >/dev/null || fail "comment/string/benign lines flagged: $(scan --diff-file "$t/ok.diff")"
for p in zenith-float-num/tests/x.rs benches/b.rs examples/e.rs tools/t.rs README.md; do
  mkdiff $p 'let y = x.clone(); unsafe {} r.unwrap(); let a: f64 = 1.0;' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null || fail "non-lib path flagged: $p"
done
printf 'diff --git a/%s b/%s\n--- a/%s\n+++ b/%s\n@@ -1,2 +1,1 @@\n-let y = x.clone();\n context\n' $LIB $LIB $LIB $LIB > "$t/rm.diff"
scan --diff-file "$t/rm.diff" >/dev/null || fail "removed line flagged"

# test code under the library dirs is counted now: tests.rs and tests/ benches/ examples/ components are NOT exempt
for p in zenith-float-num/src/ops/tests.rs zenith-float-num/src/tests/x.rs zenith-float-num/src/benches/b.rs; do
  mkdiff $p 'let y = x.clone();' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null && fail "test-code path no longer exempt, but not flagged: $p"
done

# f64 is only a calculation-path rule: compare crate is lib for clone but not for f64
mkdiff zenith-float-compare/src/a.rs 'let a: f64 = 1.0;' > "$t/c.diff"; scan --diff-file "$t/c.diff" >/dev/null || fail "f64 flagged in compare (not a calc path)"
mkdiff zenith-float-compare/src/a.rs 'let a = b.clone();' > "$t/c.diff"; scan --diff-file "$t/c.diff" >/dev/null && fail "clone in compare not flagged"

# allowlist: file:hash<TAB>reason
mkdiff $LIB 'let y = x.clone();' > "$t/d.diff"
key="$(scan --diff-file "$t/d.diff" | grep -o 'allow-key=[^ ]*' | head -1 | cut -d= -f2)"; [ -n "$key" ] || fail "no allow-key printed"
printf '%s\tproof + Jeff approval\n' "$key" > "$t/tools/rules_allowlist.txt"
scan --diff-file "$t/d.diff" >/dev/null || fail "allowlisted site still flagged"
: > "$t/tools/rules_allowlist.txt"

# git modes in a throwaway repo
g="$t/repo"; mkdir -p "$g/zenith-float-num/src" "$g/tools"; cd "$g" || exit 1
git init -q . && git config user.email a@b && git config user.name t
cp "$t/tools/rules_allowlist.txt" tools/; echo 0 > tools/clone_count_baseline.txt; echo 0 > tools/ref_sum_count_baseline.txt; echo 0 > tools/mul_ref_count_baseline.txt; echo 0 > tools/helper_copy_count_baseline.txt; echo 0 > tools/exception_count_baseline.txt
printf 'fn a() {}\n' > zenith-float-num/src/x.rs; git add -A; git commit -qm base; base=$(git rev-parse HEAD)
export RULES_SCAN_ROOT="$g"; unset RULES_SCAN_SKIP_COUNT
scan "$base..HEAD" >/dev/null || fail "empty range flagged"
printf 'fn a() {}\nfn b(x: &S) { let y = x.clone(); }\n' > zenith-float-num/src/x.rs
scan --working >/dev/null && fail "working diff clone not flagged"
git add -A; scan --staged >/dev/null && fail "staged diff clone not flagged"
scan --working >/dev/null && fail "working mode (git diff HEAD) must still see staged additions"
git commit -qm bad; scan "$base..HEAD" >/dev/null && fail "range clone not flagged"
# cfg(test) module in the tip file is exempt from the diff scan's non-clone rules (the clone COUNT does count it: see below)
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(x: &S) { let y = x.clone(); let _ = r.unwrap(); }\n}\n' > zenith-float-num/src/x.rs
git add -A; git commit -qm t; RULES_SCAN_SKIP_COUNT=1 scan "$base..HEAD" >/dev/null || fail "cfg(test) additions flagged: $(RULES_SCAN_SKIP_COUNT=1 scan "$base..HEAD")"
scan "$base..HEAD" >/dev/null && fail "clone in a cfg(test) module not counted by the range scan's clone count (baseline 0)"

# clone-count baseline: may only go down
printf 'fn a(x: &S) { let y = x.clone(); }\n' > zenith-float-num/src/x.rs; git add -A; git commit -qm c
echo 1 > tools/clone_count_baseline.txt
scan --count-only >/dev/null || fail "count == baseline should pass"
echo 0 > tools/clone_count_baseline.txt; scan --count-only >/dev/null && fail "count above baseline passed"
echo 5 > tools/clone_count_baseline.txt; scan --count-only | grep -qi "lower" || fail "no hint to lower baseline"
scan --update-baseline >/dev/null || fail "lowering baseline failed"; [ "$(tr -d '[:space:]' < tools/clone_count_baseline.txt)" = 1 ] || fail "baseline not lowered"
echo 0 > tools/clone_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising baseline allowed"; [ "$(tr -d '[:space:]' < tools/clone_count_baseline.txt)" = 0 ] || fail "baseline raised"

# add_ref/sub_ref use-count baseline: counts non-test library uses (not definitions, comments, strings, tests), may only go down
cat > zenith-float-num/src/x.rs <<'EOF'
impl N { pub fn add_ref(a: &N, b: &N) -> N { a.go(b) } pub fn sub_ref(a: &N, b: &N) -> N { a.go(b) } }
fn u(a: &N, b: &N) { let s = N::add_ref(a, b); let d = N::sub_ref(a, b); }
// N::add_ref(a, b) in a comment
fn s() { let x = "N::sub_ref(a, b)"; }
#[cfg(test)]
mod tests { fn t(a: &N, b: &N) { let s = N::add_ref(a, b); } }
EOF
git add -A; git commit -qm refs
echo 2 > tools/ref_sum_count_baseline.txt; echo 5 > tools/clone_count_baseline.txt
scan --count-only | grep -q "add_ref/sub_ref use count 2 == baseline" || fail "ref count should be 2 (definitions, comments, strings, tests excluded): $(scan --count-only)"
echo 1 > tools/ref_sum_count_baseline.txt; scan --count-only | grep -q "ABOVE" || fail "ref count above baseline not reported"
scan --count-only >/dev/null && fail "ref count above baseline passed"
echo 7 > tools/ref_sum_count_baseline.txt; scan --count-only | grep -qi "lower the baseline" || fail "no hint to lower ref baseline"
scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/ref_sum_count_baseline.txt)" = 2 ] || fail "ref baseline not lowered"
echo 0 > tools/ref_sum_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising ref baseline allowed"; [ "$(tr -d '[:space:]' < tools/ref_sum_count_baseline.txt)" = 0 ] || fail "ref baseline raised"
rm -f tools/ref_sum_count_baseline.txt; scan --count-only >/dev/null && fail "missing ref baseline passed"
# mul_ref/div_ref use-count baseline: its own file, same rules, independent of the add_ref/sub_ref count
cat > zenith-float-num/src/x.rs <<'EOF'
impl N { pub fn mul_ref(a: &N, b: &N) -> N { a.go(b) } pub fn div_ref(a: &N, b: &N) -> N { a.go(b) } }
fn u(a: &N, b: &N) { let s = N::mul_ref(a, b); let d = N::div_ref(a, b); let z = N::add_ref(a, b); }
// N::mul_ref(a, b) in a comment
fn s() { let x = "N::div_ref(a, b)"; }
impl N { pub fn sqr(&self, p: usize) -> N { self.go(self) } }
fn q(a: &N) { let s = a.sqr(1); let t = N::sqr(a, 2); }
#[cfg(test)]
mod tests { fn t(a: &N, b: &N) { let s = N::mul_ref(a, b); } }
EOF
# sqr is a plain function (no copy): it is not a mul_ref use, so the two sqr uses above leave the count at 2
echo 1 > tools/ref_sum_count_baseline.txt; echo 2 > tools/mul_ref_count_baseline.txt
scan --count-only | grep -q "mul_ref/div_ref use count 2 == baseline" || fail "mul ref count should be 2: $(scan --count-only)"
scan --count-only | grep -q "add_ref/sub_ref use count 1 == baseline" || fail "add ref count should stay 1 beside mul_ref: $(scan --count-only)"
echo 1 > tools/mul_ref_count_baseline.txt; scan --count-only | grep -q "mul_ref/div_ref use count 2 is ABOVE" || fail "mul ref above baseline not reported"
scan --count-only >/dev/null && fail "mul ref count above baseline passed"
echo 7 > tools/mul_ref_count_baseline.txt; scan --count-only | grep -qi "mul_ref/div_ref use count 2 is below baseline" || fail "no hint to lower mul ref baseline"
scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/mul_ref_count_baseline.txt)" = 2 ] || fail "mul ref baseline not lowered"
echo 0 > tools/mul_ref_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising mul ref baseline allowed"; [ "$(tr -d '[:space:]' < tools/mul_ref_count_baseline.txt)" = 0 ] || fail "mul ref baseline raised"
rm -f tools/mul_ref_count_baseline.txt; scan --count-only >/dev/null && fail "missing mul ref baseline passed"
# ---- test code is counted by the clone count (cfg(test) items, tests.rs, test dirs) ----
echo 0 > tools/ref_sum_count_baseline.txt; echo 0 > tools/mul_ref_count_baseline.txt
rm -f zenith-float-num/src/x.rs; mkdir -p zenith-float-num/src/ops zenith-float-num/tests nostd-tests/src
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(x: &S) { let y = x.clone(); }\n}\n' > zenith-float-num/src/x.rs
echo 1 > tools/clone_count_baseline.txt; scan --count-only >/dev/null || fail "cfg(test) clone not counted (baseline 1 should hold)"
echo 0 > tools/clone_count_baseline.txt; scan --count-only >/dev/null && fail "clone inside #[cfg(test)] not counted"
printf 'fn t(x: &S) { let y = x.clone(); }\n' > zenith-float-num/src/ops/tests.rs
echo 1 > tools/clone_count_baseline.txt; scan --count-only >/dev/null && fail "clone in ops/tests.rs not counted"
echo 2 > tools/clone_count_baseline.txt; scan --count-only >/dev/null || fail "tests.rs + cfg(test) clones should total 2"
printf 'fn t(x: &S) { let y = x.clone(); }\n' > zenith-float-num/tests/it.rs
echo 2 > tools/clone_count_baseline.txt; scan --count-only >/dev/null && fail "clone in integration tests dir not counted"
echo 3 > tools/clone_count_baseline.txt; scan --count-only >/dev/null || fail "three test clones should total 3"
rm -f zenith-float-num/src/ops/tests.rs zenith-float-num/tests/it.rs; printf 'fn a() {}\n' > zenith-float-num/src/x.rs; echo 0 > tools/clone_count_baseline.txt

# ---- ref ratchets stay scoped to library code: test uses of add_ref etc. do not count ----
cat > zenith-float-num/src/x.rs <<'EOF'
fn u(a: &N, b: &N) { let s = N::add_ref(a, b); let m = N::mul_ref(a, b); }
#[cfg(test)]
mod tests { fn t(a: &N, b: &N) { let s = N::add_ref(a, b); let m = N::div_ref(a, b); } }
EOF
printf 'fn t(a: &N, b: &N) { let s = N::sub_ref(a, b); }\n' > zenith-float-num/src/ops/tests.rs
echo 1 > tools/ref_sum_count_baseline.txt; echo 1 > tools/mul_ref_count_baseline.txt
scan --count-only | grep -q "add_ref/sub_ref use count 1 == baseline" || fail "ref ratchet must stay library-only: $(scan --count-only)"
scan --count-only | grep -q "add_ref/sub_ref uses incl. test code = 3" || fail "with-tests figure not reported: $(scan --count-only)"
scan --count-only | grep -q "mul_ref/div_ref use count 1 == baseline" || fail "mul/div ratchet must stay library-only"
scan --count-only >/dev/null || fail "test uses of ref functions must not trip the library-scoped ratchets"
rm -f zenith-float-num/src/ops/tests.rs; printf 'fn a() {}\n' > zenith-float-num/src/x.rs; echo 0 > tools/ref_sum_count_baseline.txt; echo 0 > tools/mul_ref_count_baseline.txt

# ---- helper-copy counter: uses count, definitions / comments / strings / use lines do not; ratchet only goes down ----
cat > zenith-float-num/src/x.rs <<'EOF'
impl N { pub(crate) fn exact_copy_ext(&self) -> N { N } pub(crate) fn exact_copy(&self) -> N { N } fn exact_operand(v: &N) -> N { N } }
fn copy_p(x: &N) -> N { N }
fn svd_copy_prec(x: &N) -> N { N }
fn lift_num(x: &N) -> Option<N> { None }
macro_rules! exact_operand { () => {} }
use crate::ext::{exact_copy, copy_p};
fn u(a: &N) { let b = a.exact_copy_ext(); let c = a.exact_copy(); let d = Self::exact_operand(a); let e = copy_p(a); let f = svd_copy_prec(a); let g = lift_num(a); }
// a.exact_copy_ext() in a comment
fn s() { let x = "exact_copy(a)"; }
#[cfg(test)]
mod tests { fn t(a: &N) { let b = a.exact_copy_ext(); } }
EOF
mkdir -p zenith-float-compare/src; printf 'fn r(c: &Rc<X>) { let a = Rc::clone(c); let b = Arc::clone(c); }\n' > zenith-float-compare/src/a.rs
printf 'fn r(c: &Rc<X>) { let a = Rc::clone(c); }\n' > zenith-float-num/src/y.rs   # only the compare crate counts handle shares
echo 0 > tools/clone_count_baseline.txt
echo 9 > tools/helper_copy_count_baseline.txt   # 6 uses + 1 in cfg(test) + Rc::clone + Arc::clone
scan --count-only | grep -q "helper-copy count 9 == baseline" || fail "helper-copy count should be 9: $(scan --count-only)"
NO_CLONE_ROOT="$g" python3 "$here/no_clone.py" --helpers | grep -q "^helper_copies: 9 " || fail "no_clone.py --helpers disagrees: $(NO_CLONE_ROOT="$g" python3 "$here/no_clone.py" --helpers)"
echo 8 > tools/helper_copy_count_baseline.txt; scan --count-only | grep -q "helper-copy count 9 is ABOVE" || fail "helper count above baseline not reported"
scan --count-only >/dev/null && fail "helper count above baseline passed"
echo 12 > tools/helper_copy_count_baseline.txt; scan --count-only | grep -qi "helper-copy count 9 is below baseline 12: lower the baseline" || fail "no hint to lower helper baseline"
scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/helper_copy_count_baseline.txt)" = 9 ] || fail "helper baseline not lowered"
echo 0 > tools/helper_copy_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising helper baseline allowed"; [ "$(tr -d '[:space:]' < tools/helper_copy_count_baseline.txt)" = 0 ] || fail "helper baseline raised"
rm -f tools/helper_copy_count_baseline.txt; scan --count-only >/dev/null && fail "missing helper baseline passed"
rm -f zenith-float-compare/src/a.rs zenith-float-num/src/y.rs; echo 0 > tools/helper_copy_count_baseline.txt

# ---- named exception for the constant caches (tools/clone_allowlist.txt) ----
nc() { NO_CLONE_ROOT="$g" python3 "$here/no_clone.py" "$@"; }
mkdir -p zenith-float-num/src/ops/consts
printf '## Exception: constant caches\n' > tools/NO_CLONE_RULES.md
CC=zenith-float-num/src/ops/consts/pi.rs
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); r } }\nfn other(x: &N) -> N { x.clone() }\n' > $CC
printf 'fn a() {}\n' > zenith-float-num/src/x.rs
entry() { printf '%s\t%s\treason\tJeff\ttools/NO_CLONE_RULES.md\t%s\n' "$1" "$2" "$3"; }
: > tools/clone_allowlist.txt
# 1. an unlisted clone fails (no_clone full scan; the diff scan flags it too, see above)
nc >/dev/null && fail "unlisted clones passed no_clone"
# 2. a listed site passes: list for_prec; the unlisted `other` clone must still fail
entry $CC for_prec 1 > tools/clone_allowlist.txt
nc --list | grep -q "pi.rs:1 \[for_prec\] clone() (excepted)" || fail "listed site not reported as excepted: $(nc --list)"
nc --list | grep -q "pi.rs:2 \[other\] clone()$" || fail "unlisted clone beside a listed one must still be a violation: $(nc --list)"
nc >/dev/null && fail "unlisted clone in the same file passed"
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); r } }\nfn other(x: &N) -> N { x }\n' > $CC
nc >/dev/null || fail "listed site did not pass: $(nc)"
# the exception only covers the listed number of sites in that fn: a second clone in it fails
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); let q = r.clone(); q } }\n' > $CC
nc >/dev/null && fail "second clone in a listed fn passed (sites=1)"
# only .clone() can be excepted: to_vec in the listed fn still fails
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); let v = s.to_vec(); r } }\n' > $CC
nc >/dev/null && fail "to_vec in a listed fn passed"
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); r } }\n' > $CC
# a listed fn name in another file is not excepted
cp $CC zenith-float-num/src/ops/consts/e.rs; nc >/dev/null && fail "same fn name in an unlisted file passed"; rm zenith-float-num/src/ops/consts/e.rs
echo 0 > tools/clone_count_baseline.txt; echo 1 > tools/clone_count_baseline.txt
echo 1 > tools/exception_count_baseline.txt; scan --count-only >/dev/null || fail "one listed site vs exception baseline 1 should pass: $(scan --count-only)"
# 3. a stale entry fails: the listed fn no longer clones
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { self.val } }\n' > $CC
nc | grep -q "STALE allowlist entry" || fail "stale entry not reported by no_clone"
nc >/dev/null && fail "stale entry passed no_clone"
echo 0 > tools/clone_count_baseline.txt
scan --count-only | grep -q "STALE clone allowlist entry" || fail "stale entry not reported by rules_scan: $(scan --count-only)"
scan --count-only >/dev/null && fail "stale entry passed rules_scan --count-only"
# a malformed entry (no doc / sites) or a doc without the exception section fails
printf '%s\t%s\treason\tJeff\n' $CC for_prec > tools/clone_allowlist.txt; nc >/dev/null && fail "4-field entry passed"
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); r } }\n' > $CC
entry $CC for_prec 1 > tools/clone_allowlist.txt; printf 'no section\n' > tools/NO_CLONE_RULES.md; nc >/dev/null && fail "entry with an undocumented exception passed"
printf '## Exception: constant caches\n' > tools/NO_CLONE_RULES.md
nc >/dev/null || fail "restored entry should pass: $(nc)"
# 4. the exception count can only go down
echo 1 > tools/clone_count_baseline.txt
echo 0 > tools/exception_count_baseline.txt; scan --count-only | grep -q "clone exceptions site count 1 is ABOVE the baseline 0" || fail "exception count above baseline not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "exception count above baseline passed"
echo 5 > tools/exception_count_baseline.txt; scan --count-only | grep -qi "clone exceptions site count 1 is below baseline 5: lower the baseline" || fail "no hint to lower exception baseline"
scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/exception_count_baseline.txt)" = 1 ] || fail "exception baseline not lowered"
echo 0 > tools/exception_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising exception baseline allowed"; [ "$(tr -d '[:space:]' < tools/exception_count_baseline.txt)" = 0 ] || fail "exception baseline raised"
# a new entry cannot hide behind the old baseline: two listed sites vs baseline 1 fails
printf 'impl P { pub(crate) fn for_prec(&mut self) -> N { let r = self.val.clone(); r } }\nfn other(x: &N) -> N { x.clone() }\n' > $CC
{ entry $CC for_prec 1; entry $CC other 1; } > tools/clone_allowlist.txt
echo 2 > tools/clone_count_baseline.txt
echo 1 > tools/exception_count_baseline.txt; scan --count-only | grep -q "clone exceptions site count 2 is ABOVE the baseline 1" || fail "adding an exception entry beyond the baseline not caught: $(scan --count-only)"
scan --count-only >/dev/null && fail "extra exception entry passed"
rm -f tools/exception_count_baseline.txt; scan --count-only >/dev/null && fail "missing exception baseline passed"

echo "rules_scan selftest: ok"
