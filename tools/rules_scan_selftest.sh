#!/usr/bin/env bash
# Self-test for tools/rules_scan.sh: fixture diffs prove every banned pattern is flagged on ADDED lines,
# that removed/context lines, comments, strings, test code and non-lib paths are not, that the
# allowlist works, that staged/working/range modes work, and that the clone-count baseline only goes down.
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
export RULES_SCAN_ROOT="$t"; mkdir -p "$t/tools"; : > "$t/tools/rules_allowlist.txt"; echo 0 > "$t/tools/clone_count_baseline.txt"; echo 0 > "$t/tools/ref_sum_count_baseline.txt"; echo 0 > "$t/tools/mul_ref_count_baseline.txt"
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

# not flagged: comment, string, removed line, context line, deny attr, unwrap_or, unsafe_code, test paths, cfg(test) not visible in diff-file mode so paths only
mkdiff $LIB '// x.clone() unsafe unwrap() f64' '/// doc: panic! #[allow]' 'let s = "x.clone() unwrap() f32";' '#![forbid(unsafe_code)]' 'let v = r.unwrap_or(0);' 'let v = r.expect_err_x;' > "$t/ok.diff"
scan --diff-file "$t/ok.diff" >/dev/null || fail "comment/string/benign lines flagged: $(scan --diff-file "$t/ok.diff")"
for p in zenith-float-num/tests/x.rs benches/b.rs examples/e.rs tools/t.rs zenith-float-num/src/ops/tests.rs README.md; do
  mkdiff $p 'let y = x.clone(); unsafe {} r.unwrap(); let a: f64 = 1.0;' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null || fail "non-lib path flagged: $p"
done
printf 'diff --git a/%s b/%s\n--- a/%s\n+++ b/%s\n@@ -1,2 +1,1 @@\n-let y = x.clone();\n context\n' $LIB $LIB $LIB $LIB > "$t/rm.diff"
scan --diff-file "$t/rm.diff" >/dev/null || fail "removed line flagged"

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
cp "$t/tools/rules_allowlist.txt" tools/; echo 0 > tools/clone_count_baseline.txt; echo 0 > tools/ref_sum_count_baseline.txt; echo 0 > tools/mul_ref_count_baseline.txt
printf 'fn a() {}\n' > zenith-float-num/src/x.rs; git add -A; git commit -qm base; base=$(git rev-parse HEAD)
export RULES_SCAN_ROOT="$g"; unset RULES_SCAN_SKIP_COUNT
scan "$base..HEAD" >/dev/null || fail "empty range flagged"
printf 'fn a() {}\nfn b(x: &S) { let y = x.clone(); }\n' > zenith-float-num/src/x.rs
scan --working >/dev/null && fail "working diff clone not flagged"
git add -A; scan --staged >/dev/null && fail "staged diff clone not flagged"
scan --working >/dev/null && fail "working mode (git diff HEAD) must still see staged additions"
git commit -qm bad; scan "$base..HEAD" >/dev/null && fail "range clone not flagged"
# cfg(test) module in the tip file is exempt
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(x: &S) { let y = x.clone(); let _ = r.unwrap(); }\n}\n' > zenith-float-num/src/x.rs
git add -A; git commit -qm t; scan "$base..HEAD" >/dev/null || fail "cfg(test) additions flagged: $(scan "$base..HEAD")"

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
echo "rules_scan selftest: ok"
