#!/usr/bin/env bash
# Self-test for tools/rules_scan.sh: fixture diffs prove every banned pattern is flagged on ADDED lines,
# that removed/context lines, comments, strings, test code and non-lib paths are not, that the
# allowlist works, that staged/working/range modes work, that the clone-count, helper-copy baselines only go down, that the borrowing ops' call graph holds no copy,
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
export RULES_SCAN_ROOT="$t"; mkdir -p "$t/tools"; : > "$t/tools/rules_allowlist.txt"; echo 0 > "$t/tools/clone_count_baseline.txt"; echo 0 > "$t/tools/helper_copy_count_baseline.txt"; echo 0 > "$t/tools/exception_count_baseline.txt"; echo 0 > "$t/tools/test_panic_count_baseline.txt"; echo 0 > "$t/tools/neg_via_mul_count_baseline.txt"
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
unwrap|let v = r.unwrap_err();
expect|let v = r.expect("boom");
expect|let v = r.expect_err("boom");
panic|panic!("no");
unreachable|unreachable!();
todo|todo!();
unimplemented|unimplemented!();
f32|let a: f32 = 1.0;
f64|let a = 1.0_f64; let b: f64 = 2.0;
CASES

# not flagged: comment, string, removed line, context line, deny attr, unwrap_or, unsafe_code, non-library paths (tools/, README); cfg(test) not visible in diff-file mode so paths only
mkdiff $LIB '// x.clone() unsafe unwrap() f64' '/// doc: panic! #[allow]' 'let s = "x.clone() unwrap() f32";' '#![forbid(unsafe_code)]' 'let v = r.unwrap_or(0);' 'let v = r.expect_err_x;' > "$t/ok.diff"
scan --diff-file "$t/ok.diff" >/dev/null || fail "comment/string/benign lines flagged: $(scan --diff-file "$t/ok.diff")"
for p in tools/t.rs README.md; do
  mkdiff $p 'let y = x.clone(); unsafe {} r.unwrap(); let a: f64 = 1.0;' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null || fail "non-lib path flagged: $p"
done
# the test dirs (no_clone.TEST_DIRS) are scanned for the clone rules and the PANIC rules (unwrap, expect, panic!, ...), not for unsafe/allow/f32/f64
for p in zenith-float-num/tests/x.rs benches/b.rs examples/e.rs tests/mod.rs nostd-tests/src/lib.rs; do
  for line in 'let y = x.clone();' 'let v = r.unwrap();' 'let v = r.expect("x");' 'panic!("p");' 'unreachable!();' 'todo!();' 'unimplemented!();'; do
    mkdiff $p "$line" > "$t/p.diff"
    scan --diff-file "$t/p.diff" >/dev/null && fail "test-dir file $p: not flagged: $line"
  done
  mkdiff $p 'unsafe {} #[allow(dead_code)] let a: f64 = 1.0; let b: f32 = 1.0;' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null || fail "test-dir file $p: unsafe/allow/f64 flagged: $(scan --diff-file "$t/p.diff")"
done
printf 'diff --git a/%s b/%s\n--- a/%s\n+++ b/%s\n@@ -1,2 +1,1 @@\n-let y = x.clone();\n context\n' $LIB $LIB $LIB $LIB > "$t/rm.diff"
scan --diff-file "$t/rm.diff" >/dev/null || fail "removed line flagged"

# test code under the library dirs is counted now: tests.rs and tests/ benches/ examples/ components are NOT exempt
for p in zenith-float-num/src/ops/tests.rs zenith-float-num/src/tests/x.rs zenith-float-num/src/benches/b.rs; do
  mkdiff $p 'let y = x.clone();' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null && fail "test-code path no longer exempt, but not flagged: $p"
done

# the PANIC rules (unwrap, expect, panic!, unreachable!, todo!, unimplemented!) apply to test files too (Jeff): a panic in test code fails the scan
for p in zenith-float-num/src/ops/tests.rs zenith-float-num/src/tests/x.rs; do
  for line in 'let v = r.unwrap();' 'let w = r.expect("x");' 'panic!("p");' 'unreachable!();' 'todo!();' 'unimplemented!();'; do
    mkdiff $p "$line" > "$t/p.diff"
    scan --diff-file "$t/p.diff" >/dev/null && fail "panic rule not applied to test file $p: $line"
  done
  # ...while the other non-clone rules stay exempt for test files
  mkdiff $p 'let a: f64 = 1.0; unsafe {} #[allow(dead_code)]' > "$t/p.diff"
  scan --diff-file "$t/p.diff" >/dev/null || fail "unsafe/allow/f64 flagged in test file $p: $(scan --diff-file "$t/p.diff")"
done
# assert! macros are fine in test code
mkdiff zenith-float-num/src/ops/tests.rs 'assert!(a == b); assert_eq!(a, b); assert_ne!(a, b); debug_assert!(a);' > "$t/p.diff"
scan --diff-file "$t/p.diff" >/dev/null || fail "assert macros flagged"

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
cp "$t/tools/rules_allowlist.txt" tools/; echo 0 > tools/clone_count_baseline.txt; echo 0 > tools/helper_copy_count_baseline.txt; echo 0 > tools/exception_count_baseline.txt; echo 0 > tools/test_panic_count_baseline.txt; echo 0 > tools/neg_via_mul_count_baseline.txt
printf 'fn a() {}\n' > zenith-float-num/src/x.rs; git add -A; git commit -qm base; base=$(git rev-parse HEAD)
export RULES_SCAN_ROOT="$g"; unset RULES_SCAN_SKIP_COUNT
scan "$base..HEAD" >/dev/null || fail "empty range flagged"
printf 'fn a() {}\nfn b(x: &S) { let y = x.clone(); }\n' > zenith-float-num/src/x.rs
scan --working >/dev/null && fail "working diff clone not flagged"
git add -A; scan --staged >/dev/null && fail "staged diff clone not flagged"
scan --working >/dev/null && fail "working mode (git diff HEAD) must still see staged additions"
git commit -qm bad; scan "$base..HEAD" >/dev/null && fail "range clone not flagged"
# cfg(test) module in the tip file: a RANGE scan leaves the panic-rule enforcement to the test-panic ratchet (checked below), but its clone IS flagged (test code is counted)
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(x: &S) { let y = x.clone(); let _ = r.unwrap(); }\n}\n' > zenith-float-num/src/x.rs
git add -A; git commit -qm t; out="$(RULES_SCAN_SKIP_COUNT=1 scan "$base..HEAD")" && fail "clone in a cfg(test) module not flagged by the diff scan"
echo "$out" | grep -q "\[clone\]" || fail "cfg(test) clone flagged with the wrong label: $out"
echo "$out" | grep -q "\[unwrap\]" && fail "range scan flags a legacy test unwrap per line (the test-panic ratchet covers range mode): $out"
scan "$base..HEAD" >/dev/null && fail "clone in a cfg(test) module not counted by the range scan's clone count (baseline 0)"

# clone-count baseline: may only go down
printf 'fn a(x: &S) { let y = x.clone(); }\n' > zenith-float-num/src/x.rs; git add -A; git commit -qm c
echo 1 > tools/clone_count_baseline.txt
scan --count-only >/dev/null || fail "count == baseline should pass"
echo 0 > tools/clone_count_baseline.txt; scan --count-only >/dev/null && fail "count above baseline passed"
echo 5 > tools/clone_count_baseline.txt; scan --count-only | grep -qi "lower" || fail "no hint to lower baseline"
scan --update-baseline >/dev/null || fail "lowering baseline failed"; [ "$(tr -d '[:space:]' < tools/clone_count_baseline.txt)" = 1 ] || fail "baseline not lowered"
echo 0 > tools/clone_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising baseline allowed"; [ "$(tr -d '[:space:]' < tools/clone_count_baseline.txt)" = 0 ] || fail "baseline raised"

# ---- borrowing ops (add_ref/sub_ref/mul_ref/div_ref): no copy token anywhere in their call graph ----
mkdir -p zenith-float-num/src/mantissa
cat > zenith-float-num/src/ext.rs <<'EOF'
impl N {
    pub fn add_ref(a: &N, b: &N) -> N { a.go(b) }
    pub fn sub_ref(a: &N, b: &N) -> N { a.go(b) }
    pub fn mul_ref(a: &N, b: &N) -> N { a.go(b) }
    pub fn div_ref(a: &N, b: &N) -> N { a.go(b) }
}
fn unrelated(a: &N) -> N { a.clone() }
EOF
kern() { printf 'impl M {\n    fn go(&self, b: &M) -> M {\n        kernel(&self.w)\n    }\n}\nfn kernel(w: &[u64]) -> M {\n    let v = %s;\n    M(v)\n}\n' "$1" > zenith-float-num/src/mantissa/base.rs; }
echo 9 > tools/clone_count_baseline.txt
kern 'w.iter().copied().sum()'
out="$(scan --count-only)" && fail "a .copied() two calls below add_ref passed: $out"
echo "$out" | grep -q 'BORROW zenith-float-num/src/mantissa/base.rs:7 in fn kernel: `.copied(`' || fail "copy site in the borrow call graph not reported: $out"
echo "$out" | grep -q "fn unrelated" && fail "a clone outside the borrow call graph was reported: $out"
for tok in 'w.to_vec()' 'w.iter().cloned().sum()' 'w.clone()'; do
  kern "$tok"
  scan --count-only >/dev/null && fail "$tok in the borrow call graph passed"
done
kern 'w.iter().sum()'; printf '// w.to_vec() in a comment\n' >> zenith-float-num/src/mantissa/base.rs
out="$(scan --count-only)"; echo "$out" | grep -q "BORROW" && fail "a copy-free borrow call graph reported a copy: $out"
echo "$out" | grep -q "0 copy sites" || fail "copy-free borrow graph not reported: $out"
sed -i 's/pub fn div_ref/pub fn div_other/' zenith-float-num/src/ext.rs
scan --count-only | grep -q "div_ref" || fail "a missing borrowing op not reported"
# ---- BORROW per-op fixtures (Harper): each op red direct, red via helper, green; edge cases; exact counts ----
# Layout: ext.rs holds the four ops (lines 2..5), each calling <op>_kern in mantissa/kern.rs, which calls <op>_step.
# Reached set when copy-free: 4 ops + 4 kern + 4 step = 12 fns.
bx_ext() { # bx_ext <op-with-direct-copy|none>
  local o
  { echo 'impl N {'
    for o in add sub mul div; do
      if [ "$o" = "$1" ]; then echo "    pub fn ${o}_ref(a: &N, b: &N) -> N { let c = a.clone(); ${o}_kern(&c, b) }"
      else echo "    pub fn ${o}_ref(a: &N, b: &N) -> N { ${o}_kern(a, b) }"; fi
    done
    echo '}'
    echo 'fn unrelated(a: &N) -> N { a.clone() }'
  } > zenith-float-num/src/ext.rs; }
bx_kern() { # bx_kern <op-with-helper-copy|none> [extra line appended inside every step body]
  local o
  { for o in add sub mul div; do
      echo "fn ${o}_kern(a: &N, b: &N) -> N { ${o}_step(&a.w, &b.w) }"
      if [ "$o" = "$1" ]; then echo "fn ${o}_step(x: &[u64], y: &[u64]) -> N { N(x.iter().copied().zip(y).count()) }"
      else echo "fn ${o}_step(x: &[u64], y: &[u64]) -> N { ${2:-} N(x.iter().zip(y).count()) }"; fi
    done
  } > zenith-float-num/src/mantissa/kern.rs; }
bx_expect() { # bx_expect <want-rc 0|1> <want-hits> <label> [grep-pattern for the BORROW line]
  local out rc; out="$(scan --count-only)"; rc=$?
  [ $rc -eq 0 ] && rc=0 || rc=1
  [ "$rc" = "$1" ] || fail "BORROW fixture $3: rc=$rc want $1: $out"
  local n; n=$(echo "$out" | grep -c 'rules_scan: BORROW ')
  [ "$n" = "$2" ] || fail "BORROW fixture $3: $n BORROW lines, want $2: $out"
  echo "$out" | grep -q "borrowing ops reach 12 fns, $2 copy sites" || fail "BORROW fixture $3: wrong reach/count summary: $out"
  [ -z "${4:-}" ] || echo "$out" | grep -qF "$4" || fail "BORROW fixture $3: message lacks [$4]: $out"
  echo "$out" | grep -q "fn unrelated" && fail "BORROW fixture $3: unrelated fn reported: $out"
  bx_n=$((bx_n + 1)); return 0
}
mkdir -p zenith-float-num/src/mantissa; rm -f zenith-float-num/src/mantissa/base.rs; echo 9 > tools/clone_count_baseline.txt; bx_n=0
ln=1; for o in add sub mul div; do ln=$((ln + 1))
  bx_ext "$o"; bx_kern none
  bx_expect 1 1 "red-direct-$o" "BORROW zenith-float-num/src/ext.rs:$ln in fn ${o}_ref: \`.clone(\` (${o}_ref)"
done
kl=0; for o in add sub mul div; do kl=$((kl + 2))
  bx_ext none; bx_kern "$o"
  bx_expect 1 1 "red-helper-$o" "BORROW zenith-float-num/src/mantissa/kern.rs:$kl in fn ${o}_step: \`.copied(\` (${o}_ref -> ${o}_kern -> ${o}_step)"
  r="$(scan --borrow-reach)"; for want in "src/ext.rs:$((kl / 2 + 1)) ${o}_ref" "mantissa/kern.rs:$((kl - 1)) ${o}_kern" "mantissa/kern.rs:$kl ${o}_step"; do
    echo "$r" | grep -qF "$want" || fail "red-helper-$o: call path link [$want] not in --borrow-reach: $r"; done
done
# green: all four copy-free (the unrelated clone in ext.rs stays and must not trip)
bx_ext none; bx_kern none; bx_expect 0 0 "green-all"
# edge: copy tokens in a comment and a string literal inside reached bodies pass
bx_ext none; bx_kern none 'let _s = "x.clone() y.to_vec()"; // z.copied() w.to_owned()'; bx_expect 0 0 "green-comment-string"
# edge: a #[cfg(test)] fn named like a reached helper, with a clone, is not part of the graph
bx_ext none; bx_kern none; printf '#[cfg(test)]\nmod tests {\n    fn add_step(x: &[u64]) -> Vec<u64> { x.to_vec() }\n}\n' >> zenith-float-num/src/mantissa/kern.rs; bx_expect 0 0 "green-cfg-test-shadow"
# edge: every banned token (and a helper-copy call) in a reached helper fails, one site each
for tt in 'x.to_vec()|.to_vec(' 'x.iter().cloned().count()|.cloned(' 'x.to_owned()|.to_owned(' 'x.clone()|.clone(' 'exact_copy(x)|exact_copy'; do
  tok=${tt%%|*}; want=${tt#*|}
  bx_ext none; bx_kern none "let _c = $tok;"; out="$(scan --count-only)" && fail "BORROW fixture token $tok in 4 reached steps passed: $out"
  echo "$out" | grep -q "borrowing ops reach 12 fns, 4 copy sites" || fail "BORROW fixture token $tok: want 4 sites: $out"
  for o in add sub mul div; do echo "$out" | grep -qF "in fn ${o}_step: \`$want\` (${o}_ref -> ${o}_kern -> ${o}_step)" || fail "BORROW fixture token $tok: no ${o}_ref chain message: $out"; done
  bx_n=$((bx_n + 1))
done
# ---- BORROW hole fixtures (Harper; from the former tools/rules_scan_borrow_holes.sh): add_ref (ext.rs:2) reaches a
# copy the old name-only rule missed; each red asserts the exact site count and the op + reach chain in the message ----
mkdir -p zenith-float-num/src/ops; echo 99 > tools/helper_copy_count_baseline.txt
hx() { # hx <name> <want-sites> <expected BORROW line or -> <add_ref body> <kern.rs content> [path:content]
  { echo 'impl N {'; echo "    pub fn add_ref(a: &N, b: &N) -> N { $4 }"
    for o in sub mul div; do echo "    pub fn ${o}_ref(a: &N, b: &N) -> N { N(0) }"; done; echo '}'; } > zenith-float-num/src/ext.rs
  printf '%b\n' "$5" > zenith-float-num/src/mantissa/kern.rs
  rm -f zenith-float-num/src/ops/h.rs; [ -z "${6:-}" ] || printf '%b\n' "${6#*:}" > "zenith-float-num/src/${6%%:*}"
  local out rc=0 n; out="$(scan --count-only 2>&1)" || rc=1
  if [ "$2" = 0 ]; then [ $rc = 0 ] || fail "BORROW hole fixture $1: copy-free graph failed: $out"
  else [ $rc = 1 ] || fail "BORROW hole fixture $1: copy reachable from add_ref passed: $out"; fi
  n=$(echo "$out" | grep -c 'rules_scan: BORROW ')
  [ "$n" = "$2" ] || fail "BORROW hole fixture $1: $n BORROW lines, want $2: $out"
  echo "$out" | grep -q ", $2 copy sites" || fail "BORROW hole fixture $1: want $2 copy sites: $out"
  [ "$3" = - ] || echo "$out" | grep -qF "rules_scan: BORROW $3" || fail "BORROW hole fixture $1: message lacks [$3]: $out"
  bx_n=$((bx_n + 1)); return 0
}
E=zenith-float-num/src/ext.rs; K=zenith-float-num/src/mantissa/kern.rs; KK='fn k(a: &N, b: &N) -> N { N(0) }'
hx ufcs-clone 1 "$E:2 in fn add_ref: \`::clone\` (add_ref)" 'let c = Clone::clone(a); k(&c, b)' "$KK"
hx space-before-paren 1 "$E:2 in fn add_ref: \`.clone(\` (add_ref)" 'let c = a.clone (); k(&c, b)' "$KK"
hx slice-to_vec-path 1 "$K:1 in fn k: \`::to_vec\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(<[u64]>::to_vec(&a.w).len()) }'
hx vec-from-slice 1 "$K:1 in fn k: \`Vec::from(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(Vec::from(&a.w[..]).len()) }'
hx deref-map-collect 1 "$K:1 in fn k: \`|x|*\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|x| *x).collect(); N(v.len()) }'
hx extend-from-slice 1 "$K:1 in fn k: \`.extend_from_slice(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); v.extend_from_slice(&a.w); N(v.len()) }'
hx fn-as-value 1 "$K:2 in fn cp: \`.clone(\` (add_ref -> k -> cp)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(a.w.iter().map(cp).count()) }\nfn cp(x: &u64) -> u64 { x.clone() }'
hx operator-dispatch 1 "$K:1 in fn add: \`.clone(\` (add_ref -> add)" 'a + b' 'impl core::ops::Add for &N { type Output = N; fn add(self, b: &N) -> N { self.clone() } }'
hx macro-body 1 "$K:1 in fn dup!: \`.clone(\` (add_ref -> k -> dup!)" 'k(a, b)' 'macro_rules! dup { ($x:expr) => { $x.clone() } }\nfn k(a: &N, b: &N) -> N { dup!(a) }'
hx helper-outside-scope 1 "zenith-float-num/src/ops/h.rs:1 in fn deep: \`.clone(\` (add_ref -> k -> deep)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { crate::ops::h::deep(a) }' 'ops/h.rs:pub fn deep(a: &N) -> N { a.clone() }'
hx multiline-call 1 "$K:3 in fn k: \`.clone(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N {\n    a\n        .clone()\n}'
hx green-multiline-block-comment 0 - 'k(a, b)' 'fn k(a: &N, b: &N) -> N {\n    /* note:\n       a.clone() is banned here */\n    N(0)\n}'
# extra variants (Harper, against 581cdfe): caught, kept as regression reds
hx v-ufcs-clone-ref 1 "$E:2 in fn add_ref: \`::clone\` (add_ref)" 'let c = Clone::clone(&a); k(&c, b)' "$KK"
hx v-qualified-as-clone 1 "$E:2 in fn add_ref: \`::clone\` (add_ref)" 'let c = <N as Clone>::clone(a); k(&c, b)' "$KK"
hx v-slice-to_owned 1 "$K:1 in fn k: \`.to_owned(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v = a.w[..].to_owned(); N(v.len()) }'
hx v-copied-collect 1 "$K:1 in fn k: \`.copied(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().copied().collect(); N(v.len()) }'
hx v-closure-in-let 1 "$K:1 in fn k: \`.clone(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let f = |x: &N| x.clone(); f(a) }'
hx v-trait-default-method 1 "$K:1 in fn dup: \`.clone(\` (add_ref -> dup)" 'a.dup()' 'trait Dup { fn dup(&self) -> N { self.clone() } }\nimpl Dup for N {}'
hx v-nested-macros 1 "$K:1 in fn inner!: \`.clone(\` (add_ref -> k -> outer! -> inner!)" 'k(a, b)' 'macro_rules! inner { ($x:expr) => { $x.clone() } }\nmacro_rules! outer { ($x:expr) => { inner!($x) } }\nfn k(a: &N, b: &N) -> N { outer!(a) }'
hx v-map-clone-path 1 "$K:1 in fn k: \`::clone\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(Clone::clone).collect(); N(v.len()) }'
# round 3 (Harper): the 8 holes of 5cf3df2's probe (tools/rules_scan_borrow_holes.sh), red against 581cdfe
hx h3-ref-pattern-closure 1 "$K:1 in fn k: \`|&\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|&x| x).collect(); N(v.len()) }'
hx h3-copy-from-slice 1 "$K:1 in fn k: \`.copy_from_slice(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = vec![0u64; a.w.len()]; v.copy_from_slice(&a.w); N(v.len()) }'
hx h3-slice-into-vec 1 "$K:1 in fn k: \`as_slice().into(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.as_slice().into(); N(v.len()) }'
hx h3-clone-into 1 "$K:1 in fn k: \`.clone_into(\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); a.w.as_slice().clone_into(&mut v); N(v.len()) }'
hx h3-vec-macro-copy 1 "$K:1 in fn k: \`vec![<operand>;n]\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v = vec![a.w[0]; 4]; N(v.len()) }'
hx h3-typed-braced-deref 1 "$K:1 in fn k: \`|x:&u64|->u64{*\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|x: &u64| -> u64 { *x }).collect(); N(v.len()) }'
hx h3-extend-ref-iter 1 "$K:1 in fn k: \`.extend(<borrowed>)\` (add_ref -> k)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); v.extend(a.w.iter()); N(v.len()) }'
hx h3-generic-turbofish 1 "$K:2 in fn g: \`.clone(\` (add_ref -> k -> g)" 'k(a, b)' 'fn k(a: &N, b: &N) -> N { g::<N>(a) }\nfn g<T: Clone>(x: &T) -> T { x.clone() }'
# named exceptions (b9812ca BORROW_EXCEPTIONS): mantissa/base.rs abs_add allows 2 `.copy_from_slice(`
B=zenith-float-num/src/mantissa/base.rs; AA='fn abs_add(a: &N, b: &N) -> N { let mut r = vec![0u64; 4]; '
printf '%s\n' "${AA}r.copy_from_slice(&a.w); r.copy_from_slice(&b.w); N(r.len()) }" > $B
hx x-exception-at-limit 0 - 'abs_add(a, b)' "$KK"
# (a) a new copy in the same fn as an exception still fails: one more of the same token (all 3 sites reported), or another token
printf '%s\n' "${AA}r.copy_from_slice(&a.w); r.copy_from_slice(&b.w); r.copy_from_slice(&a.w); N(r.len()) }" > $B
hx x-exception-same-token-over 3 "$B:1 in fn abs_add: \`.copy_from_slice(\` (add_ref -> abs_add)" 'abs_add(a, b)' "$KK"
printf '%s\n' "${AA}r.copy_from_slice(&a.w); r.copy_from_slice(&b.w); let c = a.clone(); N(c.w.len()) }" > $B
hx x-exception-other-token 1 "$B:1 in fn abs_add: \`.clone(\` (add_ref -> abs_add)" 'abs_add(a, b)' "$KK"
# ...and the same fn name in another file is not covered by the exception
rm -f $B; mkdir -p zenith-float-num/src/ops
hx x-exception-other-file 1 "zenith-float-num/src/ops/h.rs:1 in fn abs_add: \`.copy_from_slice(\` (add_ref -> abs_add)" 'abs_add(a, b)' "$KK" "ops/h.rs:${AA}r.copy_from_slice(&a.w); N(r.len()) }"
# (b) an exception whose sites shrink is reported as stale (a note; the rule does not fail on it)
printf '%s\n' "${AA}r.copy_from_slice(&a.w); N(r.len()) }" > $B
hx x-exception-stale 0 - 'abs_add(a, b)' "$KK"
scan --count-only | grep -qF "borrow exception $B abs_add \`.copy_from_slice(\` allows 2, found 1: lower it" || fail "stale abs_add exception (2 allowed, 1 found) not reported"
rm -f $B
# (c) the same-name fallback: an unresolved generic receiver follows every crate fn `frob`; copy-free ones in an
# unrelated module are reached but raise no alarm
hx x-fallback-copy-free 0 - 'k(a, b)' 'fn k<T: Fro>(t: &T, b: &N) -> N { t.frob() }\nimpl M { fn frob(&self) -> N { N(1) } }' 'ops/h.rs:impl Z { fn frob(&self) -> N { N(0) } }'
scan --borrow-reach | grep -q "ops/h.rs:1 frob" || fail "same-name fallback did not reach ops/h.rs frob"
rm -rf zenith-float-num/src/ops; echo 0 > tools/helper_copy_count_baseline.txt
[ "$bx_n" = 50 ] || fail "BORROW fixtures: ran $bx_n, want 50"
echo "rules_scan selftest: BORROW per-op fixtures ok ($bx_n)"
rm -f zenith-float-num/src/mantissa/kern.rs
rm -rf zenith-float-num/src/mantissa zenith-float-num/src/ext.rs; echo 0 > tools/clone_count_baseline.txt
# ---- BORROW holes: copies the rule must see (UFCS, spacing, slice paths, deref-collect, fn values, operators,
# macros, helpers anywhere in the crate src) and a copy-free graph with a copy token only inside a block comment ----
mkdir -p zenith-float-num/src/mantissa zenith-float-num/src/ops
echo 99 > tools/clone_count_baseline.txt; echo 99 > tools/helper_copy_count_baseline.txt
bred=""
bcase() { # bcase <want: fail|pass> <name> <add_ref body> <mantissa/kern.rs content> [path:content]
  { echo 'impl N {'; echo "    pub fn add_ref(a: &N, b: &N) -> N { $3 }"
    for o in sub mul div; do echo "    pub fn ${o}_ref(a: &N, b: &N) -> N { N(0) }"; done; echo '}'; } > zenith-float-num/src/ext.rs
  printf '%b\n' "$4" > zenith-float-num/src/mantissa/kern.rs
  rm -f zenith-float-num/src/ops/h.rs; [ -z "${5:-}" ] || printf '%b\n' "${5#*:}" > "zenith-float-num/src/${5%%:*}"
  bout="$(scan --count-only 2>&1)"; local got=pass
  echo "$bout" | grep -q "BORROW" && got=fail
  [ "$got" = "$1" ] || { bred="$bred $2"; echo "BORROW fixture red: $2 (want $1, got $got)"; }
}
bcase fail ufcs-clone 'let c = Clone::clone(a); k(&c, b)' 'fn k(a: &N, b: &N) -> N { N(0) }'
bcase fail space-before-paren 'let c = a.clone (); k(&c, b)' 'fn k(a: &N, b: &N) -> N { N(0) }'
bcase fail slice-to_vec-path 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(<[u64]>::to_vec(&a.w).len()) }'
bcase fail vec-from-slice 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(Vec::from(&a.w[..]).len()) }'
bcase fail deref-map-collect 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|x| *x).collect(); N(v.len()) }'
bcase fail extend-from-slice 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); v.extend_from_slice(&a.w); N(v.len()) }'
bcase fail fn-as-value 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(a.w.iter().map(cp).count()) }\nfn cp(x: &u64) -> u64 { x.clone() }'
bcase fail operator-dispatch 'a + b' 'impl core::ops::Add for &N { type Output = N; fn add(self, b: &N) -> N { self.clone() } }'
bcase fail macro-body 'k(a, b)' 'macro_rules! dup { ($x:expr) => { $x.clone() } }\nfn k(a: &N, b: &N) -> N { dup!(a) }'
bcase fail helper-outside-scope 'k(a, b)' 'fn k(a: &N, b: &N) -> N { crate::ops::h::deep(a) }' 'ops/h.rs:pub fn deep(a: &N) -> N { a.clone() }'
bcase pass multiline-block-comment 'k(a, b)' 'fn k(a: &N, b: &N) -> N {\n    /* note:\n       a.clone() is banned here */\n    let s = "/* not a comment"; N(s.len())\n}\nfn other(a: &N) -> N { a.clone() }'
# the failure names the op and the call chain to the copy site
bcase fail reach-path 'k(a, b)' 'fn k(a: &N, b: &N) -> N { deep(a) }\nfn deep(a: &N) -> N { a.clone() }'
echo "$bout" | grep -q 'BORROW zenith-float-num/src/mantissa/kern.rs:2 in fn deep: `.clone(` (add_ref -> k -> deep)' || { bred="$bred reach-path-message"; echo "BORROW fixture red: reach-path-message: $bout"; }
rm -rf zenith-float-num/src/mantissa zenith-float-num/src/ops/h.rs zenith-float-num/src/ext.rs
echo 0 > tools/clone_count_baseline.txt; echo 0 > tools/helper_copy_count_baseline.txt
[ -z "$bred" ] || fail "BORROW fixtures red:$bred"
# ---- test code is counted by the clone count (cfg(test) items, tests.rs, test dirs) ----
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

# ---- test-panic ratchet: unwrap/expect/panic!/unreachable!/todo!/unimplemented! in test code may only go down ----
: > tools/clone_allowlist.txt; echo 0 > tools/exception_count_baseline.txt; rm -rf zenith-float-num/src/ops/consts
echo 0 > tools/clone_count_baseline.txt
rm -rf zenith-float-num/src/x.rs zenith-float-num/src/ops zenith-float-num/tests; mkdir -p zenith-float-num/tests zenith-float-num/src/ops
cat > zenith-float-num/src/x.rs <<'EOF'
fn lib(r: R) -> u8 { assert!(true); r.unwrap_or(0) }
// r.unwrap() in a comment
fn s() { let x = "r.unwrap() panic!"; }
#[cfg(test)]
mod tests { fn t(r: R) { let a = r.unwrap(); let b = r.expect("x"); panic!("p"); } }
EOF
printf 'fn t(r: R) { let a = r.unwrap_err(); unreachable!(); todo!(); unimplemented!(); }\n' > zenith-float-num/tests/it.rs
printf 'fn t(r: R) { let a = r.unwrap(); }\n' > zenith-float-num/src/ops/tests.rs
# 3 + 4 + 1 = 8 sites (the library fn's unwrap_or, the comment and the string do not count)
echo 8 > tools/test_panic_count_baseline.txt
scan --count-only | grep -q "test-panic site count 8 == baseline" || fail "test-panic count should be 8: $(scan --count-only)"
scan --count-only >/dev/null || fail "test-panic count == baseline should pass: $(scan --count-only)"
scan --list-test-panics | grep -q "zenith-float-num/tests/it.rs:1 \[unwrap\]" || fail "--list-test-panics misses an integration-test site: $(scan --list-test-panics)"
echo 7 > tools/test_panic_count_baseline.txt; scan --count-only | grep -q "test-panic site count 8 is ABOVE the baseline 7" || fail "test-panic count above baseline not reported"
scan --count-only >/dev/null && fail "a new unwrap/panic in test code passed the ratchet"
echo 12 > tools/test_panic_count_baseline.txt; scan --count-only | grep -qi "test-panic site count 8 is below baseline 12: lower the baseline" || fail "no hint to lower the test-panic baseline"
scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/test_panic_count_baseline.txt)" = 8 ] || fail "test-panic baseline not lowered"
echo 0 > tools/test_panic_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising the test-panic baseline allowed"; [ "$(tr -d '[:space:]' < tools/test_panic_count_baseline.txt)" = 0 ] || fail "test-panic baseline raised"
rm -f tools/test_panic_count_baseline.txt; scan --count-only >/dev/null && fail "missing test-panic baseline passed"
# a library (non-test) unwrap is not a test-panic site: it is the diff scan's job (flagged there)
printf 'fn lib(r: R) -> u8 { r.unwrap() }\n' > zenith-float-num/src/x.rs; rm -f zenith-float-num/tests/it.rs zenith-float-num/src/ops/tests.rs
echo 0 > tools/test_panic_count_baseline.txt; scan --count-only >/dev/null || fail "library unwrap counted as a test-panic site"
# staged diff: a new unwrap in a cfg(test) module is flagged per line
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(r: R) { let _ = r.unwrap(); }\n}\n' > zenith-float-num/src/x.rs
git add -A; scan --staged >/dev/null && fail "unwrap in a cfg(test) module not flagged by the staged scan"
printf 'fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t(r: R) -> Result<(), E> { let _v = r?; assert!(true); Ok(()) }\n}\n' > zenith-float-num/src/x.rs
git add -A; scan --staged >/dev/null || fail "Result-returning test flagged: $(scan --staged)"

# neg_via_mul use-count baseline: counts uses in all scanned code, tests included (definitions, comments, strings and
# names that merely start with the word excluded), may only go down
printf 'impl N { pub fn neg_via_mul(&self) -> N { go(self) } }\nfn u(a: &N, b: &N) { let s = a.neg_via_mul(); let d = N::neg_via_mul(b); }\n// a.neg_via_mul() in a comment\n/// [`N::neg_via_mul`]\nfn s() { let x = "a.neg_via_mul()"; }\nfn neg_via_mul_equals_neg() {}\n' > zenith-float-num/src/x.rs
printf '#[test]\nfn t(a: &N) { let s = a.neg_via_mul(); }\n' > zenith-float-num/tests/nv.rs
rm -f tools/test_panic_count_baseline.txt; echo 0 > tools/test_panic_count_baseline.txt
echo 3 > tools/neg_via_mul_count_baseline.txt
scan --count-only | grep -q "neg_via_mul use count 3 == baseline" || fail "neg_via_mul count should be 3 (2 in src, 1 in a test dir; definition, comments, strings, longer names excluded): $(scan --count-only)"
echo 2 > tools/neg_via_mul_count_baseline.txt; scan --count-only | grep -q "neg_via_mul use count 3 is ABOVE the baseline 2" || fail "neg_via_mul count above baseline not reported"
scan --count-only >/dev/null && fail "a new neg_via_mul use passed the ratchet"
echo 9 > tools/neg_via_mul_count_baseline.txt; scan --count-only | grep -qi "neg_via_mul use count 3 is below baseline 9: lower the baseline" || fail "no hint to lower the neg_via_mul baseline"
scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/neg_via_mul_count_baseline.txt)" = 3 ] || fail "neg_via_mul baseline not lowered"
echo 0 > tools/neg_via_mul_count_baseline.txt; scan --update-baseline >/dev/null && fail "raising the neg_via_mul baseline allowed"; [ "$(tr -d '[:space:]' < tools/neg_via_mul_count_baseline.txt)" = 0 ] || fail "neg_via_mul baseline raised"
rm -f tools/neg_via_mul_count_baseline.txt; scan --count-only >/dev/null && fail "missing neg_via_mul baseline passed"
rm -f zenith-float-num/tests/nv.rs; printf 'fn a() {}\n' > zenith-float-num/src/x.rs; echo 0 > tools/neg_via_mul_count_baseline.txt

# ---- the one approved Ziv working-copy helper (exception baseline 0 -> 1): exactly one `fn ziv_pass_operand`, in
# zenith-float-num/src/num/ziv_operand.rs, with a CODEOWNERS entry; a second helper, a second definition, the wrong file or a
# missing CODEOWNERS entry fail; the definition counts 1 in the exception total, each use counts in the helper-copy counter ----
ZH=zenith-float-num/src/num/ziv_operand.rs
mkdir -p zenith-float-num/src/num .github
: > tools/clone_allowlist.txt; echo 0 > tools/helper_copy_count_baseline.txt
printf '/zenith-float-num/src/num/ziv_operand.rs @jscarr64\n' > .github/CODEOWNERS
printf 'impl N { pub(crate) fn ziv_pass_operand(&self, p: usize) -> N { go(self, p) } }\n' > $ZH
echo 1 > tools/exception_count_baseline.txt
scan --count-only | grep -q "clone exceptions site count 1 == baseline" || fail "the Ziv helper definition should count 1 in the exception total: $(scan --count-only)"
scan --count-only >/dev/null || fail "one helper, right file, CODEOWNERS entry, baseline 1 should pass: $(scan --count-only)"
echo 0 > tools/exception_count_baseline.txt
scan --count-only | grep -q "clone exceptions site count 1 is ABOVE the baseline 0" || fail "the helper without the exception baseline not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "the Ziv helper passed with exception baseline 0"
echo 1 > tools/exception_count_baseline.txt
# a second helper (any ziv + copy/operand/pass/widen/... name) is rejected even when the baseline is raised to cover it
for other in 'fn ziv_pass_operand_wide(x: &N) -> N { go(x) }' 'fn ziv_working_copy(x: &N) -> N { go(x) }' 'fn widen_for_ziv(x: &N) -> N { go(x) }' 'fn ziv_clone(x: &N) -> N { go(x) }'; do
  printf '%s\n' "$other" > zenith-float-num/src/x.rs
  echo 2 > tools/exception_count_baseline.txt
  scan --count-only | grep -q "second Ziv working-copy helper" || fail "a second helper was not rejected ($other): $(scan --count-only)"
  scan --count-only >/dev/null && fail "a second Ziv helper passed: $other"
done
printf 'fn a() {}\n' > zenith-float-num/src/x.rs; echo 1 > tools/exception_count_baseline.txt
# the same name defined twice, or in another file, is a second helper
printf 'fn ziv_pass_operand(x: &N) -> N { go(x) }\n' > zenith-float-num/src/x.rs; echo 2 > tools/exception_count_baseline.txt
scan --count-only | grep -q "second Ziv working-copy helper" || fail "a second fn ziv_pass_operand in another file not rejected: $(scan --count-only)"
scan --count-only >/dev/null && fail "ziv_pass_operand defined in a second file passed"
printf 'fn a() {}\n' > zenith-float-num/src/x.rs; printf 'impl N { fn ziv_pass_operand(&self) -> N { a() } fn ziv_pass_operand(&self, p: usize) -> N { b() } }\n' > $ZH
scan --count-only | grep -q "second Ziv working-copy helper" || fail "ziv_pass_operand defined twice in the one file not rejected: $(scan --count-only)"
scan --count-only >/dev/null && fail "ziv_pass_operand defined twice passed"
# the only definition in the wrong file fails (the file is fixed by the scan), and a missing CODEOWNERS entry fails
printf 'impl N { pub(crate) fn ziv_pass_operand(&self, p: usize) -> N { go(self, p) } }\n' > $ZH; echo 1 > tools/exception_count_baseline.txt
mv $ZH zenith-float-num/src/y.rs
scan --count-only | grep -q "ziv_pass_operand must live in zenith-float-num/src/num/ziv_operand.rs" || fail "the helper in the wrong file not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "the helper in the wrong file passed"
mv zenith-float-num/src/y.rs $ZH
: > .github/CODEOWNERS
scan --count-only | grep -q "no CODEOWNERS entry" || fail "missing CODEOWNERS entry for the helper file not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "the helper without a CODEOWNERS entry passed"
printf '/zenith-float-num/src/num/ziv_operand.rs @jscarr64\n' > .github/CODEOWNERS
scan --count-only >/dev/null || fail "restored helper setup should pass: $(scan --count-only)"
# the exception baseline can only go down: with the helper deleted the total is 0 and the baseline can be lowered, never raised
echo 5 > tools/exception_count_baseline.txt; scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/exception_count_baseline.txt)" = 1 ] || fail "exception baseline not lowered to the helper total"
rm -f $ZH; scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/exception_count_baseline.txt)" = 0 ] || fail "exception baseline not lowered to 0 once the helper is gone"
printf 'impl N { pub(crate) fn ziv_pass_operand(&self, p: usize) -> N { go(self, p) } }\n' > $ZH
scan --update-baseline >/dev/null && fail "raising the exception baseline for the helper allowed"; [ "$(tr -d '[:space:]' < tools/exception_count_baseline.txt)" = 0 ] || fail "exception baseline raised for the helper"
# every use of the helper counts in the helper-copy counter; the definition does not
printf 'fn u(m: &N) { let a = m.ziv_pass_operand(p); let b = N::ziv_pass_operand(m, q); }\n// m.ziv_pass_operand(p) in a comment\n' > zenith-float-num/src/x.rs
printf 'zenith-float-num/src/x.rs\tu\tthe kernel consumes its argument each pass\n' > tools/ziv_shape_a_sites.txt
printf '/zenith-float-num/src/num/ziv_operand.rs @jscarr64\n/tools/ziv_shape_a_sites.txt @jscarr64\n' > .github/CODEOWNERS
echo 1 > tools/exception_count_baseline.txt; echo 1 > tools/helper_copy_count_baseline.txt; echo 2 > tools/clone_helper_combined_baseline.txt
scan --count-only | grep -q "helper-copy count 2 is ABOVE the baseline 1" || fail "uses of ziv_pass_operand not counted by the helper-copy counter: $(scan --count-only)"
echo 2 > tools/helper_copy_count_baseline.txt; scan --count-only >/dev/null || fail "two uses vs helper-copy baseline 2 should pass: $(scan --count-only)"

# ---- helper-copy ruling (Jeff, 2026-10-16): combined clone+helper_copy ratchet, per-commit rule, pinned shape-A sites ----
# (a) a use outside the pinned list fails; a missing list with a use fails
printf 'fn u(m: &N) { let a = m.ziv_pass_operand(p); }\nfn w(m: &N) { let b = m.ziv_pass_operand(q); }\n' > zenith-float-num/src/x.rs
scan --count-only | grep -q "ziv_pass_operand use at zenith-float-num/src/x.rs:2 (fn w) is not on the pinned shape-A list" || fail "a use outside the pinned list not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "a ziv_pass_operand use outside the pinned list passed"
printf 'zenith-float-num/src/x.rs\tu\tr\nzenith-float-num/src/x.rs\tw\tr\n' > tools/ziv_shape_a_sites.txt
scan --count-only >/dev/null || fail "both uses pinned should pass: $(scan --count-only)"
mv tools/ziv_shape_a_sites.txt tools/pinned.bak
scan --count-only | grep -q "is not on the pinned shape-A list" || fail "a use with no pinned list not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "a use with no pinned list passed"
mv tools/pinned.bak tools/ziv_shape_a_sites.txt
# (b) the pinned list must be covered by CODEOWNERS for @jscarr64
printf '/zenith-float-num/src/num/ziv_operand.rs @jscarr64\n' > .github/CODEOWNERS
scan --count-only | grep -q "ziv_shape_a_sites.txt is not covered by CODEOWNERS for @jscarr64" || fail "uncovered pinned list not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "pinned list without CODEOWNERS passed"
printf '/zenith-float-num/src/num/ziv_operand.rs @jscarr64\n/tools/ziv_shape_a_sites.txt @someone\n' > .github/CODEOWNERS
scan --count-only >/dev/null && fail "pinned list owned by someone else passed"
printf '/zenith-float-num/src/num/ziv_operand.rs @jscarr64\n/tools/ziv_shape_a_sites.txt @jscarr64\n' > .github/CODEOWNERS
# (c) the combined clone+helper_copy count may only go down, even when each separate baseline is raised to cover it
printf 'fn u(m: &N) { let a = m.ziv_pass_operand(p); }\nfn w(m: &N) { let b = m.ziv_pass_operand(q); let c = m.clone(); }\n' > zenith-float-num/src/x.rs
echo 1 > tools/clone_count_baseline.txt
scan --count-only | grep -q "clone+helper_copy combined count 3 is ABOVE the baseline 2" || fail "combined rise not reported: $(scan --count-only)"
scan --count-only >/dev/null && fail "combined clone+helper_copy rise passed"
echo 3 > tools/clone_helper_combined_baseline.txt; scan --count-only >/dev/null || fail "combined at baseline should pass: $(scan --count-only)"
echo 9 > tools/clone_helper_combined_baseline.txt; scan --update-baseline >/dev/null; [ "$(tr -d '[:space:]' < tools/clone_helper_combined_baseline.txt)" = 3 ] || fail "combined baseline not lowered"
echo 2 > tools/clone_helper_combined_baseline.txt; scan --update-baseline >/dev/null && fail "raising the combined baseline allowed"; [ "$(tr -d '[:space:]' < tools/clone_helper_combined_baseline.txt)" = 2 ] || fail "combined baseline raised"
rm -f tools/clone_helper_combined_baseline.txt; scan --count-only >/dev/null && fail "missing combined baseline passed"
# (d) per commit: helper_copy may rise only in a commit that lowers clones by at least as many
printf 'fn u(m: &N) { let a = m.clone(); }\nfn w(m: &N) { let b = m.clone(); }\n' > zenith-float-num/src/x.rs
echo 2 > tools/clone_count_baseline.txt; echo 0 > tools/helper_copy_count_baseline.txt; echo 2 > tools/clone_helper_combined_baseline.txt
git add -A; git commit -qm r0; r0=$(git rev-parse HEAD)
printf 'fn u(m: &N) { let a = m.ziv_pass_operand(p); }\nfn w(m: &N) { let b = m.clone(); }\n' > zenith-float-num/src/x.rs
echo 1 > tools/clone_count_baseline.txt; echo 1 > tools/helper_copy_count_baseline.txt
git add -A; git commit -qm swap1
scan "$r0..HEAD" >/dev/null || fail "a clone replaced by one helper use in one commit should pass: $(scan "$r0..HEAD")"
r1=$(git rev-parse HEAD)
printf 'fn u(m: &N) { let a = m.ziv_pass_operand(p); }\nfn w(m: &N) { let b = m.clone(); }\nfn x(m: &N) { let c = exact_copy(m); }\n' > zenith-float-num/src/x.rs
echo 2 > tools/helper_copy_count_baseline.txt; echo 3 > tools/clone_helper_combined_baseline.txt
git add -A; git commit -qm rise
printf 'fn u(m: &N) { let a = m.ziv_pass_operand(p); }\nfn w(m: &N) { let b = m.ziv_pass_operand(p); }\nfn x(m: &N) { let c = exact_copy(m); }\n' > zenith-float-num/src/x.rs
printf 'zenith-float-num/src/x.rs\tu\tr\nzenith-float-num/src/x.rs\tw\tr\n' > tools/ziv_shape_a_sites.txt
echo 0 > tools/clone_count_baseline.txt; echo 3 > tools/helper_copy_count_baseline.txt
git add -A; git commit -qm swap2
out="$(scan "$r1..HEAD")" && fail "a helper_copy rise without a clone drop passed (net over the range is fine, the commit is not): $out"
echo "$out" | grep -q "raises helper_copy by 1 but lowers clones by only 0" || fail "per-commit helper rise not reported: $out"
# a documented exception (tools/helper_copy_commit_exceptions.txt: <full sha> <allowed count> <reason>) covers that
# commit only up to its allowed count; a line without a count is not an exception
rise=$(git rev-parse HEAD~1)
printf '%s 1 rise by one (selftest)\n' "$rise" > tools/helper_copy_commit_exceptions.txt
out="$(scan "$r1..HEAD")" || fail "an exception allowing the commit's rise of 1 should pass: $out"
echo "$out" | grep -q "raises helper_copy by 1 (documented exception, allowed 1" || fail "documented exception not reported: $out"
printf '%s 0 rise by one (selftest)\n' "$rise" > tools/helper_copy_commit_exceptions.txt
out="$(scan "$r1..HEAD")" && fail "a rise of 1 above an allowed count of 0 passed: $out"
echo "$out" | grep -q "raises helper_copy by 1, above its documented exception count 0" || fail "rise above the allowed count not reported: $out"
printf '%s rise by one, no count (selftest)\n' "$rise" > tools/helper_copy_commit_exceptions.txt
scan "$r1..HEAD" >/dev/null && fail "an exception line without an allowed count was honoured"
rm -f tools/helper_copy_commit_exceptions.txt
git rev-parse HEAD > tools/helper_copy_rule_base.txt; git add -A; git commit -qm base-after
scan "$(git rev-parse HEAD~1)..HEAD" >/dev/null || fail "commits before the rule base must not be checked: $(scan "$(git rev-parse HEAD~1)..HEAD")"
git rm -q tools/helper_copy_rule_base.txt; git commit -qm unbase
rm -f tools/ziv_shape_a_sites.txt tools/clone_helper_combined_baseline.txt; echo 0 > tools/clone_count_baseline.txt

rm -f $ZH; rm -rf zenith-float-num/src/num; : > .github/CODEOWNERS; printf 'fn a() {}\n' > zenith-float-num/src/x.rs
echo 0 > tools/exception_count_baseline.txt; echo 0 > tools/helper_copy_count_baseline.txt

echo "rules_scan selftest: ok"
