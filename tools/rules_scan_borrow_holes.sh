#!/usr/bin/env bash
# Probe (Harper, for Rusty): copies reachable from add_ref/sub_ref/mul_ref/div_ref that the BORROW rule should
# reject. NOT wired into rules_scan_selftest.sh: each "HOLE" line is a red fixture that the rule currently passes.
# Exit 0 = no hole, 1 = at least one hole. The last case (FP) is a copy-free graph the rule wrongly fails.
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d)"; trap 'rm -rf "$t"' EXIT
cd "$t" || exit 2; git init -q . && git config user.email a@b && git config user.name t
mkdir -p tools zenith-float-num/src/mantissa zenith-float-num/src/ops
echo 99 > tools/clone_count_baseline.txt; : > tools/rules_allowlist.txt
for f in helper_copy exception test_panic neg_via_mul; do echo 99 > "tools/${f}_count_baseline.txt"; done
export RULES_SCAN_ROOT="$t"
holes=0
probe() { # probe <want: fail|pass> <name> <ext.rs body of add_ref> <kern.rs content> [file:content extra]
  { echo 'impl N {'; echo "    pub fn add_ref(a: &N, b: &N) -> N { $3 }"
    for o in sub mul div; do echo "    pub fn ${o}_ref(a: &N, b: &N) -> N { N(0) }"; done; echo '}'; } > zenith-float-num/src/ext.rs
  printf '%b\n' "$4" > zenith-float-num/src/mantissa/kern.rs
  rm -f zenith-float-num/src/ops/h.rs; [ -z "${5:-}" ] || printf '%b\n' "${5#*:}" > "zenith-float-num/src/${5%%:*}"
  local out rc; out="$("$here/rules_scan.sh" --count-only 2>&1)"; rc=$?
  local got=pass; [ $rc -eq 0 ] || got=fail
  if [ "$got" = "$1" ]; then echo "ok    $2 ($got)"; else holes=$((holes + 1)); echo "HOLE  $2: want $1, got $got"; echo "$out" | grep BORROW | sed 's/^/        /'; fi
}
probe fail control-direct-clone 'let c = a.clone(); k(&c, b)' 'fn k(a: &N, b: &N) -> N { N(0) }'
probe fail ufcs-clone 'let c = Clone::clone(a); k(&c, b)' 'fn k(a: &N, b: &N) -> N { N(0) }'
probe fail space-before-paren 'let c = a.clone (); k(&c, b)' 'fn k(a: &N, b: &N) -> N { N(0) }'
probe fail slice-to_vec-path 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(<[u64]>::to_vec(&a.w).len()) }'
probe fail vec-from-slice 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(Vec::from(&a.w[..]).len()) }'
probe fail deref-map-collect 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|x| *x).collect(); N(v.len()) }'
probe fail extend-from-slice 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); v.extend_from_slice(&a.w); N(v.len()) }'
probe fail fn-as-value 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(a.w.iter().map(cp).count()) }\nfn cp(x: &u64) -> u64 { x.clone() }'
probe fail operator-dispatch 'a + b' 'impl core::ops::Add for &N { type Output = N; fn add(self, b: &N) -> N { self.clone() } }'
probe fail macro-body 'k(a, b)' 'macro_rules! dup { ($x:expr) => { $x.clone() } }\nfn k(a: &N, b: &N) -> N { dup!(a) }'
probe fail helper-outside-scope 'k(a, b)' 'fn k(a: &N, b: &N) -> N { crate::ops::h::deep(a) }' 'ops/h.rs:pub fn deep(a: &N) -> N { a.clone() }'
probe fail multiline-call 'k(a, b)' 'fn k(a: &N, b: &N) -> N {\n    a\n        .clone()\n}'
probe pass fp-multiline-block-comment 'k(a, b)' 'fn k(a: &N, b: &N) -> N {\n    /* note:\n       a.clone() is banned here */\n    N(0)\n}'
probe pass fp-name-collision 'k(a, b)' 'fn k(a: &N, b: &N) -> N { N(0) }' 'common_x.rs:x'
echo "borrow hole probe: $holes hole(s)"
[ "$holes" -eq 0 ]
