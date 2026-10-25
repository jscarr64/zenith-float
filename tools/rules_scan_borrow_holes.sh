#!/usr/bin/env bash
# Probe (Harper, for Rusty): NEW holes in the 581cdfe BORROW rule — copies reachable from add_ref that it passes.
# NOT wired into the selftest; once a case is caught, move it into rules_scan_selftest.sh (hx) and drop it here.
set -u
R="$(cd "$(dirname "$0")/.." && pwd)"; t="$(mktemp -d)"; trap 'rm -rf "$t"' EXIT
cd "$t" || exit 2; git init -q . && git config user.email a@b && git config user.name t
mkdir -p tools zenith-float-num/src/mantissa; echo 99 > tools/clone_count_baseline.txt; : > tools/rules_allowlist.txt
for f in helper_copy exception test_panic neg_via_mul; do echo 99 > "tools/${f}_count_baseline.txt"; done
export RULES_SCAN_ROOT="$t"; holes=0
v() { # v <name> <add_ref body> <kern.rs>
  { echo 'impl N {'; echo "    pub fn add_ref(a: &N, b: &N) -> N { $2 }"
    for o in sub mul div; do echo "    pub fn ${o}_ref(a: &N, b: &N) -> N { N(0) }"; done; echo '}'; } > zenith-float-num/src/ext.rs
  printf '%b\n' "$3" > zenith-float-num/src/mantissa/kern.rs
  out="$("$R/tools/rules_scan.sh" --count-only 2>&1)" && { holes=$((holes + 1)); echo "NEW HOLE  $1"; printf '%b\n' "    add_ref: $2\n    kern.rs:\n$3" | sed 's/^/      /'; } \
    || echo "caught    $1: $(echo "$out" | grep 'rules_scan: BORROW' | head -1 | sed 's/rules_scan: //')"
}
KK='fn k(a: &N, b: &N) -> N { N(0) }'
v ref-pattern-closure 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|&x| x).collect(); N(v.len()) }'
v copy-from-slice 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = vec![0u64; a.w.len()]; v.copy_from_slice(&a.w); N(v.len()) }'
v slice-into-vec 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.as_slice().into(); N(v.len()) }'
v clone-into 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); a.w.as_slice().clone_into(&mut v); N(v.len()) }'
v vec-macro-copy 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v = vec![a.w[0]; 4]; N(v.len()) }'
v from-iter-deref 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let v: Vec<u64> = a.w.iter().map(|x: &u64| -> u64 { *x }).collect(); N(v.len()) }'
v extend-iter-copied 'k(a, b)' 'fn k(a: &N, b: &N) -> N { let mut v = Vec::new(); v.extend(a.w.iter()); N(v.len()) }'
v generic-helper 'k(a, b)' 'fn k(a: &N, b: &N) -> N { g::<N>(a) }\nfn g<T: Clone>(x: &T) -> T { x.clone() }'
echo "variants: $holes new hole(s)"
[ "$holes" -eq 0 ]
