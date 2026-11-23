#!/usr/bin/env bash
# Red test for Jeff's ruling of 2026-11-23 07:43 CT on ROOT-CAUSE.md (clones 475 -> 1323 and test panics 3 -> 459 went
# unseen for a week because the library's rules_scan never ran): the rules are checked on every commit, push and gate.
#   1. one library scanner: only tools/rules_scan.{py,sh} (+ selftest) carry that name; it names itself first, and the
#      bound layer's scanner is zenith-float-bound/tools/bound_layer_scan.* and says it is NOT the library scan.
#   2. the gate's no-clone-progress step fails on any rise (tools/no_clone.sh --progress); a planted clone fails it.
# usage: tools/test_rules_enforcement.sh [ROOT]   (ROOT: a checkout, default this one; text and fixture checks, no cargo)
set -u
root=${1:-$(cd "$(dirname "$0")/.." && pwd)}
fail=0; bad() { echo "FAIL: $*"; fail=1; }
t=$(mktemp -d "${TMPDIR:-/tmp}/rules-enf.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/fx/tools"; : > "$t/empty.diff"

# 1. one library scanner, and each scanner names itself
names=$(cd "$root" && find . -path ./target -prune -o -path ./.git -prune -o -name __pycache__ -prune -o -name 'rules_scan*' -print | sed 's#^\./##' | LC_ALL=C sort | tr '\n' ' ')
[ "$names" = "tools/rules_scan.py tools/rules_scan.sh tools/rules_scan_selftest.sh " ] \
  || bad "files named rules_scan* other than the library scanner: '$names'"
first=$(RULES_SCAN_ROOT="$t/fx" python3 -B "$root/tools/rules_scan.py" --diff-file "$t/empty.diff" 2>&1 | head -1)
[[ $first == 'rules_scan: zenith-float LIBRARY scanner (tools/rules_scan.py'* ]] || bad "library scanner does not name itself first: '$first'"
b="$root/zenith-float-bound/tools/bound_layer_scan.py"
if [ -f "$b" ]; then
  first=$(RULES_SCAN_ROOT="$t/fx" python3 -B "$b" --diff-file "$t/empty.diff" 2>&1 | head -1)
  [[ $first == 'bound_layer_scan: zenith-float-bound LAYER scanner'*'NOT the library scan'* ]] || bad "bound scanner does not say it is not the library scan: '$first'"
else
  bad "no zenith-float-bound/tools/bound_layer_scan.py"
fi

# 2. no-clone-progress FAILS on any rise (it ended in `true` since da53be7): run the gate's own step command in a fixture
step() { grep -m1 -E "^run $1 \"" "$root/tools/full_gate.sh" | sed -E "s/^run $1 \"(.*)\"\$/\1/"; }
nc=$(step no-clone-progress)
[[ -n $nc ]] || bad "no no-clone-progress step in tools/full_gate.sh"
[[ $nc == *'true'* ]] && bad "no-clone-progress still contains \`true\`: $nc"
[[ $nc == *'no_clone.sh --progress'* ]] || bad "no-clone-progress does not run tools/no_clone.sh --progress: $nc"
f="$t/nc"; mkdir -p "$f/tools" "$f/zenith-float-num/src"
cp "$root"/tools/no_clone.py "$root"/tools/no_clone.sh "$root"/tools/no_clone_selftest.sh "$f/tools/" 2>/dev/null
printf 'clone()\t1\nto_vec()\t0\ncloned()\t0\nClone::clone\t0\nclone_from\t0\nimpl Clone\t0\nderive(Clone)\t0\ntotal\t1\nhelper_copies\t0\n' > "$f/tools/no_clone_progress_baseline.txt"
printf 'fn f(x: &S) { let y = x.clone(); }\n' > "$f/zenith-float-num/src/a.rs"
(cd "$f" && bash -c "$nc") > "$t/nc.out" 2>&1 || bad "no-clone-progress failed a tree at its baseline: $(tail -3 "$t/nc.out")"
printf 'fn f(x: &S) { let y = x.clone(); let z = x.clone(); }\n' > "$f/zenith-float-num/src/a.rs"
(cd "$f" && bash -c "$nc") > "$t/nc.out" 2>&1 && bad "no-clone-progress passed a planted clone (1 -> 2): $(tail -3 "$t/nc.out")"

[ $fail -eq 0 ] && echo "test_rules_enforcement: ok"
exit $fail
