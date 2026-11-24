#!/usr/bin/env bash
# Red test for Jeff's ruling of 2026-11-23 07:43 CT on ROOT-CAUSE.md (clones 475 -> 1323 and test panics 3 -> 459 went
# unseen for a week because the library's rules_scan never ran): the rules are checked on every commit, push and gate.
#   1. one library scanner: only tools/rules_scan.{py,sh} (+ selftest) carry that name; it names itself first, and the
#      bound layer's scanner is zenith-float-bound/tools/bound_layer_scan.* and says it is NOT the library scan.
#   2. the gate's no-clone-progress step fails on any rise (tools/no_clone.sh --progress); a planted clone fails it.
#   3. rules-scan, no-clone-progress, fmt lead the gate (then target-dir, clippy); a planted clone or a crash fails rules-scan.
#   4. the committed pre-commit hook stops a planted clone, a count rise, a crash, red fmt and red clippy.
#   5. .github/workflows/rules-scan.yml runs the library rules and no_suppress_scan on every push; its commands fail a
#      planted clone.
#   6. CONTRIBUTING.md's per-commit checklist runs the library scan and says a 'clean' claim cites its log at that commit.
# usage: tools/test_rules_enforcement.sh [ROOT]   (ROOT: a checkout, default this one; text and fixture checks, no cargo)
set -u
root=${1:-$(cd "$(dirname "$0")/.." && pwd)}
fail=0; bad() { echo "FAIL: $*"; fail=1; }
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/test-rules-enforcement"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp (Spence n3)
t=$(mktemp -d "$zf_tmp/rules-enf.XXXXXX") || exit 2; trap 'rm -rf "$t"' EXIT
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

# 3. the library rules run first: rules-scan, no-clone-progress and fmt lead the gate, then the target-dir claim (it must
#    precede every build), then clippy; a kill after the first failure can no longer skip them. The rules-scan step runs
#    the root library scanner, a planted clone fails it, and a crashing check is a named FAIL there too.
order=$(grep -oE '^run [a-z0-9_-]+' "$root/tools/full_gate.sh" | awk '{print $2}' | head -5 | tr '\n' ' ')
[ "$order" = "rules-scan no-clone-progress fmt target-dir clippy-x86_64 " ] || bad "gate does not start with the rules steps: '$order'"
rs=$(step rules-scan)
[[ $rs == *'tools/rules_scan.sh '* && $rs != *zenith-float-bound* && $rs != *'cd '* ]] || bad "rules-scan step does not run the root library scanner: $rs"
grep -q 'test_rules_enforcement.sh' "$root/tools/full_gate.sh" || bad "tools/full_gate.sh does not run tools/test_rules_enforcement.sh"
# Spence R1 on 0a1cb6c: the zero-suppression scan and the pre-push hook's selftest run in the gate's unsafe-scan step.
us=$(step unsafe-scan)
for x in tools/no_suppress_selftest.sh tools/no_suppress_scan.sh tools/pre_push_selftest.sh tools/no_tmp_selftest.sh tools/no_tmp_scan.sh tools/no_debug_assert_selftest.sh tools/no_debug_assert_scan.sh; do
  [[ $us == *"$x"* ]] || bad "the gate's unsafe-scan step does not run $x"
done
g="$t/rs"; mkdir -p "$g/tools" "$g/zenith-float-num/src"
cp "$root"/tools/rules_scan.py "$root"/tools/rules_scan.sh "$root"/tools/rules_scan_selftest.sh "$root"/tools/no_clone.py "$g/tools/"
: > "$g/tools/rules_allowlist.txt"
for b in clone_count helper_copy_count exception_count test_panic_count neg_via_mul_count; do echo 0 > "$g/tools/${b}_baseline.txt"; done
printf 'fn a() {}\n' > "$g/zenith-float-num/src/x.rs"
( cd "$g" && git init -q . && git config user.email a@b && git config user.name t && git add -A && git commit -qm base ) || bad "fixture repo"
gb=$(git -C "$g" rev-parse HEAD)
printf 'fn a() {}\nfn b() -> u32 { 1 }\n' > "$g/zenith-float-num/src/x.rs"; ( cd "$g" && git add -A && git commit -qm clean )
(cd "$g" && RATCHET_BASE=$gb bash -c "$rs") > "$t/rs.out" 2>&1 || bad "rules-scan step failed a clean commit: $(tail -5 "$t/rs.out")"
printf 'fn a() {}\nfn b(x: &S) -> S { x.clone() }\n' > "$g/zenith-float-num/src/x.rs"; ( cd "$g" && git add -A && git commit -qm planted )
(cd "$g" && RATCHET_BASE=$gb bash -c "$rs") > "$t/rs.out" 2>&1 && bad "rules-scan step passed a planted clone"
grep -q 'zenith-float-num/src/x.rs:2 \[clone\]' "$t/rs.out" || bad "planted clone not reported by the rules-scan step: $(tail -5 "$t/rs.out")"
( cd "$g" && git reset -q --hard HEAD~1 )
(cd "$g" && RULES_SCAN_FAULT=borrow RATCHET_BASE=$gb bash -c "$rs") > "$t/rs.out" 2>&1 && bad "rules-scan step passed a crashing check"
grep -q '^rules_scan: FAIL check borrow CRASHED' "$t/rs.out" || bad "crash is not a named FAIL in the rules-scan step"
grep -q '^rules_scan: checks: clone=ok helper-copy=ok exceptions=ok test-panic=ok neg_via_mul=ok combined=ok borrow=CRASH' "$t/rs.out" \
  || bad "a crash hid the other checks in the rules-scan step: $(grep 'checks:' "$t/rs.out")"

# 4. the committed pre-commit hook (.githooks/pre-commit; core.hooksPath=.githooks) runs rules_scan --working, rules_scan
#    --count-only --no-rise HEAD, cargo fmt --check and cargo clippy -D warnings: a planted clone, a count rise, a crash
#    or a red fmt/clippy stops the commit. cargo is a stub here (its args are logged).
hk="$root/.githooks/pre-commit"
if [ -x "$hk" ]; then
  h="$t/hk"; mkdir -p "$h/tools" "$h/zenith-float-num/src" "$h/.githooks" "$t/bin"
  cp "$root"/tools/rules_scan.py "$root"/tools/rules_scan.sh "$root"/tools/no_clone.py "$h/tools/"; cp "$hk" "$h/.githooks/"
  : > "$h/tools/rules_allowlist.txt"
  for b in clone_count helper_copy_count exception_count test_panic_count neg_via_mul_count; do echo 0 > "$h/tools/${b}_baseline.txt"; done
  printf 'fn a() {}\n' > "$h/zenith-float-num/src/x.rs"
  printf '#!/usr/bin/env bash\necho "$*" >> "%s/cargo.args"\n[ -n "${STUB_CARGO_FAIL:-}" ] && [[ "$*" == "$STUB_CARGO_FAIL"* ]] && exit 1\nexit 0\n' "$t" > "$t/bin/cargo"; chmod +x "$t/bin/cargo"
  hc() { (cd "$h" && PATH="$t/bin:$PATH" git commit -qm "$1") > "$t/hk.out" 2>&1; }
  ( cd "$h" && git init -q . && git config user.email a@b && git config user.name t && git add -A && git commit -qm base --no-verify \
    && git config core.hooksPath .githooks ) || bad "hook fixture repo"
  printf 'fn a() {}\nfn b() -> u32 { 1 }\n' > "$h/zenith-float-num/src/x.rs"; (cd "$h" && git add -A)
  hc clean || bad "hook refused a clean commit: $(tail -5 "$t/hk.out")"
  grep -qx 'fmt --all -- --check' "$t/cargo.args" || bad "hook did not run cargo fmt --all -- --check"
  grep -qx 'clippy --workspace --all-targets --all-features -- -D warnings' "$t/cargo.args" || bad "hook did not run cargo clippy -D warnings"
  n0=$(git -C "$h" rev-list --count HEAD)
  printf 'fn a() {}\nfn b() -> u32 { 1 }\nfn c(x: &S) -> S { x.clone() }\n' > "$h/zenith-float-num/src/x.rs"; (cd "$h" && git add -A)
  hc planted && bad "hook let a planted clone through"
  grep -q 'zenith-float-num/src/x.rs:3 \[clone\]' "$t/hk.out" || bad "hook did not report the planted clone: $(tail -5 "$t/hk.out")"
  grep -q 'clone count 1 is ABOVE the baseline 0 and RISES vs HEAD (0)' "$t/hk.out" || bad "hook did not report the count rise: $(tail -8 "$t/hk.out")"
  [ "$(git -C "$h" rev-list --count HEAD)" = "$n0" ] || bad "a commit was created although the hook failed"
  printf 'fn a() {}\nfn b() -> u32 { 1 }\n' > "$h/zenith-float-num/src/x.rs"; printf '// edit\n' >> "$h/zenith-float-num/src/x.rs"; (cd "$h" && git add -A)
  RULES_SCAN_FAULT=test-panic hc crash && bad "hook let a crashing rules_scan check through"
  grep -q 'FAIL check test-panic CRASHED' "$t/hk.out" || bad "hook did not name the crash: $(tail -5 "$t/hk.out")"
  STUB_CARGO_FAIL=fmt hc fmtred && bad "hook let a red cargo fmt through"
  STUB_CARGO_FAIL=clippy hc clippyred && bad "hook let a red clippy through"
  hc green || bad "hook refused a clean edit: $(tail -5 "$t/hk.out")"
else
  bad "no executable .githooks/pre-commit"
fi
grep -q 'core.hooksPath .githooks' "$root/CONTRIBUTING.md" || bad "CONTRIBUTING.md does not say how to enable the hook (git config core.hooksPath .githooks)"

# 5. GitHub Actions runs the library rules on every push (all branches) and PR: .github/workflows/rules-scan.yml. Its own
#    `run:` commands are executed in the part-3 fixture: they pass on the clean tree and fail on a planted clone.
wf="$root/.github/workflows/rules-scan.yml"
if [ -f "$wf" ]; then
  grep -qE "^  push:$" "$wf" && grep -qE "^    branches: \[ '\*\*' \]$" "$wf" || bad "rules-scan.yml does not run on a push to every branch"
  grep -qE "^  pull_request:$" "$wf" || bad "rules-scan.yml does not run on pull requests"
  grep -q 'fetch-depth: 0' "$wf" || bad "rules-scan.yml checkout is shallow (the range scan needs history)"
  mapfile -t runs < <(sed -nE 's/^ +- run: (.*)$/\1/p' "$wf")
  [ ${#runs[@]} -ge 3 ] || bad "rules-scan.yml has fewer than 3 run steps"
  printf '%s\n' "${runs[@]}" | grep -q 'tools/rules_scan.sh ' || bad "rules-scan.yml does not run tools/rules_scan.sh"
  printf '%s\n' "${runs[@]}" | grep -q 'no_clone.sh --progress' || bad "rules-scan.yml does not run tools/no_clone.sh --progress"
  printf '%s\n' "${runs[@]}" | grep -q 'zenith-float-bound' && bad "rules-scan.yml runs the bound scanner as the rules scan"
  cp "$root"/tools/no_clone.sh "$root"/tools/no_clone_selftest.sh "$g/tools/"
  printf '%s\n' "${runs[@]}" | grep -q 'tools/no_suppress_scan.sh' || bad "rules-scan.yml does not run tools/no_suppress_scan.sh"
  cp "$root"/tools/no_suppress_scan.sh "$root"/tools/no_suppress_scan.py "$root"/tools/no_suppress_selftest.sh "$g/tools/"
  printf '# path\tcount\tnote\n' > "$g/tools/unsafe_pending.tsv"
  printf 'clone()\t0\nto_vec()\t0\ncloned()\t0\nClone::clone\t0\nclone_from\t0\nimpl Clone\t0\nderive(Clone)\t0\ntotal\t0\nhelper_copies\t0\n' > "$g/tools/no_clone_progress_baseline.txt"
  ( cd "$g" && git add -A && git commit -qm ci-tools )
  wfrun() { local r=0 c; for c in "${runs[@]}"; do (cd "$g" && RATCHET_BASE=$gb bash -o pipefail -c "$c") >> "$t/wf.out" 2>&1 || r=1; done; return $r; }
  : > "$t/wf.out"; wfrun || bad "rules-scan.yml commands failed a clean tree: $(tail -5 "$t/wf.out")"
  printf 'fn a() {}\nfn b(x: &S) -> S { x.clone() }\n' > "$g/zenith-float-num/src/x.rs"; ( cd "$g" && git add -A && git commit -qm planted-ci )
  : > "$t/wf.out"; wfrun && bad "rules-scan.yml commands passed a planted clone"
else
  bad "no .github/workflows/rules-scan.yml"
fi

# 6. the per-commit checklist names the library rules_scan and the 'clean' evidence rule
c="$root/CONTRIBUTING.md"
grep -q '^## Per-commit checklist' "$c" || bad "CONTRIBUTING.md has no per-commit checklist"
grep -q 'tools/rules_scan.sh --count-only' "$c" || bad "the checklist does not run the library tools/rules_scan.sh --count-only"
grep -q "Every 'clean' claim cites a log from the library's own \`tools/rules_scan.sh\` at that exact commit" "$c" || bad "the checklist lacks the 'clean' evidence rule"

[ $fail -eq 0 ] && echo "test_rules_enforcement: ok"
exit $fail
