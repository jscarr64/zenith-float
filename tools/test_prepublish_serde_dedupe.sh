#!/usr/bin/env bash
# Red test (Jeff 2026-11-24 9:06 AM CT: "Go: in the gate, prepublish skips its random,serde rerun only after checking
# features-random-serde passed rc=0 on the same SHA; standalone prepublish stays complete").
# Drives the real full_gate.sh `run features-random-serde` / `run prepublish` lines and the real scripts/zenith_prepublish.sh
# in a throwaway git repo with stub cargo / fresh_test_run (no real build), and checks:
#   (a) gate mode: prepublish skips `cargo test --workspace --features random,serde` only when the gate's own serde step
#       left evidence of rc=0 on this exact HEAD with a clean tree;
#   (b) no evidence, rc != 0, another SHA, a dirty tree (at the step or at prepublish), or HEAD moving during the step
#       => prepublish runs the full serde suite;
#   (c) standalone prepublish (no gate evidence variable) always runs everything;
#   and every other prepublish part runs in every case.
# usage: tools/test_prepublish_serde_dedupe.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
fail=0
bad() { echo "FAIL: $*"; fail=1; }
SERDE='cargo test --workspace --features random,serde'

zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/test-prepublish-serde-dedupe"; mkdir -p "$zf_tmp" || exit 2  # job scratch, never /tmp (Spence n3)
w=$(mktemp -d "$zf_tmp/prepub-dedupe.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
repo="$w/repo"; mkdir -p "$repo/tools" "$repo/scripts" "$repo/zenith-float-num/src" "$w/bin" "$w/target"
cp "$root/scripts/zenith_prepublish.sh" "$repo/scripts/"
[[ -f "$root/scripts/gate_serde_evidence.sh" ]] && cp "$root/scripts/gate_serde_evidence.sh" "$repo/scripts/"
gate_serde="$(grep -E '^run features-random-serde ' "$root/tools/full_gate.sh")" || bad "no run features-random-serde line in tools/full_gate.sh"
gate_prepub="$(grep -E '^run prepublish ' "$root/tools/full_gate.sh")" || bad "no run prepublish line in tools/full_gate.sh"
[[ $fail -eq 0 ]] || { echo "test_prepublish_serde_dedupe: FAIL"; exit 1; }
# the gate names its run once (ZF_GATE_RUN_ID) so evidence from another run, or written by hand, never counts
grep -qE '^export ZF_GATE_RUN_ID=' "$root/tools/full_gate.sh" || bad "tools/full_gate.sh does not export ZF_GATE_RUN_ID"

REC="$w/rec.txt"; export REC
rec_stub() { printf '#!/usr/bin/env bash\necho "%s $*" >> "$REC"\n%s\n' "$2" "${3:-exit 0}" > "$1"; chmod +x "$1"; }
# fresh_test_run.sh stub: records the command, exits STUB_RC (default 0)
rec_stub "$repo/tools/fresh_test_run.sh" fresh_test_run '[[ -n "${STUB_COMMIT:-}" ]] && echo m > f-during && git add -A && git -c user.name=t -c user.email=t@t commit -qm during
exit "${STUB_RC:-0}"'
rec_stub "$w/bin/cargo" cargo
printf '' > "$repo/scripts/use_rust_1_99.sh"
for s in compare-bench.sh ci_no_intrinsics.sh check_precision_docs.sh check_leaves.sh; do rec_stub "$repo/scripts/$s" "$s"; done
rec_stub "$repo/scripts/prepublish_lib.py" prepublish_lib 'exit 0'
printf 'import os, sys\nopen(os.environ["REC"], "a").write("prepublish_lib %%s\\n" %% " ".join(sys.argv[1:]))\n' > "$repo/scripts/prepublish_lib.py"
echo 'pub const PROPTEST_CASES: u32 = 1000;' > "$repo/zenith-float-num/src/defs.rs"
git -C "$repo" init -q && git -C "$repo" add -A && git -C "$repo" -c user.name=t -c user.email=t@t commit -qm base

export PATH="$w/bin:$PATH" CARGO_TARGET_DIR="$w/target" RUST_TEST_THREADS=4 ZF_GATE_RUN_ID=selftest-$$-$RANDOM
T="$w/gate"; tag=t; TT="-- --test-threads=4"; mkdir -p "$T"
run() { local name=$1; shift; (cd "$repo" && bash -o pipefail -c "$*") > "$T/$tag-$name.log" 2>&1; echo "$name=$?" >> "$T/$tag.rc"; }
serde_step() { (export STUB_RC=${1:-0}; eval "$gate_serde"); }
prepub_step() { : > "$REC"; eval "$gate_prepub"; }
standalone() { : > "$REC"; (cd "$repo" && env -u ZF_GATE_RUN_ID bash scripts/zenith_prepublish.sh) > "$T/$tag-standalone.log" 2>&1; echo "standalone=$?" >> "$T/$tag.rc"; }
commit_new() { echo "$1" > "$repo/f-$1"; git -C "$repo" add -A && git -C "$repo" -c user.name=t -c user.email=t@t commit -qm "$1"; }
# expect CASE skip|run LOG: the serde rerun was (not) run, every other prepublish part ran, prepublish exited 0
expect() {
  local case=$1 want=$2 log=$3 n
  n=$(grep -c -- "fresh_test_run $SERDE" "$REC")
  if [[ $want == skip ]]; then
    [[ $n -eq 0 ]] || bad "$case: prepublish reran the random,serde suite"
    grep -q 'SKIP .*random,serde' "$log" || bad "$case: no SKIP line in the prepublish log"
  else
    [[ $n -eq 1 ]] || bad "$case: prepublish did not run the random,serde suite (ran $n times)"
    grep -q 'SKIP .*random,serde' "$log" && bad "$case: prepublish logged a SKIP"
  fi
  grep -q 'fresh_test_run cargo test -p zenith-float-num --features mpfr-tests --release' "$REC" || bad "$case: mpfr-tests not run"
  grep -q 'fresh_test_run cargo test -p zenith-float-num --test proptest_props' "$REC" || bad "$case: proptest not run"
  grep -q 'cargo build --no-default-features' "$REC" || bad "$case: no_std build not run"
  for s in compare-bench.sh ci_no_intrinsics.sh check_precision_docs.sh check_leaves.sh 'prepublish_lib version' 'prepublish_lib verify' 'prepublish_lib hygiene'; do
    grep -q "^$s" "$REC" || bad "$case: $s not run"
  done
  grep -q 'prepublish ok (11/11)' "$log" || bad "$case: prepublish did not finish 11/11"
}
plog="$T/$tag-prepublish.log"

# (a) serde step rc=0, same SHA, clean tree => skip
serde_step 0
compgen -G "$T/*evidence*" > /dev/null || bad "a: the gate serde step wrote no evidence file"
prepub_step; expect "a: rc=0 same SHA clean" skip "$plog"
grep -q "$(git -C "$repo" rev-parse HEAD)" "$plog" || bad "a: the SKIP line does not name the SHA"
# (c) standalone, even with good evidence on disk => full
standalone; expect "c: standalone" run "$T/$tag-standalone.log"
# (b) serde step failed => full
serde_step 1; prepub_step; expect "b: rc=1" run "$plog"
# (b) no evidence (gate serde step never ran / was killed before writing) => full
serde_step 0; rm -f "$T"/*.evidence "$T"/*evidence*; prepub_step; expect "b: no evidence" run "$plog"
# (b) a new commit between the serde step and prepublish => full
serde_step 0; commit_new moved; prepub_step; expect "b: other SHA" run "$plog"
# (b) dirty tree at prepublish time => full
serde_step 0; echo dirty > "$repo/f-moved"; prepub_step; expect "b: dirty at prepublish" run "$plog"
git -C "$repo" checkout -q -- f-moved
# (b) dirty tree while the serde step ran (clean again by prepublish) => full
echo x > "$repo/untracked.rs"; serde_step 0; rm -f "$repo/untracked.rs"; prepub_step; expect "b: dirty at serde step" run "$plog"
# (b) HEAD moved while the serde step ran => full
(export STUB_RC=0 STUB_COMMIT=1; eval "$gate_serde"); prepub_step; expect "b: HEAD moved during serde step" run "$plog"
# (b) evidence not written by this gate run's serde step: hand-made (no run id) or from another run (other id) => full
serde_step 0
for ev in "$T"/*evidence*; do
  [[ -f "$ev" ]] || continue
  printf 'sha=%s\nrc=0\nclean=yes\n' "$(git -C "$repo" rev-parse HEAD)" > "$ev"; prepub_step; expect "b: hand-made evidence" run "$plog"
  serde_step 0; (export ZF_GATE_RUN_ID=other-run; prepub_step); expect "b: evidence from another gate run" run "$plog"
  forged=1
done
[[ -n "${forged:-}" ]] || bad "b: no evidence file to forge"
# and a good gate run again skips (state not sticky)
serde_step 0; prepub_step; expect "a: rc=0 again" skip "$plog"

[[ $fail -eq 0 ]] && { echo "test_prepublish_serde_dedupe: ok"; exit 0; }
echo "test_prepublish_serde_dedupe: FAIL"; exit 1
