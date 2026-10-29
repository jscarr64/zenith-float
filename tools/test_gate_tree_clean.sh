#!/usr/bin/env bash
# Red test: a gate run must leave the worktree as it found it (`git status --porcelain` unchanged). The
# prepublish step's compare-bench used to rewrite the tracked doc/compare-results.tsv; it now writes to the
# run's log dir ($T), and only tools/regen_compare_results.sh writes doc/ on purpose. Runs the gate's
# prepublish and tree-clean steps on a scratch clone (this checkout's working files overlaid) with stub
# cargo/rustup; never touches this worktree.
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY CARGO_BUILD_JOBS RUST_TEST_THREADS \
  GATE_STEPS COMPARE_OUT
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
git clone -q "$root" "$tmp/repo"
git -C "$root" ls-files -m -o --exclude-standard | while IFS= read -r f; do
  [[ -f "$root/$f" ]] && { mkdir -p "$(dirname "$tmp/repo/$f")"; cp -p "$root/$f" "$tmp/repo/$f"; }
done
git -C "$tmp/repo" add -A && git -C "$tmp/repo" -c user.name=selftest -c user.email=selftest@invalid commit -q --allow-empty -m overlay
mkdir -p "$tmp/bin"
cat > "$tmp/bin/cargo" <<'SH'
#!/usr/bin/env bash
echo "$*" >> "${0%/bin/cargo}/cargo.args"
case "$*" in audit\ --version) echo "cargo-audit 0.0"; exit 0;; esac
if [[ "$1" == test && "$*" == *--message-format=json* ]]; then
  b="$CARGO_TARGET_DIR/debug/deps/stub-0"; mkdir -p "${b%/*}"; : > "$b"
  printf '{"reason":"compiler-artifact","executable":"%s","profile":{"test":true}}\n' "$b"; exit 0
fi
if [[ "$1" == run ]]; then printf 'library\tprecision_bits\ttask\tbatch_size\tbest_time_us\nstub\t132\tadd\t1\t1\n'; exit 0; fi
if [[ "$1" == test ]]; then echo "test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured"; fi
exit 0
SH
cat > "$tmp/bin/rustup" <<'SH'
#!/usr/bin/env bash
case "$*" in *rustc\ --version*) echo "rustc 1.99.1 (stub)";; esac
exit 0
SH
chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
rc=0; out=$(HOME="$tmp" PATH="$tmp/bin:$PATH" POPLOG_ROOT="$tmp/pop" T="$tmp/gate" REPO="$tmp/repo" \
  GATE_STEPS=prepublish,tree-clean bash "$tmp/repo/tools/full_gate.sh" treeclean 2>&1) || rc=$?
log=$(find "$tmp/pop" -name '*.log' 2>/dev/null | head -1 || true)
fail=0
grep -q '^run -q -p zenith-float-compare' "$tmp/cargo.args" || { echo "FAIL: prepublish never reached compare-bench (test proves nothing)"; fail=1; }
dirty=$(git -C "$tmp/repo" status --porcelain)
[[ -z "$dirty" ]] || { echo "FAIL: gate run left the tree dirty: $dirty"; fail=1; }
grep -q 'STEP END tree-clean verdict=PASS' "$log" || { echo "FAIL: no passing tree-clean step: $(grep 'tree-clean' "$log" | tail -2)"; fail=1; }
[[ -s "$tmp/gate/compare-results.tsv" ]] || { echo "FAIL: compare-bench output is not in the run's log dir \$T"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: gate run leaves git status --porcelain unchanged; compare-results.tsv in \$T"
exit "$fail"
