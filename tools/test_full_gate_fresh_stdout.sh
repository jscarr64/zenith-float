#!/usr/bin/env bash
# Red test (2026-10-23): tools/full_gate.sh run under fresh_test_run.sh (CMD... wrapper) must forward each step's
# cargo "Running ... (PATH)" and "test result:" lines to stdout, so the wrapper sees the binaries it ran instead of
# reporting "no test binary ran". Stub cargo prints a Running line for a real file inside a throwaway target dir.
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY CARGO_BUILD_JOBS RUST_TEST_THREADS
root=$(cd "$(dirname "$0")/.." && pwd)
wrap=${FRESH_TEST_RUN:-$HOME/pop-logs/bin/fresh_test_run.sh}
[[ -x "$wrap" ]] || { echo "FAIL: fresh_test_run.sh not found at $wrap (set FRESH_TEST_RUN)"; exit 1; }
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin" "$tmp/target/debug/deps"
printf '#!/bin/sh\nexit 0\n' > "$tmp/target/debug/deps/fixture-0123abcd"; chmod +x "$tmp/target/debug/deps/fixture-0123abcd"
git -C "$root" rev-parse --show-toplevel > "$tmp/target/.worktree-owner"  # claimed by this worktree, as a fresh gate dir would be
cat > "$tmp/bin/cargo" <<'SH'
#!/usr/bin/env bash
case "$*" in audit\ --version) echo "cargo-audit 0.0"; exit 0;; esac
# tools/fresh_test_run.sh (full_gate's cargo test steps, 6946fbf) builds with --no-run --message-format=json first:
# report the fixture binary inside CARGO_TARGET_DIR, then the real run prints its Running/test result lines
if [[ "$1" == test && "$*" == *--message-format=json* ]]; then
  printf '{"reason":"compiler-artifact","executable":"%s","profile":{"test":true}}\n' "$CARGO_TARGET_DIR/debug/deps/fixture-0123abcd"; exit 0
fi
if [[ "$1" == test ]]; then
  echo "     Running unittests src/lib.rs ($CARGO_TARGET_DIR/debug/deps/fixture-0123abcd)" >&2
  echo "test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out"
fi
exit 0
SH
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"; chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
out=$(cd "$root" && HOME="$tmp" PATH="$tmp/bin:$PATH" CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" \
  POPLOG_ROOT="$tmp/pop" T="$tmp/gate" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s \
  "$wrap" bash tools/full_gate.sh freshtest 2>&1) || true
fail=0; [[ -n "${DEBUG:-}" ]] && echo "$out" | tail -15
grep -q 'no test binary ran' <<<"$out" && { echo "FAIL: fresh_test_run saw no test binary (full_gate kept Running lines in \$T only)"; fail=1; }
grep -q 'FOREIGN' <<<"$out" && { echo "FAIL: wrapper reported a FOREIGN binary"; fail=1; }
grep -qE '^ *Running .*\(/.*fixture-0123abcd\)$' <<<"$out" || { echo "FAIL: no Running line forwarded to stdout"; fail=1; }
grep -q '^test result: ok. 1 passed' <<<"$out" || { echo "FAIL: no test result line forwarded to stdout"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: full_gate.sh forwards Running/test result lines to fresh_test_run.sh"
exit "$fail"
