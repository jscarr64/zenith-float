#!/usr/bin/env bash
# Red test: gate scratch never lands in /tmp or a bare mktemp, and is removed. (1) tools/fresh_test_run.sh keeps
# its build scratch inside its own target dir and removes it on a passing, failing and --no-run run (the old
# `exec` at the end dropped the cleanup trap and leaked $TMPDIR/tmp.* dirs). (2) tools/selftest_scratch.sh
# hands out dirs only under SELFTEST_SCRATCH: unset or under /tmp is an error. (3) every gate selftest takes
# its scratch from it (no bare `mktemp -d`), and full_gate.sh exports SELFTEST_SCRATCH under its log dir.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
s=${SELFTEST_SCRATCH:?SELFTEST_SCRATCH must be set (a dir outside /tmp)}
mkdir -p "$s"; t=$(mktemp -d "$s/selftest_scratch.XXXXXX"); trap 'rm -rf "$t"' EXIT
fail=0
mkdir -p "$t/bin" "$t/tmpdir"; git init -q "$t/wt"
cat > "$t/bin/cargo" <<'SH'
#!/usr/bin/env bash
if [[ "$*" == *--message-format=json* ]]; then
  b="$CARGO_TARGET_DIR/debug/deps/stub-0"; mkdir -p "${b%/*}"; : > "$b"
  printf '{"reason":"compiler-artifact","executable":"%s","profile":{"test":true}}\n' "$b"; exit 0
fi
if [[ -n "${STUB_SLEEP_PID:-}" ]]; then echo $$ > "$STUB_SLEEP_PID"; exec sleep 30; fi
echo "test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured"; exit "${STUB_RC:-0}"
SH
chmod +x "$t/bin/cargo"
for c in pass:0 fail:101 norun:0; do
  n=${c%%:*} want=${c#*:} rc=0 extra=()
  [[ $n == norun ]] && extra=(--no-run)
  (cd "$t/wt" && env PATH="$t/bin:$PATH" TMPDIR="$t/tmpdir" CARGO_TARGET_DIR="$t/target-$n" STUB_RC="$want" \
    "$root/tools/fresh_test_run.sh" cargo test "${extra[@]}" -- --x > "$t/out-$n" 2>&1) || rc=$?
  [[ $rc -eq $want ]] || { echo "FAIL[ftr-$n]: exit $rc, want $want"; sed 's/^/  | /' "$t/out-$n"; fail=1; }
  [[ -z $(ls -A "$t/tmpdir") ]] || { echo "FAIL[ftr-$n]: scratch left in \$TMPDIR: $(ls -A "$t/tmpdir" | tr '\n' ' ')"; fail=1; rm -rf "${t:?}/tmpdir/"*; }
  [[ -z $(find "$t/target-$n" -maxdepth 1 -name '.fresh.*') ]] || { echo "FAIL[ftr-$n]: .fresh.* scratch left in the target dir"; fail=1; }
done
# a step timeout (TERM) reaches the test command too: no orphan, scratch removed
rc=0; (cd "$t/wt" && env PATH="$t/bin:$PATH" TMPDIR="$t/tmpdir" CARGO_TARGET_DIR="$t/target-term" STUB_SLEEP_PID="$t/pid" \
  timeout 3 "$root/tools/fresh_test_run.sh" cargo test > "$t/out-term" 2>&1) || rc=$?
sleep 1
[[ $rc -ne 0 ]] || { echo "FAIL[ftr-term]: timed-out run exited 0"; fail=1; }
[[ -s "$t/pid" ]] && ! kill -0 "$(cat "$t/pid")" 2> /dev/null || { echo "FAIL[ftr-term]: test command outlived the timed-out fresh_test_run (orphan)"; fail=1; kill "$(cat "$t/pid")" 2> /dev/null || true; }
[[ -z $(find "$t/target-term" -maxdepth 1 -name '.fresh.*') ]] || { echo "FAIL[ftr-term]: scratch left after TERM"; fail=1; }
grep -q '^exec ' "$root/tools/fresh_test_run.sh" && { echo "FAIL[ftr-exec]: fresh_test_run.sh still ends in exec (drops its cleanup trap)"; fail=1; }
# (2) the scratch helper
sc="$root/tools/selftest_scratch.sh"
d=$(SELFTEST_SCRATCH="$t/sc" bash "$sc" probe 2>/dev/null) && [[ $d == "$t/sc/probe."* && -d $d ]] || { echo "FAIL[helper]: no dir under SELFTEST_SCRATCH: '${d:-}'"; fail=1; }
env -u SELFTEST_SCRATCH bash "$sc" probe > /dev/null 2>&1 && { echo "FAIL[helper-unset]: unset SELFTEST_SCRATCH accepted"; fail=1; }
for bad in /tmp /tmp/x /tmp/../tmp/y; do
  SELFTEST_SCRATCH=$bad bash "$sc" probe > /dev/null 2>&1 && { echo "FAIL[helper-tmp]: SELFTEST_SCRATCH=$bad accepted"; fail=1; }
done
# (3) no gate selftest makes a bare mktemp dir; full_gate exports SELFTEST_SCRATCH under $T
for f in "$root"/tools/*selftest*.sh "$root"/tools/test_*.sh; do
  [[ $f == "$root/tools/test_selftest_scratch.sh" ]] && continue
  grep -nE 'mktemp( -d)?\)|mktemp -d"|mktemp -d \$\{TMPDIR' "$f" | sed "s|^|FAIL[bare-mktemp]: ${f#"$root"/}:|" | grep . && fail=1
done
grep -q 'export SELFTEST_SCRATCH=\${SELFTEST_SCRATCH:-\$T/scratch}' "$root/tools/full_gate.sh" || { echo "FAIL[gate]: full_gate.sh does not export SELFTEST_SCRATCH under \$T"; fail=1; }
# (4) the gate's own temp files (python tempfile, mktemp) go under its log dir, and a step timeout (TERM) on
#     rules_scan.py still removes its TemporaryDirectory (it leaked one per timed-out run)
grep -q '^export TMPDIR=$T/tmp' "$root/tools/full_gate.sh" || { echo "FAIL[gate-tmpdir]: full_gate.sh does not export TMPDIR under \$T"; fail=1; }
mkdir -p "$t/py"
(cd "$root" && TMPDIR="$t/py" timeout -s TERM 4 tools/rules_scan.sh c74e834^..HEAD > "$t/out-py" 2>&1) || true
[[ -z $(ls -A "$t/py") ]] || { echo "FAIL[rules-scan-term]: TERM left rules_scan.py temp dirs: $(ls -A "$t/py" | tr '\n' ' ')"; fail=1; }
[[ $fail -eq 0 ]] && echo "ok: fresh_test_run cleans its in-target scratch (pass/fail/no-run); selftests use SELFTEST_SCRATCH, never /tmp"
exit "$fail"
