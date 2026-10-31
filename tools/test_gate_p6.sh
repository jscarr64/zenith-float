#!/usr/bin/env bash
# Red test (P6 ruling 1): (a) tools/gate_no_32bit.sh FAILs on any 32-bit target in the gate scripts
# (tools/, scripts/), passes on a clean tree, honours its allowlist; the real tree is clean;
# (b) the gate has a build-nostd-none step building zenith-float and zenith-float-num with
# --no-default-features --lib for x86_64-unknown-none, and a no-32bit-scripts step, both PASS under stub cargo;
# (c) ci_nostd.sh, ci.sh and zenith_prepublish.sh build for x86_64-unknown-none.
set -euo pipefail
unset POPLOG_FILE POPLOG_STATE_DIR POPLOG_PID POPLOG_PROGRESS_FILE POPLOG_STEP POPLOG_RUNNING POPLOG_JOB \
  POPSLOT_HELD POPSLOT_SLOT POPSLOT_FD POPSLOT_HOLDER_PID POPSLOT_AGENT POPSLOT_PRIORITY CARGO_BUILD_JOBS RUST_TEST_THREADS \
  GATE_STEPS
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$("$root/tools/selftest_scratch.sh" test_gate_p6) || exit 2; trap 'rm -rf "$tmp"' EXIT
fail=0
chk="$root/tools/gate_no_32bit.sh"
[[ -x "$chk" ]] || { echo "FAIL: tools/gate_no_32bit.sh missing"; fail=1; chk=false; }
for t in i686-unknown-linux-musl i586-unknown-linux-gnu thumbv7em-none-eabihf armv7-unknown-linux-gnueabihf wasm32-unknown-unknown riscv32imac-unknown-none-elf 'target_pointer_width = "32"'; do
  f="$tmp/dirty"; mkdir -p "$f/tools" "$f/scripts"; rm -f "$f"/tools/* "$f"/scripts/*
  printf 'cargo build --target x86_64-unknown-none\n' > "$f/tools/ok.sh"; printf '# %s\n' "$t" > "$f/scripts/bad.sh"
  rc=0; out=$("$chk" "$f" 2>&1) || rc=$?
  [[ $rc -eq 1 && "$out" == *scripts/bad.sh* ]] || { echo "FAIL[$t]: not flagged (rc=$rc)"; fail=1; }
done
c="$tmp/clean"; mkdir -p "$c/tools" "$c/scripts"
printf 'cargo build --target x86_64-unknown-none\naarch64-unknown-linux-gnu\n' > "$c/tools/a.sh"
printf 'refuse i686-unknown-linux-musl\n' > "$c/scripts/check_64bit_only.sh"
rc=0; "$chk" "$c" >/dev/null 2>&1 || rc=$?; [[ $rc -eq 0 ]] || { echo "FAIL[clean+allowlist]: rc=$rc"; fail=1; }
rc=0; out=$("$chk" "$root" 2>&1) || rc=$?; [[ $rc -eq 0 ]] || { echo "FAIL[real tree]: rc=$rc $out"; fail=1; }
for s in scripts/ci_nostd.sh scripts/ci.sh scripts/zenith_prepublish.sh; do
  grep -q -- '--target x86_64-unknown-none' "$root/$s" || { echo "FAIL[$s]: no x86_64-unknown-none build"; fail=1; }
done
mkdir -p "$tmp/bin"
printf '#!/usr/bin/env bash\necho "$*" >> "${0%%/bin/cargo}/cargo.args"\nexit 0\n' > "$tmp/bin/cargo"
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/bin/rustup"; chmod +x "$tmp/bin/cargo" "$tmp/bin/rustup"
export CARGO_TARGET_DIR="$tmp/target" POPSLOT_PREFIX="$tmp/slot" POPLOG_AGENT=gate-selftest STEP_TIMEOUT=60s
for s in build-nostd-none no-32bit-scripts; do
  : > "$tmp/cargo.args"; rc=0
  env GATE_STEPS=$s HOME="$tmp" PATH="$tmp/bin:$PATH" POPLOG_ROOT="$tmp/pop-$s" T="$tmp/gate-$s" \
    bash "$root/tools/full_gate.sh" "p6-$s" >/dev/null 2>&1 || rc=$?
  log=$(find "$tmp/pop-$s" -name '*.log' 2>/dev/null | sort | sed -n 1p || true)
  [[ -n "$log" ]] && grep -q "STEP END $s verdict=PASS" "$log" && [[ $rc -eq 0 ]] || { echo "FAIL[$s]: gate step did not PASS (rc=$rc)"; fail=1; }
  if [[ $s == build-nostd-none ]]; then
    grep -qx -- 'build --no-default-features --lib --target x86_64-unknown-none' "$tmp/cargo.args" \
      && grep -qx -- 'build -p zenith-float-num --no-default-features --lib --target x86_64-unknown-none' "$tmp/cargo.args" \
      || { echo "FAIL[$s]: cargo args: $(tr '\n' '|' < "$tmp/cargo.args")"; fail=1; }
  fi
done
[[ $fail -eq 0 ]] && echo "test_gate_p6: ok"
exit "$fail"
