#!/usr/bin/env bash
# Self-test: a full_gate.sh scan step must FAIL when its search tool is missing (exit 127)
# or errors (exit 2), must pass on the clean tree, and must fail on a planted violation.
# usage: tools/test_gate_tool_checks.sh   (runs only the scan steps and the prepublish compare step with stub cargo; no real cargo)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
t="$(mktemp -d)"
trap 'rm -rf "$t"' EXIT
bad=0
fail() { echo "FAIL: $*" >&2; bad=1; }

# the command string of `run NAME "CMD"` in full_gate.sh
step_cmd() {
  local line
  line="$(grep -E "^run $1 \"" tools/full_gate.sh)" || { echo "no step $1 in full_gate.sh" >&2; return 2; }
  line="${line#run "$1" \"}"
  printf '%s\n' "${line%\"}"
}

# stub <dir> <tool> <rc>: a fake tool that exits rc without output
stub() { mkdir -p "$1"; printf '#!/bin/sh\nexit %s\n' "$3" > "$1/$2"; chmod +x "$1/$2"; }

# step|search tool it uses|planted file|planted line
CASES=(
  "unsafe-scan|grep|zenith-float-num/src/lib.rs|pub fn planted() { unsafe { core::hint::unreachable_unchecked() } }"
  "no-intrinsics|git|zenith-float-num/src/lib.rs|pub fn planted(x: f64) -> f64 { x }"
)
for c in "${CASES[@]}"; do
  IFS='|' read -r step tool pfile pline <<< "$c"
  cmd="$(step_cmd "$step")"
  for rc in 127 2; do
    stub "$t/bin-$tool-$rc" "$tool" "$rc"
    if (cd "$root" && PATH="$t/bin-$tool-$rc:$PATH" bash -c "$cmd") >/dev/null 2>&1; then
      fail "$step passes when $tool exits $rc"
    fi
  done
  (cd "$root" && bash -c "$cmd") >/dev/null 2>&1 || fail "$step fails on the clean tree"
  copy="$t/copy-$step"
  mkdir -p "$copy" && git archive HEAD | tar -x -C "$copy"
  cp -r tools scripts "$copy/" 2>/dev/null || true
  (cd "$copy" && git init -q && git add -A >/dev/null)
  printf '%s\n' "$pline" >> "$copy/$pfile"
  (cd "$copy" && git add -A >/dev/null)
  if (cd "$copy" && bash -c "$cmd") >/dev/null 2>&1; then
    fail "$step passes with a planted violation in $pfile"
  fi
done

# scripts/ci.sh public-name check: once it lives in scripts/ci_no_name.sh run that, else the inline `if rg` block.
name_check() {
  if [[ -f "$root/scripts/ci_no_name.sh" ]]; then bash "$root/scripts/ci_no_name.sh"
  else bash -c "set -euo pipefail; $(sed -n "/^  if rg .*accumath/,/^  fi$/p" "$root/scripts/ci.sh")"; fi
}
nc="$t/namecheck"; mkdir -p "$nc"
(cd "$nc" && git init -q && echo 'fn a() {}' > ok.rs && git add ok.rs)
(cd "$nc" && name_check) >/dev/null 2>&1 || fail "name check fails on a tree without the name"
echo '// accumath' > "$nc/bad.rs"; (cd "$nc" && git add bad.rs)
(cd "$nc" && name_check) >/dev/null 2>&1 && fail "name check passes with the name planted"
for tool in rg git; do
  for rc in 127 2; do
    stub "$t/bin-n-$tool-$rc" "$tool" "$rc"
    (cd "$nc" && PATH="$t/bin-n-$tool-$rc:$PATH" name_check) >/dev/null 2>&1 && fail "name check passes with the name planted when $tool exits $rc"
  done
done

# untracked, non-ignored files count (rg semantics); .gitignore still applies
nu="$t/name-untracked"; mkdir -p "$nu/target"
(cd "$nu" && git init -q && echo 'fn a() {}' > ok.rs && printf 'target/\n' > .gitignore && git add ok.rs .gitignore)
echo '// accumath' > "$nu/target/ignored.rs"
(cd "$nu" && name_check) >/dev/null 2>&1 || fail "name check flags a gitignored file"
echo '// accumath' > "$nu/new.rs"
(cd "$nu" && name_check) >/dev/null 2>&1 && fail "name check passes with the name in an untracked, non-ignored file"

# A gate step never rewrites a tracked file (Jeff 2026-11-10; gate3 1951-fold-pr1-gate3 was NOT A GATE because the prepublish
# compare-bench step rewrote doc/compare-results.tsv). The prepublish compare step (its lines, cut from scripts/zenith_prepublish.sh)
# runs in a shared clone with stub cargo/rustup; `git status --porcelain` must not change, and the output must land outside the
# clone. A planted copy that calls compare-bench.sh without COMPARE_OUT must be caught.
compare_step() {  # FILE -> the prepublish compare step's lines
  awk '/^step "scripts\/compare-bench.sh --quick"/{p=1; next} p&&/^ok "compare-bench --quick"/{exit} p{print}' "$1"
}
cs="$t/compare"; mkdir -p "$cs/bin" "$cs/tgt"
printf '#!/usr/bin/env bash\nif [[ "$1" == run && "${2:-}" == 1.99 ]]; then echo "rustc 1.99.0 (stub)"; fi\nexit 0\n' > "$cs/bin/rustup"
printf '#!/usr/bin/env bash\necho "library\tprecision_bits\ttask\tbatch_size\tbest_time_us"\necho "zenith\t132\tadd\t1\t1"\n' > "$cs/bin/cargo"
chmod +x "$cs/bin/rustup" "$cs/bin/cargo"
compare_clean() {  # PREPUBLISH-FILE -> 0 when the compare step leaves the clone's tracked files alone
  local c="$cs/clone-$1" before after body
  rm -rf "$c"; git clone -q --shared "$root" "$c" || return 2
  cp "$root/scripts/compare-bench.sh" "$c/scripts/compare-bench.sh"
  body=$(compare_step "$2"); [ -n "$body" ] || { echo "  no compare step in $2"; return 1; }
  before=$(git -C "$c" status --porcelain)
  (cd "$c" && root=$c CARGO_TARGET_DIR="$cs/tgt/$1" PATH="$cs/bin:$PATH" bash -euo pipefail -c "$body") > "$cs/$1.out" 2>&1 \
    || { echo "  compare step failed: $(tail -n 3 "$cs/$1.out")"; return 1; }
  after=$(git -C "$c" status --porcelain)
  [ "$before" = "$after" ] || { echo "  compare step rewrote tracked files: $after"; return 1; }
  grep -q "^Wrote $cs/tgt/$1/" "$cs/$1.out" || { echo "  compare output not under CARGO_TARGET_DIR: $(grep '^Wrote' "$cs/$1.out")"; return 1; }
}
compare_clean real "$root/scripts/zenith_prepublish.sh" || fail "prepublish compare step rewrites a tracked file (or writes inside the worktree)"
sed 's#^COMPARE_OUT=.* bash "\$root/scripts/compare-bench.sh" --quick$#bash "$root/scripts/compare-bench.sh" --quick#' "$root/scripts/zenith_prepublish.sh" > "$cs/planted.sh"
grep -qx 'bash "$root/scripts/compare-bench.sh" --quick' "$cs/planted.sh" || fail "compare plant: no bare compare-bench line"
compare_clean planted "$cs/planted.sh" > /dev/null && fail "planted bare compare-bench (default doc/compare-results.tsv) was not caught"

if [[ "$bad" -ne 0 ]]; then exit 1; fi
echo "test_gate_tool_checks: ok"
