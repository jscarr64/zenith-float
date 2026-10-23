#!/usr/bin/env bash
# Self-test: a full_gate.sh scan step must FAIL when its search tool is missing (exit 127)
# or errors (exit 2), must pass on the clean tree, and must fail on a planted violation.
# usage: tools/test_gate_tool_checks.sh   (runs only the scan steps, no cargo)
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

if [[ "$bad" -ne 0 ]]; then exit 1; fi
echo "test_gate_tool_checks: ok"
