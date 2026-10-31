#!/usr/bin/env bash
# P6 (64-bit targets only): FAIL (exit 1) when a gate script under ROOT/tools or ROOT/scripts names a 32-bit
# target (triple or target_pointer_width = "32"). Allowed: scripts/check_64bit_only.sh (the refuse-32bit check
# must name the 32-bit targets it expects to refuse), this checker and its selftest.
set -euo pipefail
root=${1:-$(cd "$(dirname "$0")/.." && pwd)}
pat='(i[3-6]86-|thumbv[0-9]|armv[4-7]|arm-unknown|armebv|wasm32|riscv32|mipsel?-|powerpc-|sparc-|m68k|avr-|msp430|x86_32|pointer_width *= *"32")'
hits=$(cd "$root" && grep -rnE "$pat" tools scripts 2>/dev/null \
  | grep -vE '^(scripts/check_64bit_only\.sh|tools/gate_no_32bit\.sh|tools/test_gate_p6\.sh):' || true)
if [[ -n "$hits" ]]; then printf 'no-32bit-scripts: FAIL (32-bit target in gate scripts, P6)\n%s\n' "$hits" >&2; exit 1; fi
echo "no-32bit-scripts: ok"
