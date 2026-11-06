#!/usr/bin/env bash
# Self-test: the pinned qemu-aarch64 ZF_TIME_FACTOR must be at least the
# measured s9_erf slowdown plus a stated margin.
#
# Measurement (2026-11-05 America/Chicago, Pop, ZF_TIME_FACTOR=8 off the tight
# native budget): s9_erf --release, 3 runs each, median cargo "finished in":
#   native x86_64 11.52s (11.32, 11.52, 11.72);
#   qemu-aarch64-static 71.99s (72.94, 71.99, 71.95);
#   ratio 6.25.
# Margin 0.25 -> required floor 6.50. The pin is the value ALL-X01 / ZF-T02
# wiring should export as ZF_TIME_FACTOR under qemu-aarch64.
# Logs: ~/pop-logs/harper/2026-11-05/1236-zf-qemu-tf.log,
#       ~/pop-logs/harper/2026-11-05/1237-zf-qemu-tf.log
#
# Pin file: tools/qemu_aarch64_zf_time_factor.txt (one line, decimal in [1, 8]).
# usage: tools/test_qemu_aarch64_time_factor.sh
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
pin_file=tools/qemu_aarch64_zf_time_factor.txt
[ -f "$pin_file" ] || { echo "FAIL: $pin_file missing" >&2; exit 1; }
pin="$(tr -d " \t\r\n" < "$pin_file")"
case "$pin" in
  ""|*[!0-9.]*|*.*.*) echo "FAIL: pin $pin_file is not a plain decimal: $pin" >&2; exit 1 ;;
esac
# measured ratio 6.25 + margin 0.25 = 6.50, as milli for exact compare
floor_milli=6500
int="${pin%%.*}"
frac="${pin#*.}"
if [ "$frac" = "$pin" ]; then frac=""; fi
if [ "${#frac}" -gt 3 ]; then echo "FAIL: pin has more than 3 fraction digits: $pin" >&2; exit 1; fi
while [ "${#frac}" -lt 3 ]; do frac="${frac}0"; done
pin_milli=$((10#$int * 1000 + 10#${frac:-0}))
if [ "$pin_milli" -lt "$floor_milli" ]; then
  echo "FAIL: pinned qemu-aarch64 ZF_TIME_FACTOR=$pin (${pin_milli}/1000) < floor 6.50 (measured s9_erf ratio 6.25 + margin 0.25)" >&2
  exit 1
fi
echo "PASS: pinned qemu-aarch64 ZF_TIME_FACTOR=$pin (>= 6.50 = measured 6.25 + margin 0.25)"
