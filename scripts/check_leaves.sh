#!/usr/bin/env bash
# LIBRARY.md expr! / cexpr! complete leaf lists must match the proc-macro match arms.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
exec python3 "$root/scripts/prepublish_lib.py" leaves
