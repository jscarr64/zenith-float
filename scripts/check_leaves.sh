#!/usr/bin/env bash
# expr! / cexpr! leaf match arms exist in the macros.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
exec python3 "$root/scripts/prepublish_lib.py" leaves
