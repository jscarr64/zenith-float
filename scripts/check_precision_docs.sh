#!/usr/bin/env bash
# Every special in LIBRARY.md §15 / CAPABILITIES.md §12 has rustdoc `# Precision`.
# (Plan §19.4 still says "§12 of LIBRARY.md"; specials moved to LIBRARY §15.)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
exec python3 "$root/scripts/prepublish_lib.py" precision
