#!/usr/bin/env bash
# Every special in CAPABILITIES.md has rustdoc `# Precision`.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
exec python3 "$root/scripts/prepublish_lib.py" precision
