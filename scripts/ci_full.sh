#!/usr/bin/env bash
# Full maintainer gate (plan `ci_full.sh`). Same 12 checks as zenith_prepublish.sh.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
exec bash "$root/scripts/zenith_prepublish.sh"
