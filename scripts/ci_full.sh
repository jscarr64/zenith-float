#!/usr/bin/env bash
# Full maintainer gate (plan `ci_full.sh`). Same 12 checks as zenith_prepublish.sh.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"
exec bash "$root/scripts/zenith_prepublish.sh"
