#!/usr/bin/env bash
# Regenerate the published doc/compare-results.tsv on purpose (the gate and prepublish never write it).
# usage: tools/regen_compare_results.sh [--quick|--full]   (default --quick)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
COMPARE_OUT="$root/doc/compare-results.tsv" exec bash "$root/scripts/compare-bench.sh" "${@:---quick}"
