#!/usr/bin/env bash
# Criterion benchmarks for zenith-float-num.
#
# Quick smoke (few samples, for local/CI):  ./scripts/bench.sh
# Full run with HTML report:               ./scripts/bench.sh --full
# Single suite:                            cargo bench -p zenith-float-num --bench arithmetic
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

benches=(--bench arithmetic --bench transcendentals --bench composite --bench specials --bench linalg)

if [[ "${1:-}" == "--full" ]]; then
  cargo bench -p zenith-float-num "${benches[@]}"
else
  cargo bench -p zenith-float-num "${benches[@]}" -- --quick
fi
