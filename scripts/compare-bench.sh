#!/usr/bin/env bash
# Cross-library comparison (bigfloat-bench compatible workloads, TSV output).
#
# Quick release check (132-bit, core ops):
#   ./scripts/compare-bench.sh --quick
#
# Full matrix (132 / 1000 / 10000 bits, all tasks):
#   ./scripts/compare-bench.sh --full
#
# Output: doc/compare-results.tsv (tab-separated, no JSON)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

out="${COMPARE_OUT:-doc/compare-results.tsv}"
libs=(--lib zenith --lib astro)
features=(--features zenith-float-compare/astro)
precisions=(--precision 132)
tasks=(
  --task add --task sub --task mul --task div
  --task sqrt --task ln --task exp --task sin --task cos
)

if [[ "${1:-}" == "--full" ]]; then
  precisions=(--precision 132 --precision 1000 --precision 10000)
  tasks=(
    --task add --task sub --task mul --task div
    --task sqrt --task cbrt --task ln --task exp --task pow
    --task sin --task asin --task cos --task acos --task tan --task atan
    --task sinh --task asinh --task cosh --task acosh --task tanh --task atanh
  )
fi

if [[ "${1:-}" == "--dashu" ]]; then
  libs+=(--lib dashu)
  features=(--features "zenith-float-compare/astro,zenith-float-compare/dashu")
fi

{
  echo -e "library\tprecision_bits\ttask\tbatch_size\tbest_time_us"
  cargo run -p zenith-float-compare --release "${features[@]}" -- \
    "${libs[@]}" "${precisions[@]}" "${tasks[@]}" -n 5 \
    | tail -n +2
} > "$out"

echo "Wrote $out ($(wc -l < "$out") lines)"
