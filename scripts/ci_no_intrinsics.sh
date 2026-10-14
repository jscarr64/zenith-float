#!/usr/bin/env bash
# Fail when Rust sources use platform intrinsics, inline asm, hardware float
# tokens, or `unsafe` outside scripts/no_intrinsics_allowlist.txt.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

fail=0

if rg -n --glob '*.rs' \
  -e 'std::arch' \
  -e 'core::arch' \
  -e '\basm!' \
  -e 'global_asm!' \
  -e '#\[target_feature' \
  -e 'is_x86_feature_detected' \
  -e '\bf32\b' \
  -e '\bf64\b' \
  .; then
  echo "error: platform intrinsics, inline asm, or hardware float tokens in Rust sources" >&2
  fail=1
fi

allow=()
while IFS= read -r line || [[ -n "$line" ]]; do
  case "$line" in
    ''|'#'*) continue ;;
  esac
  allow+=("${line%%[[:space:]]*}")
done < "$root/scripts/no_intrinsics_allowlist.txt"

while IFS= read -r hit; do
  [[ -z "$hit" ]] && continue
  file=${hit%%:*}
  file=${file#./}
  permitted=0
  for prefix in "${allow[@]}"; do
    if [[ "$file" == "$prefix" || "$file" == "$prefix"* ]]; then
      permitted=1
      break
    fi
  done
  if [[ "$permitted" -eq 0 ]]; then
    echo "error: unsafe outside the allowlist: $hit" >&2
    fail=1
  fi
done < <(rg -n --glob '*.rs' '\bunsafe\b' . || true)

if [[ "$fail" -ne 0 ]]; then
  exit 1
fi

echo "ci_no_intrinsics: ok"
