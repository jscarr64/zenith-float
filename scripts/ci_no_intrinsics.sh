#!/usr/bin/env bash
# Fail when Rust sources use platform intrinsics, inline asm, hardware float
# tokens, or `unsafe` outside scripts/no_intrinsics_allowlist.txt.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if ! command -v git >/dev/null 2>&1; then
  echo "error: git not found" >&2
  exit 1
fi

fail=0

# git grep exits 0 on a match, 1 on no match, and 2 or more on an error.
status=0
git grep -nE \
  -e 'std::arch' \
  -e 'core::arch' \
  -e '\basm!' \
  -e 'global_asm!' \
  -e '#\[target_feature' \
  -e 'is_x86_feature_detected' \
  -e '\bf32\b' \
  -e '\bf64\b' \
  -- '*.rs' || status=$?
if [[ "$status" -ge 2 ]]; then
  echo "error: git grep failed with status $status" >&2
  exit 1
fi
if [[ "$status" -eq 0 ]]; then
  echo "error: platform intrinsics, inline asm, or hardware float tokens in Rust sources" >&2
  fail=1
fi

status=0
unsafe_hits="$(git grep -nE '\bunsafe\b' -- '*.rs')" || status=$?
if [[ "$status" -ge 2 ]]; then
  echo "error: git grep failed with status $status" >&2
  exit 1
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
done <<< "$unsafe_hits"

if [[ "$fail" -ne 0 ]]; then
  exit 1
fi

echo "ci_no_intrinsics: ok"
