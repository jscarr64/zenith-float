#!/usr/bin/env bash
# Fail when the tracked tree mentions the private workspace name (case-insensitive).
# git grep exits 0 on a match, 1 on no match, 2+ on an error; an error (or a missing git) fails, it is not "no match".
set -euo pipefail
command -v git >/dev/null 2>&1 || { echo "error: git not found" >&2; exit 2; }
status=0
git grep -n -i 'accumath' -- . ':!scripts/ci.sh' ':!scripts/ci_no_name.sh' ':!tools/test_gate_tool_checks.sh' ':!CHANGELOG.md' || status=$?
if [[ "$status" -ge 2 ]]; then
  echo "error: git grep failed with status $status" >&2
  exit 2
fi
if [[ "$status" -eq 0 ]]; then
  echo "error: this public crate must not mention that name" >&2
  exit 1
fi
echo "ci_no_name: ok"
