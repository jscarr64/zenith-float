#!/usr/bin/env bash
# Self-test: arm-tests.yml must run the workspace tests with --release, matching the
# build-release step of the local gate. Fails closed if the workflow or step is missing.
# usage: tools/test_arm_tests_release.sh   (reads the workflow only, no cargo)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
wf=.github/workflows/arm-tests.yml
[ -f "$wf" ] || { echo "FAIL: $wf missing" >&2; exit 1; }
line="$(grep -E -- "run: cargo test .*--workspace" "$wf")" || { echo "FAIL: no workspace cargo test step in $wf" >&2; exit 1; }
for flag in --release --locked --workspace; do
  case " $line " in
    *" $flag "*) ;;
    *) echo "FAIL: workspace test step lacks $flag: $line" >&2; exit 1 ;;
  esac
done
echo "PASS: $line"
