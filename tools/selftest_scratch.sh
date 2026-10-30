#!/usr/bin/env bash
# usage: tools/selftest_scratch.sh NAME  -> prints a fresh dir $SELFTEST_SCRATCH/NAME.XXXXXX (the caller removes it).
# Gate selftests take their scratch only from here: SELFTEST_SCRATCH unset, or under /tmp, is an error (exit 2).
# full_gate.sh sets SELFTEST_SCRATCH=$T/scratch when the caller has not.
set -euo pipefail
name=${1:?usage: tools/selftest_scratch.sh NAME}
s=${SELFTEST_SCRATCH:-}
[[ -n $s ]] || { echo "selftest_scratch: SELFTEST_SCRATCH is unset (set it to a dir outside /tmp)" >&2; exit 2; }
real=$(readlink -m -- "$s")
case $real in /tmp | /tmp/*) echo "selftest_scratch: SELFTEST_SCRATCH=$s is under /tmp" >&2; exit 2;; esac
mkdir -p -- "$real"
mktemp -d "$real/$name.XXXXXX"
