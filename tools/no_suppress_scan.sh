#!/usr/bin/env bash
# Whole-tree hard fail on any lint suppression or unsafe code. See the header of no_suppress_scan.py.
exec python3 -B "$(dirname "$0")/no_suppress_scan.py" "$@"
