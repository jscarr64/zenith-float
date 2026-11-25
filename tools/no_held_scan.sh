#!/usr/bin/env bash
# Whole-tree check that no held, xfail, skip or inverted-expectation list is left. See the header of no_held_scan.py.
exec python3 -B "$(dirname "$0")/no_held_scan.py" "$@"
