#!/usr/bin/env bash
# Whole-tree check that no script keeps work files in /tmp. See the header of no_tmp_scan.py.
exec python3 -B "$(dirname "$0")/no_tmp_scan.py" "$@"
