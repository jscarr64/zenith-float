#!/usr/bin/env bash
# Whole-tree check that no Rust code uses debug_assert!. See the header of no_debug_assert_scan.py.
exec python3 -B "$(dirname "$0")/no_debug_assert_scan.py" "$@"
