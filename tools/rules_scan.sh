#!/usr/bin/env bash
# Scan ADDED diff lines for banned constructs in library code. See the header of rules_scan.py.
exec python3 -B "$(dirname "$0")/rules_scan.py" "$@"
