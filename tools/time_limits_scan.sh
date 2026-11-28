#!/usr/bin/env bash
# usage: tools/time_limits_scan.sh  (see tools/time_limits_scan.py)
exec python3 "$(cd "$(dirname "$0")" && pwd)/time_limits_scan.py" "$@"
