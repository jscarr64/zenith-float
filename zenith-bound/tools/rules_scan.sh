#!/usr/bin/env bash
# Scan ADDED diff lines for banned constructs in Rust code. See the header of rules_scan.py.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap rules_scan "$@"   # run log (tools/runlog.sh)
python3 "$(dirname "$0")/rules_scan.py" "$@"; exit $?
