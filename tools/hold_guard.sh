#!/usr/bin/env bash
# HOLD files change only with a recorded Jeff ruling. See the header of hold_guard.py.
exec python3 -B "$(dirname "$0")/hold_guard.py" "$@"
