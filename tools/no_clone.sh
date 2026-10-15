#!/usr/bin/env bash
# Zero-clone gate: exit 1 on any .clone()/.to_vec()/.cloned()/Clone derive or impl in non-test library code.
# Usage: tools/no_clone.sh [--list] [--files]
exec python3 "$(dirname "$0")/no_clone.py" "$@"
