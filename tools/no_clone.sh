#!/usr/bin/env bash
# Zero-clone gate: exit 1 on any .clone()/.to_vec()/.cloned()/Clone derive or impl in library AND test code (see tools/NO_CLONE_RULES.md).
# Usage: tools/no_clone.sh [--list] [--files] [--helpers]
exec python3 "$(dirname "$0")/no_clone.py" "$@"
