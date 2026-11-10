#!/usr/bin/env bash
# Sourced by the tools: zb_scratch prints (and creates) the parent dir for temp files and dirs. Never /tmp: a reboot wipes it and lost
# a pinned worktree and a WIP commit (2026-10-23, Jeff's rule). ZB_SCRATCH overrides; default $HOME/.cache/zenith-bound/scratch.
zb_scratch() { local d=${ZB_SCRATCH:-$HOME/.cache/zenith-bound/scratch}; mkdir -p "$d" || return 2; printf '%s\n' "$d"; }
