#!/usr/bin/env bash
# usage: tools/toolchain_check.sh [rust-toolchain.toml]   (gate step `toolchain`)
# Prints which rustc/cargo the gate runs and their versions, and fails (exit 1, "TOOLCHAIN MISMATCH") unless both report the version
# pinned by the toolchain file's `channel` ("1.99" matches 1.99 and 1.99.x, not 1.990). A missing file, a channel that is not a pinned
# version (e.g. "stable"), or a missing rustc/cargo also fails. Why: on Pop /usr/bin/rustc is 1.95 and a PATH without ~/.cargo/bin first
# silently builds with it.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap toolchain_check "$@"   # run log (tools/runlog.sh)
set -u
f=${1:-rust-toolchain.toml}
[ -r "$f" ] || { echo "toolchain: cannot read $f"; exit 2; }
want=$(sed -n 's/^[[:space:]]*channel[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "$f")
[[ "$want" =~ ^[0-9]+\.[0-9]+(\.[0-9]+)?$ ]] || { echo "toolchain: channel '$want' in $f is not a pinned version"; exit 2; }
rc=0
for tool in rustc cargo; do
  path=$(command -v "$tool") || { echo "MISSING TOOL: $tool (step cannot run)"; exit 127; }
  have=$("$tool" --version 2>&1) || { echo "toolchain: $tool --version failed: $have"; exit 2; }
  echo "$tool: $path -> $have (pinned channel $want)"
  ver=${have#"$tool "}; ver=${ver%% *}
  if [ "$ver" != "$want" ] && [[ "$ver" != "$want".* ]]; then echo "TOOLCHAIN MISMATCH: $tool $ver is not channel $want"; rc=1; fi
done
exit "$rc"
