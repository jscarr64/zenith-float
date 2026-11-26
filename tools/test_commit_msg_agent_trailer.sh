#!/usr/bin/env bash
# Red test (Spence B3 on 82e4afb: 15 commits after the 5bfe97b pre-commit hook still lack an `Agent:` trailer, because no
# hook reads the message). .githooks/commit-msg must refuse a commit whose message has no `Agent: <name>` trailer, and
# accept one that has it. Fixture repo with core.hooksPath on a copy of the repo's .githooks/commit-msg only (the
# pre-commit hook is tools/test_rules_enforcement.sh part 4).
# usage: tools/test_commit_msg_agent_trailer.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
fail=0
bad() { echo "FAIL: $*"; fail=1; }
w=$(mktemp -d "${TMPDIR:-/tmp}/commit-msg-trailer.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
repo="$w/repo"; mkdir -p "$repo/.githooks"
[[ -f "$root/.githooks/commit-msg" ]] && cp "$root/.githooks/commit-msg" "$repo/.githooks/" && chmod +x "$repo/.githooks/commit-msg"
g() { git -C "$repo" -c user.name=t -c user.email=t@t "$@"; }
g init -q && g config core.hooksPath .githooks
n=0
# try WANT MESSAGE: WANT is ok or refused
try() {
  local want=$1 msg=$2 before after
  n=$((n + 1)); echo "$n" > "$repo/f$n"; g add -A
  before=$(g rev-list --count --all 2>/dev/null || echo 0)
  printf '%s' "$msg" > "$w/msg"
  g commit -q -F "$w/msg" > "$w/out" 2>&1
  after=$(g rev-list --count --all 2>/dev/null || echo 0)
  if [[ $want == ok ]]; then
    [[ $after -eq $((before + 1)) ]] || bad "refused a good message: $(head -1 "$w/msg") ($(tr '\n' ' ' < "$w/out"))"
  else
    [[ $after -eq $before ]] || bad "accepted a message without an Agent trailer: $(tr '\n' '|' < "$w/msg")"
    g reset -q
  fi
}
try ok $'base\n\nAgent: paul\n'
try refused $'no trailer at all\n'
try refused $'body only\n\nSome text that says Agent: paul inside a sentence.\n'
try refused $'empty name\n\nAgent:\n'
try refused $'trailer not at the end\n\nAgent: paul\n\nmore body after it\n'
try refused $'lower-case key\n\nagent: paul\n'
try ok $'with other trailers\n\nbody\n\nCo-authored-by: x <x@x>\nAgent: rusty\n'
try ok $'comment lines are stripped\n\nAgent: harper\n# Please enter the commit message\n'
[[ $fail -eq 0 ]] && { echo "test_commit_msg_agent_trailer: ok"; exit 0; }
echo "test_commit_msg_agent_trailer: FAIL"; exit 1
