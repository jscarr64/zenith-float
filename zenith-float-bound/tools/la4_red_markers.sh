#!/usr/bin/env bash
# Red-marker gate (Spence 2026-10-24): every `red: not implemented` test marker in tests/**/*.rs names its CONTRACT item
# (`CONTRACT <SLICE>.<n>`, e.g. CONTRACT LA4.1), and slice LA4 may not be marked done (a `Status: done` line in the
# docs/CONTRACT.md "### Slice LA4" section) while any marker naming an LA4 item remains. usage: tools/la4_red_markers.sh [ROOT]
# Every CONTRACT id a marker names must have a **<id>** entry in docs/CONTRACT.md. `red: CONTRACT 16.<SLICE>.<n>: ...` markers are checked too. Exit 0 ok, 1 violation, 2 setup error (missing tool or CONTRACT.md).
set -u
export LC_ALL=C
root=${1:-$(cd "$(dirname "$0")/.." && pwd)}
for c in grep awk; do command -v "$c" > /dev/null || { echo "la4_red_markers: MISSING TOOL $c"; exit 2; }; done
contract="$root/docs/CONTRACT.md"
[ -f "$contract" ] || { echo "la4_red_markers: no docs/CONTRACT.md under $root"; exit 2; }
markers=$(grep -rn --include='*.rs' -F 'red: not implemented' "$root/tests" 2> /dev/null)
bad=$(printf '%s\n' "$markers" | grep -vE 'CONTRACT [A-Z]+[0-9]*\.[0-9]+' | grep -v '^$')
if [ -n "$bad" ]; then echo "la4_red_markers: FAIL red marker without its CONTRACT <SLICE>.<n> item:"; echo "$bad"; exit 1; fi
# every named CONTRACT id must be a **<id>** entry in CONTRACT.md (a prose mention is not an entry)
missing=""
for id in $(printf '%s\n' "$markers" | grep -oE 'CONTRACT [A-Z]+[0-9]*\.[0-9]+' | awk '{ print $2 }' | sort -u); do
  grep -qF "**$id**" "$contract" || missing="$missing $id"
done
if [ -n "$missing" ]; then
  for id in $missing; do echo "la4_red_markers: FAIL CONTRACT $id has no entry in docs/CONTRACT.md (want a **$id** line):"; printf '%s\n' "$markers" | grep -F "CONTRACT $id:"; done; exit 1
fi
# Harper (2026-10-29): `red: CONTRACT <id>: ...` markers name a slice clause id 16.<SLICE>.<n> (e.g. 16.S2.1). The id must be well-formed
# and have a `- **<id> <title>.**` entry in docs/CONTRACT.md (the form the 16.x clauses use; a prose mention is not an entry).
cmarkers=$(grep -rn --include='*.rs' -F 'red: CONTRACT ' "$root/tests" 2> /dev/null)
cbad=$(printf '%s\n' "$cmarkers" | grep -vE 'red: CONTRACT 16\.[A-Z]+[0-9]*\.[0-9]+: ' | grep -v '^$')
if [ -n "$cbad" ]; then echo "la4_red_markers: FAIL malformed 'red: CONTRACT' marker id (want 'red: CONTRACT 16.<SLICE>.<n>: ...'):"; echo "$cbad"; exit 1; fi
cmissing=""
for m in $(printf '%s\n' "$cmarkers" | grep -oE 'red: CONTRACT 16\.[A-Z]+[0-9]*\.[0-9]+:' | awk '{ print $3 }' | sort -u); do
  id=${m%:}
  grep -qE "^- \*\*${id//./\\.}[ *]" "$contract" || cmissing="$cmissing $id"
done
if [ -n "$cmissing" ]; then
  for id in $cmissing; do echo "la4_red_markers: FAIL CONTRACT $id has no entry in docs/CONTRACT.md (want a '- **$id <title>.**' line):"; printf '%s\n' "$cmarkers" | grep -F "CONTRACT $id:"; done; exit 1
fi
la4=$(printf '%s\n' "$markers" | grep -E 'CONTRACT LA4\.[0-9]+')
n=$(printf '%s\n' "$la4" | grep -c .)
done_line=$(awk '/^### Slice LA4([^0-9]|$)/ { f = 1; next } /^##/ { f = 0 } f && tolower($0) ~ /status:?[*_ ]*done/' "$contract")
if [ -n "$done_line" ] && [ "$n" -gt 0 ]; then
  echo "la4_red_markers: FAIL LA4 is marked done ($done_line) but $n 'red: not implemented' LA4 test(s) remain:"; echo "$la4"; exit 1
fi
echo "la4_red_markers: ok ($n LA4 red marker(s); LA4 marked done: $([ -n "$done_line" ] && echo yes || echo no))"
