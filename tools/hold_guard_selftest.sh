#!/usr/bin/env bash
# Selftest for tools/hold_guard.sh (Paul 2026-11-25 9:31 AM CT): fixture git repos with a held file and HOLD.md. An edit,
# an unpinned held file, a deletion or a taken-live file without a manifest line citing Jeff, and a ruling that does
# not name Jeff must FAIL; the pinned tree and changes recorded with a Jeff ruling must PASS. Then the real tree.
# usage: tools/hold_guard_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
guard="$root/tools/hold_guard.sh"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
zf_tmp="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/hold-guard-selftest"; mkdir -p "$zf_tmp" || exit 2  # job scratch
w=$(mktemp -d "$zf_tmp/run.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
held='// HOLD H99: held, not compiled.
// // fn f() -> u8 { 1 }'
sha() { sha256sum "$1" | cut -d' ' -f1; }
fixture() {  # a repo with src/h.rs held and HOLD.md, both pinned under a Jeff ruling
  n=$((n + 1)); r="$w/r$n"; mkdir -p "$r/src" "$r/tools"
  printf '%s\n' "$held" > "$r/src/h.rs"; printf '# HOLD index\n' > "$r/HOLD.md"
  printf '%s  src/h.rs  Jeff 2026-11-10 HOLD method\n%s  HOLD.md  Jeff 2026-11-10 HOLD method\n' \
    "$(sha "$r/src/h.rs")" "$(sha "$r/HOLD.md")" > "$r/tools/hold_manifest.txt"
}
commit() { git -C "$r" add -A && git -C "$r" -c user.name=t -c user.email=t@t commit -qm f; }
expect() {  # expect pass|fail [pattern]
  "$guard" --root "$r" > "$w/out$n" 2>&1; local rc=$?
  if [[ $1 == pass ]]; then
    [[ $rc -eq 0 ]] || bad "case $n want PASS, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-220)"
  else
    [[ $rc -eq 1 ]] && grep -q -- "$2" "$w/out$n" || bad "case $n want FAIL /$2/, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-220)"
  fi
}
git_init() { git -C "$r" init -q; }
fixture; git_init; commit; expect pass
fixture; git_init; printf '// // fn g() {}\n' >> "$r/src/h.rs"; commit; expect fail 'src/h.rs: CHANGED'
fixture; git_init; printf '# HOLD index\n| H99 | x |\n' > "$r/HOLD.md"; commit; expect fail 'HOLD.md: CHANGED'
fixture; git_init; printf '%s\n' "$held" > "$r/src/h2.rs"; commit; expect fail 'src/h2.rs: held file not in'
fixture; git_init; rm "$r/src/h.rs"; commit; expect fail 'src/h.rs: held .* but missing'
fixture; git_init; printf 'fn f() -> u8 { 1 }\n' > "$r/src/h.rs"; commit; expect fail 'HOLD header is gone'
fixture; git_init; sed -i 's/Jeff 2026-11-10 HOLD method/Paul call/' "$r/tools/hold_manifest.txt"; commit; expect fail 'does not name `Jeff'
# recorded with Jeff's ruling: an edit, a deletion, taking a file live
fixture; git_init; printf '// // fn g() {}\n' >> "$r/src/h.rs"
sed -i "s|^[0-9a-f]*  src/h.rs .*|$(sha "$r/src/h.rs")  src/h.rs  Jeff 2026-12-01 edit approved|" "$r/tools/hold_manifest.txt"; commit; expect pass
fixture; git_init; rm "$r/src/h.rs"; sed -i 's|^[0-9a-f]*  src/h.rs .*|deleted  src/h.rs  Jeff 2026-12-01 delete H99|' "$r/tools/hold_manifest.txt"; commit; expect pass
fixture; git_init; printf 'fn f() -> u8 { 1 }\n' > "$r/src/h.rs"; sed -i 's|^[0-9a-f]*  src/h.rs .*|live  src/h.rs  Jeff 2026-12-01 take H99 live|' "$r/tools/hold_manifest.txt"; commit; expect pass
fixture; git_init; sed -i 's|^[0-9a-f]*  src/h.rs .*|live  src/h.rs  Jeff 2026-12-01 take H99 live|' "$r/tools/hold_manifest.txt"; commit; expect fail 'still has a HOLD header'
r="$root"; n=$((n + 1)); expect pass
[[ $fail -eq 0 ]] && echo "hold_guard_selftest: ok ($((n - 1)) cases + tree)" || exit 1
