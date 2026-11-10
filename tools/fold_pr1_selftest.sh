#!/usr/bin/env bash
# Fold PR 1 selftest: the bound layer joins zenith-float as the workspace crate zenith-float-bound
# (RELEASE-2.0-CHECKLIST 4.1, 4.3, 4.5, 4.6; Rusty fold plan 2026-11-10 section 2.4, rows A..Y plus NOSTR, NOCYCLE, HOLD and HOLDIDX).
# usage: tools/fold_pr1_selftest.sh [BASE_REF]      (BASE_REF defaults to origin/main)
# Every row runs its real check on the COMMITTED HEAD tree (git archive into a scratch dir, never the working tree) and
# must PASS; then the row's planted defect is applied to a second copy and the same check must FAIL on it (DETECTED).
# A check that cannot go red proves nothing. Exit 0 only if every row is PASS + DETECTED, else 1; 2 on setup errors.
# Scratch: $FOLD_PR1_SCRATCH (default $HOME/.cache/zenith-float/fold-pr1). Cargo runs --offline; builds use
# $CARGO_TARGET_DIR (default under the scratch dir). FOLD_PR1_NO_BUILD=1 skips the three cargo build/test steps.
set -u
export LC_ALL=C
repo=$(git rev-parse --show-toplevel) || exit 2
cd "$repo" || exit 2
BASE=${1:-origin/main}
git rev-parse -q --verify "$BASE^{commit}" > /dev/null || { echo "fold-pr1: FAIL base $BASE is not a commit"; exit 2; }
ZB_TIP=fcaf135ee1787afecb03d90b3a87b22a6d5a04b5
ZB_COMMITS=612
PFX=zenith-float-bound
# The old crate name, as a regex and as text, built so that this file never holds it literally.
OLD_PAT='zenith[-_]bound'
OLD_NAME="zenith""-bound"
REMOVED="rust-toolchain.toml .github/CODEOWNERS"
HEAD_SHA=$(git rev-parse HEAD) || exit 2
NO_BUILD=${FOLD_PR1_NO_BUILD:-0}
scratch=${FOLD_PR1_SCRATCH:-$HOME/.cache/zenith-float/fold-pr1}
mkdir -p "$scratch" || exit 2
work=$(mktemp -d "$scratch/selftest.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
real="$work/real"; plant="$work/plant"
mkdir -p "$real" "$plant" || exit 2
git archive "$HEAD_SHA" | tar -x -C "$real" || exit 2
git archive "$HEAD_SHA" | tar -x -C "$plant" || exit 2
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-$scratch/target}
echo "fold-pr1: head=$HEAD_SHA base=$(git rev-parse "$BASE") zb_tip=$ZB_TIP target=$CARGO_TARGET_DIR no_build=$NO_BUILD"

# Restore one planted file from the real copy (or remove it when the real copy has none).
restore() { if [ -e "$real/$1" ]; then cp "$real/$1" "$plant/$1"; else rm -f "$plant/$1"; fi; }

# --- git rows (refs) ---------------------------------------------------------------------------------------------
check_A() { [ "$(git cat-file -t "$1" 2> /dev/null)" = commit ]; }
find_import() {
  local c
  for c in $(git rev-list --first-parent --merges "$BASE..$1" 2> /dev/null); do
    if [ "$(git rev-parse -q --verify "$c^2" 2> /dev/null)" = "$ZB_TIP" ]; then echo "$c"; return 0; fi
  done
  return 1
}
check_Y() {
  check_A "$ZB_TIP" || return 1
  git merge-base --is-ancestor "$ZB_TIP" "$1" || return 1
  [ "$(git rev-list --count "$ZB_TIP")" = "$ZB_COMMITS" ] || return 1
  find_import "$1" > /dev/null
}
# B (adapted to the new prefix): the import merge holds the zb tip tree byte-for-byte under $PFX/, and at the checked
# ref every zb tip path is still under $PFX/ except exactly the removed files ($2), which are gone. Later commits may add
# paths (topic-tip merges, the crate rustfmt.toml); a held file (HOLD.md) keeps its path, so it still counts.
check_B() {
  local m want have missing f
  m=$(find_import "$1") || return 1
  [ "$(git rev-parse "$m:$PFX" 2> /dev/null)" = "$(git rev-parse "$ZB_TIP^{tree}")" ] || return 1
  want=$(git ls-tree -r --name-only "$ZB_TIP" | grep -vxF -f <(printf '%s\n' $2)) || return 1
  have=$(git ls-tree -r --name-only "$1:$PFX" 2> /dev/null) || return 1
  missing=$(comm -23 <(printf '%s\n' "$want" | sort) <(printf '%s\n' "$have" | sort))
  [ -z "$missing" ] || { printf '  B missing: %s\n' "$missing" | head -n 20; return 1; }
  for f in $2; do
    if printf '%s\n' "$have" | grep -qxF "$f"; then echo "  B: removed file $f is back"; return 1; fi
  done
}
check_P() {
  local c last commits
  commits=$(git rev-list "$BASE..$1" "^$ZB_TIP" 2> /dev/null) || return 1
  [ -n "$commits" ] || return 1
  for c in $commits; do
    last=$(git log -1 --format=%B "$c" | sed '/^[[:space:]]*$/d' | tail -n 1)
    [ "$last" = "Agent: rusty" ] || { echo "  P: $c has no final 'Agent: rusty' line"; return 1; }
  done
}
check_O() {
  local f found=0
  for f in $(git diff --name-only --diff-filter=A "$BASE" "$1" -- changelog.d); do
    git show "$1:$f" | grep -qF "$PFX" && found=1
  done
  [ "$found" = 1 ]
}

# --- tree rows (directory of a committed tree) -------------------------------------------------------------------
check_D() {
  local m="$1/$PFX/Cargo.toml" line
  [ -f "$m" ] || return 1
  grep -qxF "name = \"$PFX\"" "$m" || return 1
  [ "$(grep -cE '^zenith-float-num[[:space:]]*=' "$m")" = 1 ] || return 1
  line=$(grep -E '^zenith-float-num[[:space:]]*=' "$m")
  printf '%s\n' "$line" | grep -qF 'path = "../zenith-float-num"' || return 1
  ! printf '%s\n' "$line" | grep -qE '(git|rev|branch|tag|version)[[:space:]]*='
}
check_T() {
  local hits
  hits=$(find "$1" -name Cargo.toml -not -path '*/target/*' -print0 | xargs -0 grep -nE '627fedd|git[[:space:]]*=[[:space:]]*"https://github.com/jscarr64/zenith')
  [ -z "$hits" ] || { printf '  T: %s\n' "$hits"; return 1; }
}
check_C() {
  local out
  out=$(cd "$1" && cargo tree --offline -p "$PFX" -e normal 2>&1) || { echo "  C: cargo tree failed: $(printf '%s' "$out" | tail -n 1)"; return 1; }
  printf '%s\n' "$out" | grep -q 'zenith-float-num v' || return 1
  ! printf '%s\n' "$out" | grep -qE '(^|[[:space:]])zenith-float v'
}
check_NOCYCLE() {
  local out
  out=$(cd "$1" && cargo tree --offline -p zenith-float-num -e normal,build,dev 2>&1) || return 1
  printf '%s\n' "$out" | grep -qF "$PFX v" && return 1
  grep -rqF --include='*.rs' 'zenith_float_bound' "$1/zenith-float-num" && return 1
  ! grep -qF "$PFX" "$1/zenith-float-num/Cargo.toml"
}
check_E() {
  local out
  out=$(cd "$1" && cargo metadata --offline --no-deps --format-version 1 2> /dev/null) || return 1
  printf '%s\n' "$out" | grep -qE "\"workspace_members\":\[[^]]*#$PFX@"
}
check_F() {
  local out
  out=$(cd "$1" && cargo tree --offline -p zenith-float --no-default-features -e normal 2>&1) || return 1
  printf '%s\n' "$out" | grep -qF "$PFX v" && return 1
  [ "$2" = build ] || return 0
  (cd "$1" && cargo build --offline --no-default-features) > "$work/F.build.log" 2>&1 || { tail -n 5 "$work/F.build.log"; return 1; }
}
check_G() {
  local out
  out=$(cd "$1" && cargo tree --offline -p zenith-float -e normal 2>&1) || return 1
  printf '%s\n' "$out" | grep -qF "$PFX v" || return 1
  grep -qE '^default[[:space:]]*=.*"bound"' "$1/Cargo.toml" || return 1
  grep -qE "^bound[[:space:]]*=.*\"dep:$PFX\"" "$1/Cargo.toml" || return 1
  [ "$2" = build ] || return 0
  (cd "$1" && cargo build --offline -p zenith-float) > "$work/G.build.log" 2>&1 || { tail -n 5 "$work/G.build.log"; return 1; }
}
check_L() {
  local lib="$1/src/lib.rs" out
  grep -qxF 'pub use zenith_float_num::*;' "$lib" || return 1
  grep -qxF 'pub use zenith_float_bound as bound;' "$lib" || return 1
  grep -B 3 -xF 'pub use zenith_float_bound as bound;' "$lib" | grep -qxF '#[cfg(feature = "bound")]' || return 1
  [ "$(grep -vE '^[[:space:]]*//' "$lib" | grep -cF 'zenith_float_bound')" = 1 ] || return 1
  [ "$2" = build ] || return 0
  out=$(cd "$1" && cargo test --offline -p zenith-float --test fold_pr1_facade 2>&1) || { printf '%s\n' "$out" | tail -n 5; return 1; }
  printf '%s\n' "$out" | grep -qF 'test fold_pr1_facade_ball_still_num_ball ... ok' || return 1
  printf '%s\n' "$out" | grep -qF 'test fold_pr1_bound_ball_is_bound_crate_ball ... ok'
}
check_Q() {
  [ ! -e "$1/$PFX/rust-toolchain.toml" ] && [ ! -e "$1/$PFX/rust-toolchain" ] || return 1
  [ "$(find "$1" -name 'rust-toolchain*' -not -path '*/target/*' | wc -l)" = 1 ] || return 1
  grep -qxF 'channel = "1.99"' "$1/rust-toolchain.toml"
}
check_R() {
  local co="$1/.github/CODEOWNERS" p
  [ -f "$co" ] || return 1
  for p in "/$PFX/" "/$PFX/Cargo.toml" "/$PFX/tools/full_gate.sh" "/$PFX/tools/digests.sh" "/$PFX/tools/digest_expected.txt" \
    "/$PFX/tools/digest_targets.txt" "/$PFX/tools/tightness_report.py" "/$PFX/tools/tightness_baseline.txt"; do
    awk -v p="$p" '$1 == p && $2 == "@jscarr64" { f = 1 } END { exit !f }' "$co" || { echo "  R: missing $p"; return 1; }
  done
  [ "$(find "$1" -name CODEOWNERS -not -path '*/target/*' | wc -l)" = 1 ]
}
lockfacts() {
  awk '
    function flush(  i) { if (name != "") { print "P " name " " ver " " src; for (i = 1; i <= nd; i++) print "D " name " " deps[i] } }
    /^\[\[package\]\]/ { flush(); name = ""; ver = ""; src = "-"; indeps = 0; nd = 0; next }
    /^name = /         { v = $3; gsub(/"/, "", v); name = v; next }
    /^version = /      { v = $3; gsub(/"/, "", v); ver = v; next }
    /^source = /       { v = $3; gsub(/"/, "", v); src = v; next }
    /^dependencies = \[$/ { indeps = 1; next }
    indeps && /^\]/    { indeps = 0; next }
    indeps             { d = $0; gsub(/[",]/, "", d); gsub(/^[[:space:]]+|[[:space:]]+$/, "", d); deps[++nd] = d; next }
    END { flush() }' | sort
}
check_S() {
  local added removed want
  git show "$BASE:Cargo.lock" | lockfacts > "$work/lock.base" || return 1
  lockfacts < "$1/Cargo.lock" > "$work/lock.head" || return 1
  added=$(comm -13 "$work/lock.base" "$work/lock.head")
  removed=$(comm -23 "$work/lock.base" "$work/lock.head")
  want=$(printf '%s\n' "D zenith-float $PFX" "D $PFX zenith-float-num" "P $PFX 2.0.0 -" | sort)
  [ -z "$removed" ] || { printf '  S removed: %s\n' "$removed"; return 1; }
  [ "$added" = "$want" ] || { printf '  S added: %s\n' "$added"; return 1; }
}
check_NOSTR() {
  local hits
  hits=$(grep -rIniE --exclude-dir=target "$OLD_PAT" "$1" | sed "s|^$1/||")
  [ -z "$hits" ] || { printf '  NOSTR: %s\n' "$hits" | head -n 20; return 1; }
}

# HOLD (Jeff 2026-11-10 15:36): held code lives commented out under a `// HOLD <id>:` first line, indexed in the one
# repo-root HOLD.md (table rows `| <id> | `path`, ... | ...`). Every held file's id has an index row naming that file,
# every index path exists and starts with its row's header, and no live (uncommented) `mod` line or Cargo target path
# names a held file.
check_HOLD() {
  local d=$1 idx ids id paths f hid stem dir hits
  idx="$d/HOLD.md"
  [ -f "$idx" ] || { echo "  HOLD: no HOLD.md at the repo root"; return 1; }
  ids=$(grep -oE '^\| H[0-9]+ \|' "$idx" | grep -oE 'H[0-9]+')
  [ -n "$ids" ] || { echo "  HOLD: HOLD.md has no index rows"; return 1; }
  [ "$(printf '%s\n' "$ids" | sort | uniq -d)" = "" ] || { echo "  HOLD: duplicate index ids"; return 1; }
  for id in $ids; do
    paths=$(grep -E "^\| $id \|" "$idx" | awk -F'|' '{print $3}' | grep -oE '`[^`]+`' | tr -d '`')
    [ -n "$paths" ] || { echo "  HOLD: $id names no path"; return 1; }
    for f in $paths; do
      [ -f "$d/$f" ] || { echo "  HOLD: $id path $f missing"; return 1; }
      head -n 1 "$d/$f" | grep -qE "^// HOLD $id:" || { echo "  HOLD: $f has no '// HOLD $id:' first line"; return 1; }
    done
  done
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    hid=$(head -n 1 "$d/$f" | grep -oE '^// HOLD H[0-9]+' | awk '{print $3}')
    grep -E "^\| $hid \|" "$idx" | awk -F'|' '{print $3}' | grep -qF "\`$f\`" || { echo "  HOLD: $f ($hid) has no index row"; return 1; }
    case $f in
      *.rs)
        stem=$(basename "$f" .rs); dir=$(dirname "$f")
        hits=$(grep -rnE --include='*.rs' "^[[:space:]]*(pub(\([a-z]+\))?[[:space:]]+)?mod[[:space:]]+$stem[[:space:]]*;" "$d/$dir" "$d/$(dirname "$dir")" 2> /dev/null | sed "s|^$d/||")
        [ -z "$hits" ] || { echo "  HOLD: live mod line for held $f: $hits"; return 1; }
        hits=$(grep -rnE --include='*.rs' "#\[path[[:space:]]*=.*$(basename "$f")" "$d" 2> /dev/null | sed "s|^$d/||")
        [ -z "$hits" ] || { echo "  HOLD: #[path] names held $f: $hits"; return 1; }
        ;;
    esac
    hits=$(grep -rnE --include='Cargo.toml' "^[^#]*path[[:space:]]*=.*$(basename "$f")" "$d" 2> /dev/null | sed "s|^$d/||")
    [ -z "$hits" ] || { echo "  HOLD: Cargo target names held $f: $hits"; return 1; }
  done < <(cd "$d" && grep -rlE --exclude-dir=target --exclude=HOLD.md '^// HOLD H[0-9]+:' . | sed 's|^\./||' | while IFS= read -r g; do head -n 1 "$g" | grep -qE '^// HOLD H[0-9]+:' && echo "$g"; done)
}
plant_hold_orphan() {  # a held-looking file with no index row
  local f="$PFX/src/fold_pr1_plant.rs" rc
  printf '// HOLD H999: plant\n// pub fn x() {}\n' > "$plant/$f"
  check_HOLD "$plant"; rc=$?
  restore "$f"; return $rc
}

# --- driver ------------------------------------------------------------------------------------------------------
bad=0
row() {   # row ID REAL_CMD... -- PLANT_CMD...   (plant side runs with the planted copy in place)
  local id=$1 r p; shift
  local -a realcmd=() plantcmd=(); local side=real
  for a in "$@"; do if [ "$a" = "--" ]; then side=plant; elif [ $side = real ]; then realcmd+=("$a"); else plantcmd+=("$a"); fi; done
  if "${realcmd[@]}"; then r=PASS; else r=FAIL; fi
  if "${plantcmd[@]}"; then p=MISSED; else p=DETECTED; fi
  echo "ROW $id real=$r plant=$p"
  [ $r = PASS ] && [ $p = DETECTED ] || bad=1
}
# Plants: each edits $plant, runs the row's check there, then restores the file.
plant_file() {  # FILE SED_EXPR CHECK ARGS...
  local f=$1 e=$2 rc; shift 2
  [ -f "$plant/$f" ] || return 1   # nothing to plant into: the check cannot be exercised, count it as MISSED
  sed -i "$e" "$plant/$f"
  "$@"; rc=$?
  restore "$f"; return $rc
}
plant_add_line() {  # FILE LINE CHECK ARGS...
  local f=$1 l=$2 rc; shift 2
  [ -f "$plant/$f" ] || return 1
  printf '%s\n' "$l" >> "$plant/$f"
  "$@"; rc=$?
  restore "$f"; return $rc
}
plant_new_file() {  # FILE CHECK ARGS...
  local f=$1 rc; shift
  mkdir -p "$(dirname "$plant/$f")"; printf 'x\n' > "$plant/$f"
  "$@"; rc=$?
  restore "$f"; return $rc
}
plant_P() {
  local c
  c=$(git commit-tree "$HEAD_SHA^{tree}" -p "$HEAD_SHA" -m "fold-pr1 plant: commit without the agent trailer") || return 1
  check_P "$c"
}
B=build; [ "$NO_BUILD" = 1 ] && B=nobuild

row A check_A "$ZB_TIP" -- check_A 0000000000000000000000000000000000000001
row B check_B "$HEAD_SHA" "$REMOVED" -- check_B "$HEAD_SHA" "rust-toolchain.toml"
row C check_C "$real" -- plant_file "$PFX/Cargo.toml" '/^zenith-float-num[[:space:]]*=/a zenith-float = { path = ".." }' check_C "$plant"
row D check_D "$real" -- plant_file "$PFX/Cargo.toml" 's|path = "../zenith-float-num"|git = "https://example.invalid/zf", rev = "0"|' check_D "$plant"
row E check_E "$real" -- plant_file Cargo.toml "s|^exclude = \[\"nostd-tests\"\]|exclude = [\"nostd-tests\", \"$PFX\"]|" check_E "$plant"
row F check_F "$real" "$B" -- plant_file Cargo.toml "/^$PFX[[:space:]]*=/s/, optional = true//" check_F "$plant" nobuild
row G check_G "$real" "$B" -- plant_file Cargo.toml 's/^default = .*/default = ["std"]/' check_G "$plant" nobuild
row K check_C "$real" -- plant_file "$PFX/Cargo.toml" '/^zenith-float-num[[:space:]]*=/a zenith-float = { path = "..", default-features = false }' check_C "$plant"
row L check_L "$real" "$B" -- plant_add_line src/lib.rs 'pub use zenith_float_bound::ball::Ball as BoundBall;' check_L "$plant" nobuild
row O check_O "$HEAD_SHA" -- check_O "$BASE"
row P check_P "$HEAD_SHA" -- plant_P
row Q check_Q "$real" -- plant_new_file "$PFX/rust-toolchain.toml" check_Q "$plant"
row R check_R "$real" -- plant_file .github/CODEOWNERS "\\|^/$PFX/tools/digests.sh|d" check_R "$plant"
row S check_S "$real" -- plant_add_line Cargo.lock $'\n[[package]]\nname = "fold-pr1-plant"\nversion = "0.0.1"' check_S "$plant"
row T check_T "$real" -- plant_file "$PFX/Cargo.toml" 's|^zenith-float-num[[:space:]]*=.*|&\n# zenith-float-num rev 627fedd11d5176a4299b0972e5aa723019e6377d|' check_T "$plant"
row Y check_Y "$HEAD_SHA" -- check_Y "$BASE"
row NOSTR check_NOSTR "$real" -- plant_add_line README.md "Formerly $OLD_NAME." check_NOSTR "$plant"
row HOLD check_HOLD "$real" -- plant_add_line "$PFX/src/lib.rs" 'pub mod carlson;' check_HOLD "$plant"
row HOLDIDX check_HOLD "$real" -- plant_hold_orphan
row NOCYCLE check_NOCYCLE "$real" -- plant_add_line zenith-float-num/src/lib.rs 'use zenith_float_bound as _;' check_NOCYCLE "$plant"

if [ $bad = 0 ]; then echo "fold-pr1: RESULT PASS (every row PASS + DETECTED) head=$HEAD_SHA"; exit 0; fi
echo "fold-pr1: RESULT FAIL head=$HEAD_SHA"; exit 1
