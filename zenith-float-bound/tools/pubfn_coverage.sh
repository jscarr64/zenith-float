#!/usr/bin/env bash
# Public-API coverage list (TEST_PLAN section 11). grep and awk only (no rg on Pop), deterministic, no network.
#
# usage: tools/pubfn_coverage.sh [--list]   (checks the checkout this script lives in; --list prints the items and exits 0)
#   env: PUBFN_ROOT=DIR     repo root to check (default: the parent of tools/)
#        PUBFN_EXEMPT=FILE  exemption list (default: $PUBFN_ROOT/tools/pubfn_coverage_exempt.txt)
#
# ITEMS (the public API), from src/**/*.rs and, with the prefix `common::`, tests/common/**/*.rs:
#   every `pub fn` / `pub const fn`     as  Type::name  (enclosing impl), or  mod::path::name  (free fn; `crate::name` at the src root)
#   every `impl Trait for Type`         as  Trait for Type
#   every variant of a `pub enum`       as  Enum::Variant
#   (`pub(crate) fn` and friends are not public; consts and statics are not listed.)
# COVERS: every .rs file under tests/ (tests/common excluded) may carry header lines  `//! COVERS: name, name, ...`.
# FAILS (exit 1, one line each) when:
#   UNCOVERED  an item is in no COVERS line and not in the exemption list
#   STALE      a COVERS name is not an item (renamed or removed)
#   STALE-EXEMPT  an exemption-list entry is not an item (renamed or removed)
#   NOTEST     the declaring file has no #[test] fn whose text mentions the name (word match on the last `::` segment,
#              or on the Trait / Type of an `impl Trait for Type` name)
# Parser limits (rustfmt layout assumed): an impl/mod/enum block opens with `{` and its members start on later lines; macro-generated items are not seen.
# Name coverage is deliberately weak; the mutation matrix (tools/mutate.sh) is the strong check.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap pubfn_coverage "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
root=${PUBFN_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
exempt=${PUBFN_EXEMPT:-$root/tools/pubfn_coverage_exempt.txt}
cd "$root" || { echo "pubfn_coverage: cannot cd to $root"; exit 2; }
tmp=$(mktemp -d); runlog_on_exit "rm -rf \"$tmp\""

# ---- 1. items -------------------------------------------------------------------------------------------------------
extract() { # extract <file> <modpath> <prefix>
awk -v modpath="$2" -v prefix="$3" '
function strip(l) { sub(/\/\/.*/, "", l); gsub(/"([^"\\]|\\.)*"/, "\"\"", l); return l }
function count(s, ch,   n, i) { n = 0; for (i = 1; i <= length(s); i++) if (substr(s, i, 1) == ch) n++; return n }
function lastseg(s,   n, a) { n = split(s, a, "::"); return a[n] }
function ident(s) {  # first type-ish identifier of s
  sub(/^[ \t&]+/, "", s); sub(/^(dyn|mut)[ \t]+/, "", s)
  sub(/[<{ \t(;].*$/, "", s); return lastseg(s)
}
function join(a, b) { return (a == "") ? b : (b == "") ? a : a "::" b }
function emit(n) { print join(prefix, n) }
function modname(   i, m) { m = modpath; for (i = 1; i <= sp; i++) if (kind[i] == "mod") m = join(m, nm[i]); return m }
BEGIN { depth = 0; sp = 0; pend = "" }
{
  line = strip($0)
  d = depth
  # --- open a context on impl / mod / enum lines
  if (pend == "") {
    if (line ~ /^[ \t]*(pub(\([a-z]+\))?[ \t]+)?impl[ \t<]/) {
      s = line; sub(/^[ \t]*(pub(\([a-z]+\))?[ \t]+)?impl[ \t]*/, "", s)
      if (substr(s, 1, 1) == "<") {  # skip balanced generics
        lvl = 0
        for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); if (c == "<") lvl++; else if (c == ">") { lvl--; if (lvl == 0) { s = substr(s, i + 1); break } } }
      }
      sub(/[ \t]*where[ \t].*$/, "", s); sub(/[ \t]*\{.*$/, "", s)
      if (match(s, / for /)) {
        tr = ident(substr(s, 1, RSTART - 1)); ty = ident(substr(s, RSTART + 5))
        if (d == (sp ? dep[sp] + 1 : 0)) emit(tr " for " ty)
        pend = "impl"; pname = ty
      } else { pend = "impl"; pname = ident(s) }
      ctxpub = 1
    } else if (line ~ /^[ \t]*pub[ \t]+enum[ \t]/) {
      s = line; sub(/^[ \t]*pub[ \t]+enum[ \t]+/, "", s); pend = "enum"; pname = ident(s)
    } else if (line ~ /^[ \t]*(pub(\([a-z]+\))?[ \t]+)?mod[ \t]/ && line !~ /;[ \t]*$/) {
      s = line; sub(/^[ \t]*(pub(\([a-z]+\))?[ \t]+)?mod[ \t]+/, "", s); pend = "mod"; pname = ident(s)
    }
  }
  # --- members of the innermost context (depth before this line)
  incontext = (sp > 0 && d == dep[sp] + 1)
  if (line ~ /^[ \t]*pub[ \t]+(const[ \t]+)?fn[ \t]+[A-Za-z_]/ && (incontext || (sp == 0 && d == 0))) {
    s = line; sub(/^[ \t]*pub[ \t]+(const[ \t]+)?fn[ \t]+/, "", s); sub(/[^A-Za-z0-9_].*$/, "", s)
    if (incontext && kind[sp] == "impl") emit(nm[sp] "::" s)
    else if (incontext && kind[sp] == "mod") emit(join(modname(), s))
    else if (!incontext) { m = modname(); emit(m == "" ? (prefix == "" ? "crate::" s : s) : m "::" s) }
  } else if (incontext && kind[sp] == "enum" && ispub[sp] && line ~ /^[ \t]*[A-Z][A-Za-z0-9_]*[ \t]*([,({=]|$)/) {
    s = line; sub(/^[ \t]*/, "", s); sub(/[^A-Za-z0-9_].*$/, "", s); emit(nm[sp] "::" s)
  }
  # --- push a pending context at the first `{`, then update depth and pop closed contexts
  if (pend != "" && index(line, "{") > 0) { sp++; kind[sp] = pend; nm[sp] = pname; dep[sp] = d; ispub[sp] = 1; pend = "" }
  else if (pend != "" && line ~ /;[ \t]*$/) pend = ""
  depth += count(line, "{") - count(line, "}")
  while (sp > 0 && depth <= dep[sp]) sp--
}' "$1"
}

: > "$tmp/items.raw"
if [ -d src ]; then
  find src -name '*.rs' | sort | while read -r f; do
    rel=${f#src/}; rel=${rel%.rs}
    case "$rel" in lib|main) mp="";; */mod) mp=${rel%/mod};; mod) mp="";; *) mp=$rel;; esac
    extract "$f" "$(printf '%s' "$mp" | sed 's#/#::#g')" ""
  done >> "$tmp/items.raw"
fi
if [ -d tests/common ]; then
  find tests/common -name '*.rs' | sort | while read -r f; do
    rel=${f#tests/common/}; rel=${rel%.rs}
    case "$rel" in mod) mp="";; */mod) mp=${rel%/mod};; *) mp=$rel;; esac
    extract "$f" "$(printf '%s' "$mp" | sed 's#/#::#g')" "common"
  done >> "$tmp/items.raw"
fi
sort -u "$tmp/items.raw" > "$tmp/items"
if [ "${1:-}" = "--list" ]; then cat "$tmp/items"; exit 0; fi

# ---- 2. COVERS declarations: "<file>\t<name>" ---------------------------------------------------------------------------
: > "$tmp/covers"
if [ -d tests ]; then
  find tests -name '*.rs' -not -path 'tests/common/*' | sort | while read -r f; do
    awk -v file="$f" '
      /^\/\/![ \t]*COVERS:/ {
        s = $0; sub(/^\/\/![ \t]*COVERS:[ \t]*/, "", s); n = split(s, a, ",")
        for (i = 1; i <= n; i++) { g = a[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print file "\t" g }
      }' "$f"
  done | sort -u >> "$tmp/covers"
fi
cut -f2 "$tmp/covers" | sort -u > "$tmp/covered_names"

# ---- 3. exemptions -------------------------------------------------------------------------------------------------------
# Entries are the first tab field of each non-comment, non-blank line. awk, not grep -v: a list with no entries (the shipped
# state) is a normal case, and grep exits 1 on zero selected lines, which under pipefail made the load "fail" (ERR line 119).
# awk exits 0 on zero entries and non-zero only on a real read error, which stops the check.
if [ -f "$exempt" ]; then
  awk -F "\t" '$0 !~ /^[[:space:]]*(#|$)/ { print $1 }' "$exempt" | sort -u > "$tmp/exempt" \
    || { echo "pubfn_coverage: cannot read exempt list $exempt"; exit 2; }
else
  : > "$tmp/exempt"
fi
nex=$(wc -l < "$tmp/exempt" | tr -d ' ')
if [ "$nex" = 0 ]; then echo "pubfn_coverage: exempt list is empty ($exempt)"; else echo "pubfn_coverage: exempt list: $nex items ($exempt)"; fi

fail=0
# UNCOVERED
while IFS= read -r it; do
  grep -qxF -- "$it" "$tmp/covered_names" && continue
  grep -qxF -- "$it" "$tmp/exempt" && continue
  echo "FAIL UNCOVERED: $it is in no COVERS line"; fail=1
done < "$tmp/items"
# STALE-EXEMPT: an exempt entry that names no public item (renamed or removed) would otherwise be silently ignored
while IFS= read -r ex; do
  grep -qxF -- "$ex" "$tmp/items" && continue
  echo "FAIL STALE-EXEMPT: $ex is in the exempt list but is not a public item"; fail=1
done < "$tmp/exempt"
# STALE and NOTEST
while IFS="$(printf '\t')" read -r file name; do
  if ! grep -qxF -- "$name" "$tmp/items"; then
    echo "FAIL STALE: $file names $name, which is not a public item"; fail=1; continue
  fi
  case "$name" in
    *" for "*) tok1=${name%% for *}; tok2=${name##* for }; tok1=${tok1##*::}; tok2=${tok2##*::} ;;
    *) tok1=${name##*::}; tok2="" ;;
  esac
  ok=$(awk -v t1="$tok1" -v t2="$tok2" '
    function scan(s,   n, a, i) { gsub(/[^A-Za-z0-9_]/, " ", s); n = split(s, a, " ")
      for (i = 1; i <= n; i++) if (a[i] == t1 || (t2 != "" && a[i] == t2)) found = 1 }
    function strip(l) { sub(/\/\/.*/, "", l); gsub(/"([^"\\]|\\.)*"/, "\"\"", l); return l }
    function count(s, ch,   n, i) { n = 0; for (i = 1; i <= length(s); i++) if (substr(s, i, 1) == ch) n++; return n }
    BEGIN { want = 0; intest = 0; seen = 0; dd = 0; found = 0 }
    { raw = $0; line = strip(raw)
      if (!intest && line ~ /^[ \t]*#\[test\]/) { want = 1; next }
      if (want && !intest && line ~ /fn[ \t]/) { intest = 1; dd = 0; seen = 0; want = 0 }
      if (intest) { scan(raw); dd += count(line, "{") - count(line, "}"); if (index(line, "{") > 0) seen = 1; if (seen && dd <= 0) intest = 0 }
    }
    END { print found ? "yes" : "no" }' "$file")
  if [ "$ok" != yes ]; then
    echo "FAIL NOTEST: $file declares $name but has no #[test] fn mentioning $tok1${tok2:+ or $tok2}"; fail=1
  fi
done < "$tmp/covers"

n=$(wc -l < "$tmp/items" | tr -d ' ')
if [ "$fail" = 0 ]; then echo "pubfn_coverage: ok ($n public items)"; else echo "pubfn_coverage: FAILED ($n public items)"; fi
exit "$fail"
