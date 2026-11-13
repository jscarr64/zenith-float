#!/usr/bin/env bash
# Format-1 preservation check for the format-2 LA golds (Harper 2026-10-29, Spence's review of harper/zb-la-golds d9cfacbe).
# For la_logdet.txt and la_solve.txt, compares the working-tree file with the format-1 file at BASE (default 8cd0298):
#   data rows: same count and order; a refusal row (ERR ...) is byte-identical; any other row is the format-1 row byte for byte,
#              then one space and exactly six more tokens (cond_row k_row scale_row containment_only reason rref);
#   comment lines: every format-1 comment line is kept, in order, byte-identical, except the '# columns' line, which may
#              only gain the six cond column names; added comment lines must be cond/format docs (allowed prefixes below).
# Prints every changed and added comment line. usage: tools/la_format1_preserved.sh [BASE [ROOT]]. Exit 0 ok, 1 violation, 2 setup.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap la_format1_preserved "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
base=${1:-8cd0298}
root=${2:-$(cd "$(dirname "$0")/.." && pwd)}
. "${BASH_SOURCE[0]%/*}/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/la-format1.XXXXXX") || exit 2; runlog_on_exit "rm -rf \"$t\""
rc=0
for p in tests/goldens/slices/cholesky/la_logdet.txt tests/goldens/slices/lu/la_solve.txt; do
  git -C "$root" show "$base:$p" > "$t/old" || { echo "la_format1_preserved: cannot read $base:$p"; exit 2; }
  [ -f "$root/$p" ] || { echo "la_format1_preserved: missing $p"; exit 2; }
  awk -v name="$p" '
    function allowed(s) { return s ~ /^# (format: 2 |cond columns: |cond_row: |k_row: |scale_row: |containment_only reason: |rref: |tightness \(replay\): |arb cross-check: )/ }
    FNR == NR { if ($0 ~ /^#/) oc[++no] = $0; else od[++nd] = $0; next }
    /^#/ { nc[++mc] = $0; next }
    { m++; if (m > nd) { printf "%s: FAIL extra data row %d: %s\n", name, m, $1; bad++; next }
      o = od[m]
      if (o ~ / ERR /) { if ($0 != o) { printf "%s: FAIL refusal row %d (%s) changed\n", name, m, $1; bad++ } ; next }
      if (substr($0, 1, length(o) + 1) != o " ") { printf "%s: FAIL row %d (%s): format-1 prefix not byte-identical\n", name, m, $1; bad++; next }
      k = split(substr($0, length(o) + 2), x, " "); if (k != 6) { printf "%s: FAIL row %d (%s): %d added tokens, want 6\n", name, m, $1, k; bad++ }
      else ext++ }
    END {
      if (m != nd) { printf "%s: FAIL %d data rows, format 1 had %d\n", name, m, nd; bad++ }
      j = 1
      for (i = 1; i <= mc; i++) {
        if (j <= no && nc[i] == oc[j]) { j++; continue }
        if (j <= no && oc[j] ~ /^# columns/ && nc[i] ~ /^# columns/) {
          s = oc[j]; sub(/ \| ERR /, " cond_row k_row scale_row containment_only reason rref | ERR ", s)
          if (nc[i] == s) { printf "%s: changed comment (cond names added): %s\n", name, nc[i]; chg++; j++; continue }
          printf "%s: FAIL columns line changed beyond the cond names: %s\n", name, nc[i]; bad++; j++; continue }
        if (allowed(nc[i])) { printf "%s: added comment: %s\n", name, substr(nc[i], 1, 60); add++; continue }
        printf "%s: FAIL comment line not in format 1 and not a cond/format doc: %s\n", name, nc[i]; bad++ }
      if (j <= no) { printf "%s: FAIL format-1 comment line dropped or reordered: %s\n", name, oc[j]; bad++ }
      printf "%s: %d data rows (%d with 6 cond columns, %d refusals unchanged), %d comment lines changed, %d added, %d violations\n", name, m, ext, nd - ext, chg, add, bad + 0
      exit (bad > 0) }' "$t/old" "$root/$p" || rc=1
done
[ "$rc" = 0 ] && echo "la_format1_preserved: ok (base $base)" || echo "la_format1_preserved: FAIL (base $base)"
exit "$rc"
