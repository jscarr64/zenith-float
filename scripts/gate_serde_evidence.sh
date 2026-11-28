#!/usr/bin/env bash
# Gate evidence for the features-random-serde step (Jeff 2026-11-24 9:06 AM CT: "Go: in the gate, prepublish skips its
# random,serde rerun only after checking features-random-serde passed rc=0 on the same SHA; standalone prepublish stays
# complete"). Selftest: tools/test_prepublish_serde_dedupe.sh.
#   run FILE CMD...  (tools/full_gate.sh features-random-serde step) remove FILE, run CMD, then write FILE: this gate run's
#                    ZF_GATE_RUN_ID, the HEAD SHA, CMD's rc, and clean=yes only if HEAD did not move and the tree was clean
#                    (git status --porcelain empty) before and after CMD. Exits with CMD's rc.
#   check FILE       (scripts/zenith_prepublish.sh, gate mode) exit 0 and print why only if FILE was written by this gate
#                    run (run id matches), rc=0, clean=yes, its SHA is HEAD now and the tree is clean now; else exit 1 and
#                    print why the serde suite must run.
set -u
clean() { [[ -z "$(git status --porcelain 2>/dev/null)" ]] && git rev-parse HEAD >/dev/null 2>&1; }
field() { sed -n "s/^$1=//p" "$2" | head -1; }
case "${1:-}" in
  run)
    file=${2:?usage: gate_serde_evidence.sh run FILE CMD...}; shift 2
    rm -f "$file"
    sha0=$(git rev-parse HEAD 2>/dev/null || echo none); c0=no; clean && c0=yes
    "$@"; rc=$?
    sha1=$(git rev-parse HEAD 2>/dev/null || echo none); c1=no; clean && c1=yes
    ok=no; [[ $c0 == yes && $c1 == yes && $sha0 == "$sha1" ]] && ok=yes
    printf 'step=features-random-serde\nrun_id=%s\nsha=%s\nrc=%s\nclean=%s\n' "${ZF_GATE_RUN_ID:-}" "$sha0" "$rc" "$ok" > "$file"
    exit "$rc" ;;
  check)
    file=${2:?usage: gate_serde_evidence.sh check FILE}
    [[ -f "$file" ]] || { echo "no features-random-serde evidence ($file)"; exit 1; }
    id=$(field run_id "$file"); sha=$(field sha "$file"); rc=$(field rc "$file"); c=$(field clean "$file")
    head=$(git rev-parse HEAD 2>/dev/null || echo none)
    [[ -n "${ZF_GATE_RUN_ID:-}" && "$id" == "$ZF_GATE_RUN_ID" ]] || { echo "evidence is not from this gate run (run_id '$id')"; exit 1; }
    [[ "$rc" == 0 ]] || { echo "features-random-serde rc=$rc"; exit 1; }
    [[ "$c" == yes ]] || { echo "tree dirty or HEAD moved during features-random-serde"; exit 1; }
    [[ "$sha" == "$head" ]] || { echo "features-random-serde ran on $sha, HEAD is $head"; exit 1; }
    clean || { echo "tree is dirty now"; exit 1; }
    echo "features-random-serde passed rc=0 on $sha (clean tree, gate run $id)"; exit 0 ;;
  *) echo "usage: gate_serde_evidence.sh run FILE CMD... | check FILE" >&2; exit 2 ;;
esac
