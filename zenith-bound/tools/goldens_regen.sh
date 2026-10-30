#!/usr/bin/env bash
# Golden regeneration gate (TEST_PLAN section 9): the goldlib self-test, then every generator with `--check DIR` into a fresh temp dir
# under $(zb_scratch) (byte-compared with tests/goldens), then `sha256sum -c tests/goldens/SHA256SUMS`. Uses the oracle venv ~/.local/oracle-venv/bin/python (scripts/goldlib.py require_oracle). No network.
# Temp files go under $(zb_scratch). `--slices-only` runs only the per-slice checks (no oracle; tools/goldens_regen_selftest.sh uses it).
# Per-slice SHA256SUMS: entries in strict sha256sum format, or, with no entries, an explicit `# EMPTY` marker line (never both).
# Untracked *_heavy.txt goldens (.gitignored; hashes pinned elsewhere) are not demanded in SHA256SUMS. Nothing under the repo is written. To add a generator, add a line to the GENERATORS list.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap goldens_regen "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
# PYTHONNOUSERSITE=1: nothing from ~/.local outside the venv (its 'flint' is an unrelated tool). GOLDENS_PYTHON may override the interpreter,
# but goldlib's startup check (four exact pins, sys.prefix is the venv, no ~/.local on sys.path) still applies, so only the venv passes.
export PYTHONNOUSERSITE=1
root=${GOLDENS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
PY=${GOLDENS_PYTHON:-$HOME/.local/oracle-venv/bin/python}
GENERATORS="scripts/gen_goldens.py scripts/gen_ball_goldens.py scripts/gen_ball_enclosure.py"
# stdout generators (no --check): name and committed file
STDOUT_GENERATORS="scripts/gen_review_goldens.py:mag_ops_review.txt"
self="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
slices_only=0; regen=0
for arg in "$@"; do case $arg in --slices-only) slices_only=1 ;; --regen) regen=1 ;; *) echo "goldens-regen: usage: tools/goldens_regen.sh [--slices-only] [--regen]"; exit 2 ;; esac; done
# Modes (Paul, 2026-10-24: oracle python ran flat out 10+ min in every gate). The key is a sha256 of every scripts/*.py, every
# committed golden input (tests/goldens/**: *.txt except untracked *_heavy.txt, SHA256SUMS, HEAVY.SHA256, GENERATORS), this script,
# and the oracle (interpreter path and venv *.dist-info names). It is stamped under $GOLDENS_KEY_DIR (default $(zb_scratch)/goldens-regen-keys)
# after a passing python run.
#   cached      key equals the stamp: python is skipped entirely; only the pinned SHA256SUMS / HEAVY.SHA256 are verified.
#   regen-gate  key changed (or no stamp): every generator runs --check; heavy rows are checked by HEAVY.SHA256 pin only (--skip-heavy).
#   regen-full  --regen: every generator runs --check and every heavy row is recomputed against its pin (about 14 min).
# ZB_GATE_MODE=merge (merge points) REFUSES (exit 2) without --regen, and with --slices-only (the merge gate runs every generator;
# Spence LA4 nit, 2026-10-25). ZB_GATE_MODE is recorded first thing, then the mode and its reason are logged.
# The key also covers scripts/requirements.lock (the oracle pins).
echo "goldens-regen: ZB_GATE_MODE=${ZB_GATE_MODE:-regular}"   # recorded first, also before an unknown-mode refusal
case "${ZB_GATE_MODE:-regular}" in regular|merge) ;; *) echo "goldens-regen: REFUSING: unknown ZB_GATE_MODE=${ZB_GATE_MODE} (want regular or merge)"; exit 2 ;; esac
gate_mode=${ZB_GATE_MODE:-regular}
if [ "$gate_mode" = merge ] && [ "$regen" != 1 ]; then echo "goldens-regen: REFUSING: ZB_GATE_MODE=merge needs --regen (the merge gate must regenerate and recompute every heavy row)"; exit 2; fi
if [ "$gate_mode" = merge ] && [ "$slices_only" = 1 ]; then echo "goldens-regen: REFUSING: ZB_GATE_MODE=merge does not allow --slices-only (the merge gate must run every generator, not only the per-slice checks)"; exit 2; fi
cd "$root" || exit 2
tmp=$(mktemp -d "$(zb_scratch)"/zb-goldens-regen.XXXXXX) || { echo "goldens-regen: mktemp failed"; exit 2; }; runlog_on_exit "rm -rf \"$tmp\""
# Pop direct-run guard (~/pop-logs/bin/popguard.sh, 11:49 2026-10-24: a bare 58 GB generator run global-OOMed Pop): before the
# first generator, refuse (exit 3) outside a capped popslot/poprun memory scope. Hosts without popguard.sh (CI) run unguarded, with a WARNING.
POPGUARD=${POPGUARD:-$HOME/pop-logs/bin/popguard.sh}; pg_ok=0
# POPGUARD pointing at nothing must not switch the guard off (Spence, 2026-10-24): on a host that has ~/pop-logs/bin/popguard.sh
# it REFUSES (exit 2); only a host with no popguard.sh at all runs unguarded, with a loud WARNING on stderr.
pg_require() {
  [ "$pg_ok" = 1 ] && return 0
  if [ -x "$POPGUARD" ]; then "$POPGUARD" goldens_regen || exit 3
  elif [ -e "$HOME/pop-logs/bin/popguard.sh" ]; then
    echo "goldens-regen: REFUSING: POPGUARD=$POPGUARD is missing or not executable, but this host has $HOME/pop-logs/bin/popguard.sh (the memory guard cannot be switched off)" >&2; exit 2
  else echo "goldens-regen: WARNING: no popguard.sh on this host ($POPGUARD): running WITHOUT the Pop memory guard" >&2; fi
  pg_ok=1
}
key_inputs() {
  { find scripts -type f -name '*.py' 2> /dev/null
    find scripts -maxdepth 1 -type f -name requirements.lock 2> /dev/null
    find tests/goldens -type f \( -name '*.txt' -o -name SHA256SUMS -o -name HEAVY.SHA256 -o -name GENERATORS \) ! -name '*_heavy.txt' 2> /dev/null
  } | sort | while IFS= read -r f; do sha256sum "$f"; done
  echo "tool $(sha256sum < "$self")"; echo "slices_only=$slices_only"; echo "py $(readlink -f "$PY" 2> /dev/null || echo "$PY")"
  ls -d "$(dirname "$(dirname "$PY")")"/lib/python*/site-packages/*.dist-info 2> /dev/null | sed 's|.*/||' | sort
}
key=$(key_inputs | sha256sum | cut -c1-64); short=$(printf '%s' "$key" | cut -c1-16)
keydir=${GOLDENS_KEY_DIR:-$(zb_scratch)/goldens-regen-keys}; stamp="$keydir/$(printf '%s' "$root" | sha256sum | cut -c1-16).key"
if [ "$regen" = 1 ]; then mode=full; heavy_args=""
  if [ "$gate_mode" = merge ]; then echo "goldens-regen: mode=regen-full (--regen; gate mode=merge); heavy rows recomputed"; else echo "goldens-regen: mode=regen-full (--regen); heavy rows recomputed"; fi
elif [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$key" ]; then mode=cached
  echo "goldens-regen: mode=cached (key $short unchanged since the last passing run); python skipped, pinned sums verified only"
else mode=gate; heavy_args="--skip-heavy"
  if [ -f "$stamp" ]; then why="inputs differ from the last passing run"; else why="no passing run recorded"; fi
  echo "goldens-regen: mode=regen-gate (key $short changed: $why); heavy rows by HEAVY.SHA256 pin only"
fi
rc=0
if [ "$slices_only" = 0 ] && [ "$mode" = cached ]; then
(cd tests/goldens && sha256sum --strict -c SHA256SUMS) || { echo "goldens-regen: sha256sum -c FAILED"; rc=1; }
elif [ "$slices_only" = 0 ]; then
pg_require
runlog_step goldlib-selftest -- "$PY" scripts/goldlib.py --selftest || rc=1
for g in $GENERATORS; do
  runlog_step "$(basename "$g" .py)" -- "$PY" "$g" --check "$tmp/$(basename "$g" .py)" || { echo "goldens-regen: $g FAILED"; rc=1; }
done
for e in $STDOUT_GENERATORS; do
  g=${e%%:*}; f=${e#*:}
  { "$PY" "$g" > "$tmp/$f" && cmp "$tmp/$f" "tests/goldens/$f" && echo "check: OK $f"; } || { echo "goldens-regen: $g FAILED"; rc=1; }
done
(cd tests/goldens && sha256sum --strict -c SHA256SUMS) || { echo "goldens-regen: sha256sum -c FAILED"; rc=1; }
fi
# per-slice golden dirs (S0, 2026-10-24): tests/goldens/slices/<slice>/ with its own GENERATORS and SHA256SUMS, so slices never share a file
for d in tests/goldens/slices/*/; do
  [ -d "$d" ] || continue
  s=$(basename "$d")
  [ -f "$d/SHA256SUMS" ] && [ -f "$d/GENERATORS" ] || { echo "goldens-regen: slice $s lacks SHA256SUMS or GENERATORS"; rc=1; continue; }
  if [ -f "$d/HEAVY.SHA256" ]; then   # pinned heavy sums: strict format in every mode (cheap; python checks names in regen modes)
    hbad=$(grep -v -E '^[[:space:]]*(#|$)' "$d/HEAVY.SHA256" | grep -v -E '^[0-9a-f]{64}  [^/[:space:]]+$' || true)
    [ -z "$hbad" ] || { echo "goldens-regen: slice $s: malformed HEAVY.SHA256 line(s): $(printf '%s\n' "$hbad" | head -3)"; rc=1; }
  fi
  [ "$mode" = cached ] || for g in $(grep -v -E '^[[:space:]]*(#|$)' "$d/GENERATORS"); do
    pg_require
    runlog_step "$s-$(basename "$g" .py)" -- "$PY" "$g" --check "$tmp/slice-$s-$(basename "$g" .py)" $heavy_args || { echo "goldens-regen: $g FAILED"; rc=1; }
  done
  entries=$(grep -v -E '^[[:space:]]*(#|$)' "$d/SHA256SUMS" || true)   # no entries is a valid state (rc 1): no ERR trap line
  marked=0; if grep -q -E '^# EMPTY([: ]|$)' "$d/SHA256SUMS"; then marked=1; fi
  if [ -z "$entries" ] && [ "$marked" = 1 ]; then
    echo "goldens-regen: slice $s: SHA256SUMS marked '# EMPTY', skipped"
  elif [ -z "$entries" ]; then
    echo "goldens-regen: slice $s: SHA256SUMS has neither entries nor a '# EMPTY' marker line"; rc=1
  else
    [ "$marked" = 0 ] || { echo "goldens-regen: slice $s: SHA256SUMS is marked '# EMPTY' but lists entries"; rc=1; }
    # grep -v rc 1 = every line valid (normal; must not trip ERR); rc 2 = grep itself failed (loud).
    grc=0; bad=$(printf '%s\n' "$entries" | grep -v -E '^[0-9a-f]{64}  [^/[:space:]]+$') || grc=$?
    case $grc in 0|1) ;; *) echo "goldens-regen: slice $s: grep failed (rc $grc) checking SHA256SUMS lines"; rc=1 ;; esac
    bad=$(printf '%s\n' "$bad" | head -3)
    [ -z "$bad" ] || { echo "goldens-regen: slice $s: malformed SHA256SUMS line(s): $bad"; rc=1; }
    (cd "$d" && printf '%s\n' "$entries" | sha256sum --strict -c -) || { echo "goldens-regen: slice $s sha256sum -c FAILED"; rc=1; }
  fi
  for f in "$d"*.txt; do
    [ -f "$f" ] || continue
    case "$f" in *_heavy.txt) continue ;; esac   # heavy goldens: .gitignored, pinned outside SHA256SUMS
    if [ "$marked" = 1 ]; then echo "goldens-regen: slice $s: SHA256SUMS marked '# EMPTY' but has $(basename "$f")"; rc=1; continue; fi
    grep -q -E "^[0-9a-f]{64}  $(basename "$f")\$" "$d/SHA256SUMS" || { echo "goldens-regen: slice $s: $(basename "$f") is not in its SHA256SUMS"; rc=1; }
  done
done
# every committed golden file must be listed in SHA256SUMS
[ "$slices_only" = 1 ] || for f in tests/goldens/*.txt; do
  grep -q -E "^[0-9a-f]{64}  $(basename "$f")\$" tests/goldens/SHA256SUMS || { echo "goldens-regen: $(basename "$f") is not in SHA256SUMS"; rc=1; }
done
if [ "$rc" = 0 ] && [ "$mode" != cached ]; then mkdir -p "$keydir" && printf '%s\n' "$key" > "$stamp.tmp" && mv "$stamp.tmp" "$stamp" && echo "goldens-regen: key $short stamped (next unchanged run is cached)"; fi
[ "$rc" = 0 ] && echo "goldens-regen: ok" || echo "goldens-regen: FAILED"
exit "$rc"
