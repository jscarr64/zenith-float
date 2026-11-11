#!/usr/bin/env bash
# K11 (ACCEPTANCE-pr39; Paul's ruling 2026-11-10 23:15 CT): shipped scripts carry no path on one machine, and the oracle
# generators find their interpreter with `#!/usr/bin/env python3` and then FAIL CLOSED unless it is the pinned oracle env.
#   1. static: every tracked python script with a shebang uses `#!/usr/bin/env python3`; no tracked script (*.py *.sh *.rs
#      *.toml) names /home/jsc64 or /mnt/hub (Spence's K11 regex).
#   2. forged env: a throwaway team env named like the pinned one (oracle-2026-10-30) but with mpmath 1.2.1, first on PATH.
#      Every zenith-float-bound/scripts generator, executed through its shebang, must exit nonzero with the pin-mismatch
#      message, before any argument parsing or work.
#   3. control: the same forged env with the four pinned versions passes `goldlib.py --team-env` (so step 2 fails on the
#      version, not on the fixture).
# Nothing is written outside a mktemp dir; HOME points at it, so no real team env, popguard or ~/.local is consulted.
set -u
repo="$(cd "$(dirname "$0")/.." && pwd)"; cd "$repo" || exit 1
t="$(mktemp -d)"; trap 'rm -rf "$t"' EXIT
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }

# 1. static
while IFS= read -r f; do
  first=$(head -1 "$f")
  case "$first" in '#!'*python*) n=$((n + 1)); [ "$first" = "#!/usr/bin/env python3" ] || bad "$f shebang is '$first', want '#!/usr/bin/env python3'";; esac
done < <(git ls-files '*.py')
hits=$(git grep -nE '/home/jsc64|/mnt/hub' -- '*.py' '*.sh' '*.rs' '*.toml' || true)
[ -z "$hits" ] || bad "machine paths in shipped scripts:
$hits"

# forged team env (HOME=$t): $t/team-env/oracle-2026-10-30 is a real venv with stub oracle modules
mkenv() {  # mkenv MPMATH_VERSION [VENV_DIR]  (default: the pinned team env name under $t/team-env)
  rm -rf "$t/team-env" "$t/.local"; mkdir -p "$t/team-env" "$t/pop-logs/bin"
  printf '#!/bin/sh\nexit 0\n' > "$t/team-env/verify_env.sh"; printf '#!/bin/sh\nexit 0\n' > "$t/pop-logs/bin/popguard.sh"
  chmod +x "$t/team-env/verify_env.sh" "$t/pop-logs/bin/popguard.sh"
  v="${2:-$t/team-env/oracle-2026-10-30}"; mkdir -p "$(dirname "$v")"; ENVBIN="$v/bin"; /usr/bin/python3 -m venv --without-pip "$v" || { bad "cannot create the forged venv"; return 1; }
  sp=$("$v/bin/python3" -c 'import site; print(site.getsitepackages()[0])')
  printf '__version__ = "%s"\niv = None\n' "$1" > "$sp/mpmath.py"
  printf '__version__ = "1.14.0"\n' > "$sp/sympy.py"
  printf 'def version():\n    return "2.3.2"\n__version__ = "2.3.2"\nmpq = object\n' > "$sp/gmpy2.py"
  printf '__version__ = "0.9.0"\narb = object\n' > "$sp/flint.py"
}
run_forged() { env -i HOME="$t" PATH="$ENVBIN:/usr/bin:/bin" PYTHONNOUSERSITE=1 "$@" 2>&1; }

sweep() {  # sweep REGEX LABEL: every generator must exit nonzero with output matching REGEX
for s in zenith-float-bound/scripts/*.py; do
  first=$(head -1 "$s")
  case "$first" in '#!'*) ;; *) continue ;; esac
  # never execute a script whose shebang could reach a real oracle env (it would do real work); step 1 already failed it
  [ "$first" = "#!/usr/bin/env python3" ] || continue
  n=$((n + 1)); arg=--no-such-flag; [ "$s" = zenith-float-bound/scripts/goldlib.py ] && arg=--team-env
  rel="./${s#zenith-float-bound/}"; how=()   # executable: through its shebang; otherwise `python3 FILE` (the same PATH lookup)
  [ -x "$s" ] || how=(python3)
  out=$(cd zenith-float-bound && run_forged "${how[@]}" "$rel" "$arg"); rc=$?
  [ "$rc" != 0 ] || bad "$s ran under a forged env ($2) with rc 0: $out"
  grep -qE "$1" <<< "$out" || bad "$s: no $2 refusal (rc $rc): $(tail -2 <<< "$out")"
done
}

# 2. wrong pin: every generator refuses
mkenv 1.2.1 || exit 1
sweep 'mpmath is 1\.2\.1, want 1\.3\.0|oracle pin mismatch: .*mpmath 1\.2\.1' "pin-mismatch (forged mpmath 1.2.1)"

# 4. right pins, but not a team env (Jeff 2026-10-30: never ~/.local/oracle-venv): every generator refuses
mkenv 1.3.0 "$t/.local/oracle-venv" || exit 1
sweep 'REFUSED\[(oracle-venv|not-team-env)\]' "team-env (pinned versions in \$HOME/.local/oracle-venv)"

# 3. control: right pins in the team env pass the env check
mkenv 1.3.0 || exit 1
out=$(cd zenith-float-bound && run_forged python3 ./scripts/goldlib.py --team-env); rc=$?
n=$((n + 1)); [ "$rc" = 0 ] && grep -q '^team-env=oracle-2026-10-30$' <<< "$out" || bad "control: pinned forged env refused (rc $rc): $(tail -2 <<< "$out")"

[ "$fail" = 0 ] && echo "oracle_env_selftest: ok ($n checks)" || { echo "oracle_env_selftest: FAILED"; exit 1; }
