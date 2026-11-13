#!/usr/bin/env bash
# No tool hardcodes /tmp (Jeff 2026-10-24: a reboot wiped /tmp and with it work dirs, a dependency pin and a commit). Temp files and
# dirs go under "$(zb_scratch)" (tools/scratch.sh). Scans the code lines (not comments) of tools/*.sh for `/tmp`; the one exception is
# the shared cross-agent locks /tmp/pop-deep.lock and /tmp/pop-gate.lock (a lock file, not a work dir; moving it would split the lock between agents).
# Then proves the scan can fail (a fixture with /tmp is flagged) and does not flag a comment. usage: tools/no_tmp_selftest.sh
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap no_tmp_selftest "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
here="$(cd "$(dirname "$0")" && pwd)"; . "$here/scratch.sh"
fail() { echo "no-tmp selftest FAIL: $*"; exit 1; }
t=$(mktemp -d "$(zb_scratch)/no-tmp-selftest.XXXXXX") || fail "mktemp failed"; runlog_on_exit "rm -rf \"$t\""; . "$here/selftest_cases.sh"; st_init "$t"
scan() { # scan DIR: prints file:line for each code line naming /tmp
  grep -HnE '(^|[ "'"'"'=:(-])/tmp([/ "'"'"'$}]|$)' "$1"/*.sh | grep -vE '^[^:]+:[0-9]+: *#' | grep -vE 'pop-(deep|gate)\.lock|/no_tmp_selftest\.sh:'
}
ck; hits=$(scan "$here") && fail "tools hardcode /tmp (use \"\$(zb_scratch)\" from tools/scratch.sh):
$hits"
mkdir -p "$t/fx"; printf '#!/usr/bin/env bash\n# a comment naming /tmp is fine\nd=$CARGO_TARGET_DIR/tmp/x\n' > "$t/fx/ok.sh"
ck; scan "$t/fx" > /dev/null && fail "the scan flagged a comment line or a \$dir/tmp path"
printf 'x=$(mktemp -d /tmp/foo.XXXXXX)\n' > "$t/fx/bad.sh"
ck; grep -q '^'"$t"'/fx/bad.sh:1:' <<< "$(scan "$t/fx")" || fail "the scan missed a hardcoded /tmp"
printf 'y=$(mktemp -d "${TMPDIR:-/tmp}/foo.XXXXXX")\n' > "$t/fx/bad2.sh"   # e83631d goldens_regen: /tmp then } was not matched
ck; grep -q '^'"$t"'/fx/bad2.sh:1:' <<< "$(scan "$t/fx")" || fail "the scan missed \${TMPDIR:-/tmp}"
ck; [ "$(ZB_SCRATCH="$t/s" zb_scratch)" = "$t/s" ] && [ -d "$t/s" ] || fail "zb_scratch ignores ZB_SCRATCH or does not create it"
ck; case "$(env -u ZB_SCRATCH HOME="$t/h" bash -c ". '$here/scratch.sh'; zb_scratch")" in "$t/h/"*) ;; *) fail "zb_scratch default is not under \$HOME" ;; esac
st_done no_tmp 6   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "no-tmp selftest: ok"
