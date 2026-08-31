#!/usr/bin/env bash
# Compare Criterion bench output against doc/bench-baselines.tsv (tab-separated, no JSON).
#
#   ./scripts/bench-compare.sh              # run --quick benches, compare to baseline
#   ./scripts/bench-compare.sh --update     # refresh doc/bench-baselines.tsv from a bench run
#   ./scripts/bench-compare.sh --update-from LOG  # refresh baseline from an existing log
#   ./scripts/bench-compare.sh --check FILE # compare an existing log to baseline
#
# Environment:
#   REGRESSION_PCT   allowed slowdown in wall time (default 10)
#   REGRESSION_FILE  baseline path (default doc/bench-baselines.tsv)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

baseline="${REGRESSION_FILE:-doc/bench-baselines.tsv}"
regression_pct="${REGRESSION_PCT:-10}"
tmp_log="$(mktemp)"
tmp_tsv="$(mktemp)"
trap 'rm -f "$tmp_log" "$tmp_tsv"' EXIT

usage() {
  sed -n '2,8p' "$0"
  exit "${1:-0}"
}

# Parse Criterion text lines into TSV: benchmark_id, time_ns, throughput_per_s (- if absent).
parse_bench_log() {
  local log="${1:?log file required}"
  awk '
    function to_ns(val, unit,    scale) {
      if (unit == "ns") scale = 1
      else if (unit == "µs" || unit == "us") scale = 1000
      else if (unit == "ms") scale = 1000000
      else if (unit == "s") scale = 1000000000
      else scale = 1
      return int(val * scale + 0.5)
    }
    function to_elems(val, unit,    scale) {
      if (unit ~ /^Melem/) scale = 1000000
      else if (unit ~ /^Kelem/) scale = 1000
      else scale = 1
      return int(val * scale + 0.5)
    }
    function parse_triple(line,    s, n) {
      sub(/^.*\[/, "", line)
      sub(/\].*$/, "", line)
      n = split(line, s, /[[:space:]]+/)
      if (n < 4) return ""
      return s[3] " " s[4]
    }
    function is_sample_time(line) {
      return (line ~ /time:/) && (line !~ /%/) && (line ~ /(ns|µs|us|ms| s)/)
    }
    {
      if ($0 ~ /^[a-z0-9][a-z0-9._/-]*$/ && $0 !~ /Benchmarking/ && $0 !~ /Analyzing/) {
        pending = $0
      } else if (is_sample_time($0)) {
        if ($2 ~ /^time:$/) {
          id = $1
        } else if ($1 ~ /^time:$/ && pending != "") {
          id = pending
        } else {
          next
        }
        pair = parse_triple($0)
        if (pair == "") next
        split(pair, t, /[[:space:]]+/)
        time_ns[id] = to_ns(t[1], t[2])
        if (!(id in thrpt)) thrpt[id] = "-"
        pending = ""
      } else if ($0 ~ /thrpt:/ && $0 !~ /%/ && id != "") {
        pair = parse_triple($0)
        if (pair != "") {
          split(pair, t, /[[:space:]]+/)
          thrpt[id] = to_elems(t[1], t[2])
        }
      }
    }
    END {
      for (k in time_ns) {
        printf "%s\t%d\t%s\n", k, time_ns[k], thrpt[k]
      }
    }
  ' "$log" | { printf 'benchmark_id\ttime_ns\tthroughput_per_s\n'; sort -t "$(printf '\t')" -k1,1; }
}

run_benches() {
  cargo bench -p zenith-float-num \
    --bench arithmetic --bench transcendentals --bench composite \
    --bench specials --bench linalg \
    -- --quick 2>&1
}

mode=compare
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage 0
elif [[ "${1:-}" == "--update" ]]; then
  mode=update
elif [[ "${1:-}" == "--update-from" ]]; then
  if [[ -z "${2:-}" ]]; then
    echo "error: --update-from requires a bench log file path" >&2
    exit 1
  fi
  cp "$2" "$tmp_log"
  mode=update-only
elif [[ "${1:-}" == "--check" ]]; then
  if [[ -z "${2:-}" ]]; then
    echo "error: --check requires a bench log file path" >&2
    exit 1
  fi
  cp "$2" "$tmp_log"
  mode=check-only
fi

if [[ "$mode" == "compare" || "$mode" == "update" ]]; then
  echo "Running Criterion benches (--quick) ..." >&2
  run_benches | tee "$tmp_log"
fi

parse_bench_log "$tmp_log" > "$tmp_tsv"
count="$(tail -n +2 "$tmp_tsv" | wc -l | tr -d ' ')"
echo "Parsed $count benchmark rows." >&2

if [[ "$mode" == "update" || "$mode" == "update-only" ]]; then
  cp "$tmp_tsv" "$baseline"
  echo "Updated $baseline"
  exit 0
fi

if [[ ! -f "$baseline" ]]; then
  echo "error: missing baseline $baseline (run: ./scripts/bench-compare.sh --update)" >&2
  exit 1
fi

awk -F '\t' -v pct="$regression_pct" '
  BEGIN {
    thresh = 1 + pct / 100
    thr_thresh = 1 - pct / 100
    fails = 0
    missing = 0
    checked = 0
  }
  FNR == NR && FNR > 1 {
    base_time[$1] = $2
    base_thr[$1] = $3
    next
  }
  FNR > 1 {
    if (!($1 in base_time)) {
      missing++
      next
    }
    checked++
    bt = base_time[$1] + 0
    nt = $2 + 0
    if (bt > 0 && nt > bt * thresh) {
      printf "REGRESSION time  %s: %d ns (baseline %d ns, +%.1f%%)\n", $1, nt, bt, (nt / bt - 1) * 100
      fails++
    }
    bs = base_thr[$1]
    ns = $3
    if (bs != "-" && ns != "-" && bs + 0 > 0 && ns + 0 < (bs + 0) * thr_thresh) {
      printf "REGRESSION thrpt %s: %s elem/s (baseline %s elem/s, %.1f%%)\n", $1, ns, bs, (ns / (bs + 0) - 1) * 100
      fails++
    }
  }
  END {
    if (checked == 0) {
      print "error: no overlapping benchmarks between run and baseline" > "/dev/stderr"
      exit 2
    }
    if (missing > 0) {
      printf "note: %d benchmarks in run are not in baseline (new benches? run --update)\n", missing > "/dev/stderr"
    }
    if (fails > 0) {
      printf "%d regression(s) over %g%% threshold\n", fails, pct > "/dev/stderr"
      exit 1
    }
    printf "ok: %d benchmarks within %g%% of baseline\n", checked, pct > "/dev/stderr"
  }
' "$baseline" "$tmp_tsv"
