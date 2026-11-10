//! Tightness bookkeeping (S12, TEST_PLAN section 9 'Metric'): the ratio `rad / R_ref` of a case, exact integers only, and the results file
//! `tools/tightness_report.py` reads (format in that script's header). The checks themselves (the tight and loose assertions, the factors
//! `tau` and `lambda`) are decision D2 (Rusty) and stay with the callers: `check_tight` takes `tau` as an argument.
//! A ratio is `ceil(rad * 65536 / r_ref)` ("q16", 65536 = ratio 1), `rad` and `r_ref` being integer multiples of one common unit that the
//! caller aligned (never floats). Every failure is an `Err`; nothing here panics.

use std::collections::BTreeMap;
use std::path::Path;

/// Ratio 1 in q16.
pub const Q16_ONE: u64 = 65_536;

/// `ceil(rad * 65536 / r_ref)`; `Err` for `r_ref == 0` or a ratio that does not fit `u64`.
pub fn slack_q16(rad: u128, r_ref: u128) -> Result<u64, String> {
    if r_ref == 0 {
        return Err(String::from("slack_q16: R_ref is zero"));
    }
    let scaled = rad
        .checked_mul(u128::from(Q16_ONE))
        .ok_or_else(|| String::from("slack_q16: rad * 65536 overflows u128"))?;
    let q = scaled / r_ref + u128::from(scaled % r_ref != 0);
    u64::try_from(q).map_err(|_| String::from("slack_q16: ratio does not fit u64"))
}

/// The tight check of the plan: `rad <= tau * R_ref`, exactly. `Err` names the case; overflow is an `Err` too.
pub fn check_tight(rad: u128, r_ref: u128, tau: u64) -> Result<(), String> {
    let bound = r_ref
        .checked_mul(u128::from(tau))
        .ok_or_else(|| String::from("check_tight: tau * R_ref overflows u128"))?;
    if rad <= bound {
        Ok(())
    } else {
        Err(format!(
            "tight check failed: rad {} > tau {} * R_ref {} = {}",
            rad, tau, r_ref, bound
        ))
    }
}

/// The zero-based index of the p99 element among `n` sorted samples: `ceil(99 n / 100) - 1`. `Err` for `n == 0` (no sample) or when `99 n` overflows.
pub fn p99_rank(n: usize) -> Result<usize, String> {
    n.checked_mul(99)
        .ok_or_else(|| format!("p99_rank: 99 * {} overflows usize", n))?
        .div_ceil(100)
        .checked_sub(1)
        .ok_or_else(|| String::from("p99_rank: no samples"))
}

/// The statistics of one `(op, prec)` cell, all ratios in q16.
#[derive(Debug, PartialEq, Eq)]
pub struct CellStat {
    pub op: String,
    pub prec: u32,
    pub count: u64,
    pub min: u64,
    pub median: u64,
    pub p99: u64,
    pub max: u64,
}

/// The ratios recorded per `(op, prec)`.
pub struct TightnessLog {
    cells: BTreeMap<(String, u32), Vec<u64>>,
}

const HEADER: &str = "# zenith-bound tightness results (S12), written by the tightness sets; ratio = rad / R_ref in units of 1/65536 (q16), rounded up\n# columns: OP PREC COUNT MIN MEDIAN P99 MAX\n";

impl TightnessLog {
    pub fn empty() -> TightnessLog {
        TightnessLog {
            cells: BTreeMap::new(),
        }
    }

    /// Records one case of `op` at `prec` (an op name is lowercase letters and digits, `prec >= 1`).
    pub fn record(&mut self, op: &str, prec: u32, rad: u128, r_ref: u128) -> Result<(), String> {
        if op.is_empty()
            || !op
                .chars()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit())
        {
            return Err(format!("record: bad op name {:?}", op));
        }
        if prec == 0 {
            return Err(String::from("record: prec 0"));
        }
        let q = slack_q16(rad, r_ref)?;
        self.cells
            .entry((String::from(op), prec))
            .or_default()
            .push(q);
        Ok(())
    }

    /// The per-cell statistics, sorted by `(op, prec)`. Median is the lower median, p99 is the element at rank `ceil(0.99 n)`.
    pub fn stats(&self) -> Result<Vec<CellStat>, String> {
        let mut out = Vec::new();
        for ((op, prec), ratios) in self.cells.iter() {
            let mut sorted: Vec<u64> = Vec::new();
            for q in ratios {
                sorted.push(*q);
            }
            sorted.sort_unstable();
            let n = sorted.len();
            let at = |rank: usize| match sorted.get(rank) {
                Some(&q) => Ok(q),
                None => Err(format!("stats: rank {} outside {} samples", rank, n)),
            };
            let last = n
                .checked_sub(1)
                .ok_or_else(|| String::from("stats: empty cell"))?;
            let p99_at = p99_rank(n)?;
            out.push(CellStat {
                op: String::from(op.as_str()),
                prec: *prec,
                count: u64::try_from(n).map_err(|_| String::from("stats: count"))?,
                min: at(0)?,
                median: at(last / 2)?,
                p99: at(p99_at)?,
                max: at(last)?,
            });
        }
        Ok(out)
    }

    /// The results file text (header, then one row per cell).
    pub fn render(&self) -> Result<String, String> {
        let mut s = String::from(HEADER);
        for c in self.stats()? {
            s.push_str(&format!(
                "{} {} {} {} {} {} {}\n",
                c.op, c.prec, c.count, c.min, c.median, c.p99, c.max
            ));
        }
        Ok(s)
    }

    /// Writes `render()` to `path`.
    pub fn write_results(&self, path: &Path) -> Result<(), String> {
        std::fs::write(path, self.render()?).map_err(|e| format!("{:?}: {}", path, e))
    }
}

/// Parses a results text (the inverse of `render`): strict, sorted, no duplicates, `MIN <= MEDIAN <= P99 <= MAX`.
pub fn parse_results(text: &str) -> Result<Vec<CellStat>, String> {
    if text.contains('\r') || (!text.is_empty() && !text.ends_with('\n')) {
        return Err(String::from(
            "results need LF line endings and a final newline",
        ));
    }
    let mut out: Vec<CellStat> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split(' ').collect();
        if f.len() != 7 {
            return Err(format!("{:?}: want 7 fields, got {}", line, f.len()));
        }
        let get = |i: usize| f.get(i).ok_or_else(|| format!("{:?}: field {}", line, i));
        let n = |i: usize| -> Result<u64, String> {
            let tok = *get(i)?;
            let v: u64 = tok
                .parse()
                .map_err(|_| format!("{:?}: {:?} is not a number", line, tok))?;
            if v.to_string() != tok {
                return Err(format!("{:?}: {:?} is not canonical", line, tok));
            }
            Ok(v)
        };
        let op = String::from(*get(0)?);
        let prec = u32::try_from(n(1)?).map_err(|_| format!("{:?}: prec", line))?;
        let (count, min, median, p99, max) = (n(2)?, n(3)?, n(4)?, n(5)?, n(6)?);
        if prec == 0 || count == 0 || !(min <= median && median <= p99 && p99 <= max) {
            return Err(format!(
                "{:?}: prec and count must be positive and MIN <= MEDIAN <= P99 <= MAX",
                line
            ));
        }
        if let Some(last) = out.last() {
            if (last.op.as_str(), last.prec) >= (op.as_str(), prec) {
                return Err(format!(
                    "{:?}: rows must be strictly sorted by (OP, PREC)",
                    line
                ));
            }
        }
        out.push(CellStat {
            op,
            prec,
            count,
            min,
            median,
            p99,
            max,
        });
    }
    Ok(out)
}
