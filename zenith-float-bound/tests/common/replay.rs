//! S2/S4 failure report and replay (TEST_PLAN section 7). A failing property returns an `Err` whose text is ONE replay line:
//!
//!   `suite=<name> seed=<hex> case=<n> op=<name> prec=<p> rm=<mode> a=<enc> b=<enc>`
//!
//! `seed` is the BASE seed of the run (16 hex digits) and `case` the case index; the case stream is `case_seed(seed, case)`.
//! Environment (read through `std::env::var`, no unsafe; a malformed value is an `Err`, never a panic or a silent default):
//!   ZB_FUZZ_SEED=<hex, optional 0x>   base seed                         (fallback: ZB_PROP_SEED, the builder's name, `0x...` hex or decimal)
//!   ZB_FUZZ_CASES=<count>             number of cases, at least 1       (fallback: ZB_PROP_CASES)
//!   ZB_REPLAY=<seed>:<case>           run exactly that case: seed hex (optional 0x), case decimal; overrides seed and cases
//! The ZB_FUZZ_* names win when both sets are present. Counts are counts, never wall-clock time (PLAN Part J).

/// The seed of case `case` of a run with base seed `base`. Mirrors the builder's `case_seed` in tests/model/rng.rs (case 0 is the base itself).
pub fn case_seed(base: u64, case: u64) -> u64 {
    if case == 0 {
        base
    } else {
        super::rng::splitmix64(base ^ super::modular::mul_u64(case, 0x9E37_79B9_7F4A_7C15))
    }
}

/// The fields of one replay line. Text fields contain no whitespace.
#[derive(Debug, PartialEq, Eq)]
pub struct ReplayLine {
    pub suite: String,
    pub seed: u64,
    pub case: u64,
    pub op: String,
    pub prec: usize,
    pub rm: String,
    pub a: String,
    pub b: String,
}

fn no_space(label: &str, text: &str) -> Result<(), String> {
    if text.is_empty() || text.chars().any(char::is_whitespace) {
        Err(format!(
            "{} must be non-empty and contain no whitespace: {:?}",
            label, text
        ))
    } else {
        Ok(())
    }
}

impl ReplayLine {
    /// The replay line. `Err` when a text field is empty or contains whitespace (it would not parse back).
    pub fn format_line(&self) -> Result<String, String> {
        no_space("suite", &self.suite)?;
        no_space("op", &self.op)?;
        no_space("rm", &self.rm)?;
        no_space("a", &self.a)?;
        no_space("b", &self.b)?;
        Ok(format!(
            "suite={} seed={:016x} case={} op={} prec={} rm={} a={} b={}",
            self.suite, self.seed, self.case, self.op, self.prec, self.rm, self.a, self.b
        ))
    }

    /// Parses a replay line: exactly the eight keys, in this order, each once, nothing else.
    pub fn parse_line(text: &str) -> Result<ReplayLine, String> {
        const KEYS: [&str; 8] = ["suite", "seed", "case", "op", "prec", "rm", "a", "b"];
        let mut values: Vec<&str> = Vec::new();
        let mut tokens = text.split_whitespace();
        for key in KEYS.iter() {
            let token = tokens
                .next()
                .ok_or_else(|| format!("replay line is missing {}=", key))?;
            let (k, v) = token
                .split_once('=')
                .ok_or_else(|| format!("replay token {:?} has no '='", token))?;
            if k != *key {
                return Err(format!(
                    "replay line has {:?} where {}= is expected",
                    k, key
                ));
            }
            if v.is_empty() {
                return Err(format!("replay line has an empty {}=", key));
            }
            values.push(v);
        }
        if let Some(extra) = tokens.next() {
            return Err(format!("replay line has an extra token {:?}", extra));
        }
        let mut it = values.into_iter();
        let mut next = || {
            it.next()
                .ok_or_else(|| String::from("replay line is short"))
        };
        let suite = String::from(next()?);
        let seed = parse_hex_u64("seed", next()?)?;
        let case = next()?
            .parse::<u64>()
            .map_err(|_| String::from("case is not a decimal number"))?;
        let op = String::from(next()?);
        let prec = next()?
            .parse::<usize>()
            .map_err(|_| String::from("prec is not a decimal number"))
            .and_then(|p| match p <= zenith_bound::MAX_PRECISION {
                true => Ok(p),
                false => Err(String::from("prec is above MAX_PRECISION")),
            })?;
        let rm = String::from(next()?);
        let a = String::from(next()?);
        let b = String::from(next()?);
        Ok(ReplayLine {
            suite,
            seed,
            case,
            op,
            prec,
            rm,
            a,
            b,
        })
    }
}

/// A hexadecimal `u64`, with or without a `0x` prefix.
pub fn parse_hex_u64(label: &str, text: &str) -> Result<u64, String> {
    let digits = text.strip_prefix("0x").unwrap_or(text);
    if digits.is_empty() {
        return Err(format!("{} is empty", label));
    }
    u64::from_str_radix(digits, 16)
        .map_err(|_| format!("{} is not a hex number: {:?}", label, text))
}

/// The case a replay asks for.
#[derive(Debug, PartialEq, Eq)]
pub struct ReplayTarget {
    pub seed: u64,
    pub case: u64,
}

/// The configuration of one run.
#[derive(Debug, PartialEq, Eq)]
pub struct FuzzConfig {
    pub seed: u64,
    pub cases: u64,
    pub replay: Option<ReplayTarget>,
}

impl FuzzConfig {
    /// The case indices to run: just the replayed case, or `0..cases`.
    pub fn case_range(&self) -> Result<core::ops::Range<u64>, String> {
        match &self.replay {
            Some(t) => {
                let end = t
                    .case
                    .checked_add(1)
                    .ok_or_else(|| String::from("replay case overflows"))?;
                Ok(t.case..end)
            }
            None => Ok(0..self.cases),
        }
    }

    /// The seed of the stream of case `case`: from the replay seed when replaying, else from the run seed.
    pub fn seed_of(&self, case: u64) -> u64 {
        match &self.replay {
            Some(t) => case_seed(t.seed, case),
            None => case_seed(self.seed, case),
        }
    }
}

/// Where a variable is read from: `Ok(None)` when it is not set, `Err` when it is set but not valid unicode.
pub type Lookup<'a> = dyn Fn(&str) -> Result<Option<String>, String> + 'a;

/// Parses the configuration from `lookup` (see the module header). `default_seed` and `default_cases` apply when nothing is set.
pub fn parse_config(
    lookup: &Lookup,
    default_seed: u64,
    default_cases: u64,
) -> Result<FuzzConfig, String> {
    let seed = match lookup("ZB_FUZZ_SEED")? {
        Some(text) => parse_hex_u64("ZB_FUZZ_SEED", &text)?,
        None => match lookup("ZB_PROP_SEED")? {
            Some(text) => match text.strip_prefix("0x") {
                Some(hex) => parse_hex_u64("ZB_PROP_SEED", hex)?,
                None => text
                    .parse::<u64>()
                    .map_err(|_| format!("ZB_PROP_SEED is not a number: {:?}", text))?,
            },
            None => default_seed,
        },
    };
    let cases_text = match lookup("ZB_FUZZ_CASES")? {
        Some(text) => Some(("ZB_FUZZ_CASES", text)),
        None => lookup("ZB_PROP_CASES")?.map(|text| ("ZB_PROP_CASES", text)),
    };
    let cases = match cases_text {
        Some((name, text)) => {
            let n = text
                .parse::<u64>()
                .map_err(|_| format!("{} is not a number: {:?}", name, text))?;
            if n == 0 {
                return Err(format!("{} must be at least 1", name));
            }
            n
        }
        None => default_cases,
    };
    let replay = match lookup("ZB_REPLAY")? {
        Some(text) => {
            let (seed_text, case_text) = text
                .split_once(':')
                .ok_or_else(|| format!("ZB_REPLAY must be <seed>:<case>, got {:?}", text))?;
            let seed = parse_hex_u64("ZB_REPLAY seed", seed_text)?;
            let case = case_text
                .parse::<u64>()
                .map_err(|_| format!("ZB_REPLAY case is not a decimal number: {:?}", case_text))?;
            Some(ReplayTarget { seed, case })
        }
        None => None,
    };
    Ok(FuzzConfig {
        seed,
        cases,
        replay,
    })
}

/// `parse_config` over the process environment (`std::env::var`).
pub fn config_from_env(default_seed: u64, default_cases: u64) -> Result<FuzzConfig, String> {
    let lookup = |name: &str| -> Result<Option<String>, String> {
        match std::env::var(name) {
            Ok(text) => Ok(Some(text)),
            Err(std::env::VarError::NotPresent) => Ok(None),
            Err(std::env::VarError::NotUnicode(_)) => Err(format!("{} is not valid unicode", name)),
        }
    };
    parse_config(&lookup, default_seed, default_cases)
}
