//! Deterministic seeded PRNG (splitmix64, no dependency) and the property-test configuration: case counts and seeds come from the environment so that fast fixed-seed runs go in the gate and wider runs are on demand.

/// splitmix64 generator (Steele, Lea, Flood), a deterministic 64-bit PRNG.
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Creates a new SplitMix64 PRNG with the given seed.
    pub fn new(seed: u64) -> SplitMix64 {
        SplitMix64 { state: seed }
    }

    /// Returns the next 64-bit random number.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Returns the next 64-bit random number as an i64.
    pub fn next_i64(&mut self) -> i64 {
        let bytes = self.next_u64().to_ne_bytes();
        i64::from_ne_bytes(bytes)
    }

    /// Returns the next random byte.
    pub fn byte(&mut self) -> u8 {
        let [b0, ..] = self.next_u64().to_le_bytes();
        b0
    }

    /// Returns a random number in range [0, n).
    pub fn below(&mut self, n: u8) -> u8 {
        if n == 0 {
            0
        } else {
            self.byte() % n
        }
    }
}

/// Base seed of every property test unless ZB_PROP_SEED overrides it (each test mixes in the hash of its own name).
pub const DEFAULT_SEED: u64 = 0x5EED_2B0D_C0FF_EE01;

/// Number of cases and base seed of one property run.
#[derive(Debug)]
pub struct PropConfig {
    pub seed: u64,
    pub cases: u64,
}

/// Reads property test configuration from environment variables.
pub fn config(default_cases: u64) -> Result<PropConfig, String> {
    let cases = match std::env::var("ZB_PROP_CASES") {
        Ok(text) => text
            .parse::<u64>()
            .map_err(|_| format!("ZB_PROP_CASES is not a number: {}", text))?,
        Err(std::env::VarError::NotPresent) => default_cases,
        Err(_) => return Err(String::from("ZB_PROP_CASES is not a number")),
    };

    let seed = match std::env::var("ZB_PROP_SEED") {
        Ok(text) => match text.strip_prefix("0x") {
            Some(hex) => u64::from_str_radix(hex, 16)
                .map_err(|_| format!("ZB_PROP_SEED is not a hex number: {}", text))?,
            None => text
                .parse::<u64>()
                .map_err(|_| format!("ZB_PROP_SEED is not a number: {}", text))?,
        },
        Err(std::env::VarError::NotPresent) => DEFAULT_SEED,
        Err(_) => return Err(String::from("ZB_PROP_SEED is not a number")),
    };

    Ok(PropConfig { seed, cases })
}

/// Computes the seed for a specific test case.
pub fn case_seed(base: u64, case: u64) -> u64 {
    if case == 0 {
        base
    } else {
        SplitMix64::new(base ^ case.wrapping_mul(0x9E37_79B9_7F4A_7C15)).next_u64()
    }
}

/// FNV-1a hash of a test name; it is mixed into the base seed so that each property test draws its own cases.
pub fn name_hash(name: &str) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325;
    for b in name.bytes() {
        h = (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// Returns a hint for replaying a failing test case: `ZB_PROP_SEED` is the case seed with the test name hash undone,
/// so that running only the test `name` with one case reproduces it.
pub fn replay_hint(name: &str, seed: u64) -> String {
    format!(
        "replay: ZB_PROP_SEED={:#x} ZB_PROP_CASES=1 cargo test {}",
        seed ^ name_hash(name),
        name
    )
}

/// Runs `body` once per case of the configured run; a failure names the case, its seed and the replay command.
/// The base seed is mixed with the hash of `name`, so two tests with the same `ZB_PROP_SEED` draw different cases.
pub fn run_cases(
    name: &str,
    default_cases: u64,
    body: &dyn Fn(&mut SplitMix64) -> Result<(), String>,
) -> Result<(), String> {
    let cfg = config(default_cases)?;
    let base = cfg.seed ^ name_hash(name);
    for case in 0..cfg.cases {
        let seed = case_seed(base, case);
        let mut rng = SplitMix64::new(seed);
        body(&mut rng).map_err(|e| {
            format!(
                "{}: case {} seed {:#x}: {}; {}",
                name,
                case,
                seed,
                e,
                replay_hint(name, seed)
            )
        })?;
    }
    Ok(())
}
