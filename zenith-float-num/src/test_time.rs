//! Wall-clock budgets for the library's unit tests on emulated runners.
//!
//! A unit test that bounds a wall-clock time states the native budget. Under an emulator (qemu-aarch64
//! in `.github/workflows/arm-tests.yml` and `tools/full_gate.sh`) the same work takes `ZF_TIME_FACTOR`
//! times longer (pin `tools/qemu_aarch64_zf_time_factor.txt`, measured 6.25 + 0.25 margin), so the
//! budget is the native one times that factor: the same work allowance, not a looser one. Native runs
//! leave `ZF_TIME_FACTOR` unset (factor 1). Same grammar and range as `tests/s9_specials_golds.rs`.

use core::time::Duration;

/// Largest accepted factor, in thousandths (8.000), as in `tests/s9_specials_golds.rs`.
const MAX_TIME_FACTOR_MILLI: u64 = 8000;

/// `ZF_TIME_FACTOR` in thousandths: unset is 1000; otherwise ASCII digits with at most 3 fraction
/// digits and a value in [1, 8]. Anything else is `Err`, so a typo fails the test instead of
/// silently changing the budget.
pub(crate) fn time_factor_milli(raw: Option<&str>) -> Result<u64, String> {
    let Some(s) = raw else { return Ok(1000) };
    let bad =
        || format!("ZF_TIME_FACTOR={s:?}: want a decimal in [1, 8] with at most 3 fraction digits");
    let (int, frac, dot) = match s.split_once('.') {
        Some((i, f)) => (i, f, true),
        None => (s, "", false),
    };
    let digits = |t: &str| !t.is_empty() && t.bytes().all(|b| b.is_ascii_digit());
    if !digits(int) || (dot && !digits(frac)) || frac.len() > 3 {
        return Err(bad());
    }
    let whole: u64 = int.parse().map_err(|_| bad())?;
    let mut part: u64 = if dot { frac.parse().map_err(|_| bad())? } else { 0 };
    for _ in frac.len()..3 {
        part = part.checked_mul(10).ok_or_else(bad)?;
    }
    let milli = whole
        .checked_mul(1000)
        .and_then(|m| m.checked_add(part))
        .ok_or_else(bad)?;
    if !(1000..=MAX_TIME_FACTOR_MILLI).contains(&milli) {
        return Err(bad());
    }
    Ok(milli)
}

/// `native` times the factor (exact: native ms times the milli-factor is µs).
pub(crate) fn scaled(native: Duration, milli: u64) -> Result<Duration, String> {
    let ms = u64::try_from(native.as_millis()).map_err(|e| format!("budget {native:?}: {e}"))?;
    ms.checked_mul(milli)
        .map(Duration::from_micros)
        .ok_or_else(|| format!("budget {native:?} x {milli}/1000 overflows"))
}

/// The budget for this run: `native` times `ZF_TIME_FACTOR` (unset: `native`).
pub(crate) fn budget(native: Duration) -> Result<Duration, String> {
    let milli = match std::env::var("ZF_TIME_FACTOR") {
        Err(std::env::VarError::NotPresent) => time_factor_milli(None)?,
        Err(e) => return Err(format!("ZF_TIME_FACTOR: {e}")),
        Ok(s) => time_factor_milli(Some(&s))?,
    };
    scaled(native, milli)
}

/// `elapsed` is under this run's budget for `native` (see [`budget`]). Every wall-clock bound in the
/// library's unit tests goes through here or [`within_ms`] (Jeff 2026-11-22 ~11:58 PM CT: qemu limits
/// scale with the time factor). An invalid `ZF_TIME_FACTOR` counts as over budget, so it fails the test.
pub(crate) fn within(elapsed: Duration, native: Duration) -> bool {
    budget(native).is_ok_and(|b| elapsed < b)
}

/// [`within`] for a test that keeps its time as `Duration::as_millis()`.
pub(crate) fn within_ms(elapsed_ms: u128, native_ms: u64) -> bool {
    u64::try_from(elapsed_ms)
        .is_ok_and(|ms| within(Duration::from_millis(ms), Duration::from_millis(native_ms)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_grammar_and_scaling() -> Result<(), String> {
        assert_eq!(time_factor_milli(None)?, 1000);
        assert_eq!(time_factor_milli(Some("6.50"))?, 6500);
        assert_eq!(time_factor_milli(Some("8"))?, 8000);
        for bad in ["", "0.5", "8.001", "6.5e0", "-2", "inf", "6.5000", "x"] {
            assert!(time_factor_milli(Some(bad)).is_err(), "{bad:?} accepted");
        }
        // qemu pin 6.50: the 30 s ln_gamma and 1 s APRCL budgets become 195 s and 6.5 s.
        assert_eq!(
            scaled(Duration::from_secs(30), 6500)?,
            Duration::from_secs(195)
        );
        assert_eq!(
            scaled(Duration::from_secs(1), 6500)?,
            Duration::from_millis(6500)
        );
        assert_eq!(
            scaled(Duration::from_secs(1), 1000)?,
            Duration::from_secs(1)
        );
        Ok(())
    }

    #[test]
    fn within_is_strict_and_matches_the_old_native_comparisons() {
        // Natively (factor 1) `d.as_secs() < n` and `ms < n` read exactly as before.
        let one = Duration::from_secs(1);
        if std::env::var_os("ZF_TIME_FACTOR").is_none() {
            assert!(within(Duration::from_millis(999), one));
            assert!(!within(one, one));
            assert!(within_ms(199, 200));
            assert!(!within_ms(200, 200));
        }
        assert!(!within_ms(u128::MAX, 200));
        assert!(!within(Duration::MAX, one));
    }
}
