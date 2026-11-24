//! Wall-clock budgets for the library's unit tests on emulated runners.
//!
//! A unit test that bounds a wall-clock time states the native budget. Under an emulator (qemu-aarch64
//! in `.github/workflows/arm-tests.yml` and `tools/full_gate.sh`) the same work takes `ZF_TIME_FACTOR`
//! times longer (pin `tools/qemu_aarch64_zf_time_factor.txt`, measured 6.25 + 0.25 margin), so the
//! budget is the native one times that factor: the same work allowance, not a looser one. Native runs
//! leave `ZF_TIME_FACTOR` unset (factor 1). Same grammar and range as `tests/s9_specials_golds.rs`.
//!
//! Release step only (Jeff 2026-11-23 23:04 CT: "timing limits are enforced in the release step only;
//! debug still runs every test for correctness, and the debug step limit goes to 60 min"): in a debug
//! build (`debug_assertions`) [`within`], [`within_ms`] and [`limit_holds`] pass, so a test that misses
//! its wall-clock limit still runs the rest of its body and every correctness assert. The budgets keep
//! their values (pinned in `tools/time_limits.tsv`). A hang watchdog that kills a child or stops waiting
//! (tier1 ZF-78) uses [`under_budget`], which applies in both builds.

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

/// Wall-clock limits fail a test only in a release build (Jeff 2026-11-23 23:04 CT).
pub(crate) const LIMITS_ENFORCED: bool = !cfg!(debug_assertions);

/// The result of a wall-clock comparison for this build: `ok` in a release build, `true` in a debug
/// build. Every timing assert that is not a [`within`] / [`within_ms`] budget goes through here.
pub(crate) fn limit_holds(ok: bool) -> bool {
    ok || !LIMITS_ENFORCED
}

/// `elapsed` is under this run's budget for `native` (see [`budget`]), in every build. For hang
/// watchdogs only; a timing assert uses [`within`]. An invalid `ZF_TIME_FACTOR` counts as over budget.
pub(crate) fn under_budget(elapsed: Duration, native: Duration) -> bool {
    budget(native).is_ok_and(|b| elapsed < b)
}

/// The timing assert: [`under_budget`] in a release build, `true` in a debug build ([`limit_holds`]).
/// Every wall-clock bound in the library's unit tests goes through here or [`within_ms`] (Jeff
/// 2026-11-22 ~11:58 PM CT: qemu limits scale with the time factor).
pub(crate) fn within(elapsed: Duration, native: Duration) -> bool {
    limit_holds(under_budget(elapsed, native))
}

/// [`within`] for a test that keeps its time as `Duration::as_millis()`.
pub(crate) fn within_ms(elapsed_ms: u128, native_ms: u64) -> bool {
    limit_holds(
        u64::try_from(elapsed_ms).is_ok_and(|ms| {
            under_budget(Duration::from_millis(ms), Duration::from_millis(native_ms))
        }),
    )
}

/// A debug build waits this many times a row's release limit before it calls the row hung: the factor
/// `tests/specials_oracle.rs` already gave its unoptimised rows, now shared (pinned in tools/time_limits.tsv).
pub(crate) const DEBUG_WATCHDOG_FACTOR: u32 = 10;

/// Why a row's result did not arrive.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum RowWait {
    /// Release step: the row missed its wall-clock limit.
    Limit,
    /// Debug build: no result after the hang watchdog (`limit` x [`DEBUG_WATCHDOG_FACTOR`]).
    Hang,
    /// The worker stopped without sending.
    Gone,
}

/// Waits for one row's result. Release step: the row fails at `limit`, as before ([`RowWait::Limit`]). Debug
/// build: the limit is not enforced (Jeff 2026-11-23 23:04 CT); the wait stops only at the hang watchdog
/// `limit` x [`DEBUG_WATCHDOG_FACTOR`] ([`RowWait::Hang`]), so a slow row still has its result checked.
pub(crate) fn wait_row<T>(
    rx: &std::sync::mpsc::Receiver<T>,
    limit: Duration,
) -> Result<T, RowWait> {
    use std::sync::mpsc::RecvTimeoutError;
    let (wait, missed) = if LIMITS_ENFORCED {
        (limit, RowWait::Limit)
    } else {
        (limit.saturating_mul(DEBUG_WATCHDOG_FACTOR), RowWait::Hang)
    };
    rx.recv_timeout(wait).map_err(|e| match e {
        RecvTimeoutError::Timeout => missed,
        RecvTimeoutError::Disconnected => RowWait::Gone,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A row that takes 3x its limit fails in the release step only; in a debug build its result arrives and
    /// is checked. A row that never answers is a hang in both (debug: at the watchdog). Spence R1 on 82e4afb.
    #[test]
    fn slow_row_fails_its_limit_in_release_only() {
        let limit = Duration::from_millis(20);
        let (tx, rx) = std::sync::mpsc::channel();
        let slow = std::thread::spawn(move || {
            std::thread::sleep(Duration::from_millis(60));
            let _sent = tx.send(7u8);
        });
        let got = wait_row(&rx, limit);
        let _joined = slow.join();
        if LIMITS_ENFORCED {
            assert_eq!(
                got,
                Err(RowWait::Limit),
                "release: a slow row must miss its limit"
            );
        } else {
            assert_eq!(
                got,
                Ok(7),
                "debug: a slow row is not a failure; its result must arrive"
            );
        }
        let (keep, rx) = std::sync::mpsc::channel::<u8>();
        let start = std::time::Instant::now();
        let hung = wait_row(&rx, Duration::from_millis(5));
        let waited = start.elapsed();
        drop(keep);
        if LIMITS_ENFORCED {
            assert_eq!(hung, Err(RowWait::Limit));
        } else {
            assert_eq!(
                hung,
                Err(RowWait::Hang),
                "debug: no answer is a hang at the watchdog"
            );
            assert!(
                waited >= Duration::from_millis(50),
                "debug watchdog is 10 x the limit, waited {waited:?}"
            );
        }
        let (gone, rx) = std::sync::mpsc::channel::<u8>();
        drop(gone);
        assert_eq!(wait_row(&rx, limit), Err(RowWait::Gone));
    }

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
        // Natively (factor 1) `d.as_secs() < n` and `ms < n` read exactly as before; in a debug
        // build the timing assert passes (release step only, Jeff 2026-11-23 23:04 CT).
        let one = Duration::from_secs(1);
        let release = LIMITS_ENFORCED;
        assert_eq!(release, !cfg!(debug_assertions));
        if std::env::var_os("ZF_TIME_FACTOR").is_none() {
            assert!(under_budget(Duration::from_millis(999), one));
            assert!(!under_budget(one, one));
            assert!(within(Duration::from_millis(999), one));
            assert_eq!(within(one, one), !release);
            assert!(within_ms(199, 200));
            assert_eq!(within_ms(200, 200), !release);
        }
        assert_eq!(within_ms(u128::MAX, 200), !release);
        assert_eq!(within(Duration::MAX, one), !release);
        assert!(!under_budget(Duration::MAX, one));
        assert!(limit_holds(true));
        assert_eq!(limit_holds(false), !release);
    }

    /// Jeff 2026-11-23 23:04 CT: "timing limits are enforced in the release step only; debug still
    /// runs every test for correctness". A run 10x over its budget passes the wall-clock check in a
    /// debug build (debug_assertions) and fails it in a release build, whatever ZF_TIME_FACTOR is
    /// (at most 8).
    #[test]
    fn debug_build_waives_wall_clock_limits_release_enforces_them() {
        let debug = cfg!(debug_assertions);
        assert_eq!(
            within(Duration::from_secs(10), Duration::from_secs(1)),
            debug
        );
        assert_eq!(within_ms(10_000, 1_000), debug);
        assert_eq!(within(Duration::MAX, Duration::from_secs(1)), debug);
    }

    /// A missed wall-clock limit must not end a debug test early: the test body after the timed
    /// section still runs, so its correctness asserts still fire. In a release build the timing
    /// assert itself fails first, as before.
    #[test]
    fn correctness_asserts_still_fire_after_a_missed_limit() {
        let ran = std::panic::catch_unwind(|| {
            assert!(
                within(Duration::from_secs(10), Duration::from_secs(1)),
                "TIMING limit missed"
            );
            assert_eq!(1 + 1, 3, "CORRECTNESS assert");
        });
        let payload: &str = match &ran {
            Ok(()) => "no panic",
            Err(e) => e
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| e.downcast_ref::<&str>().copied())
                .unwrap_or(""),
        };
        let want = if cfg!(debug_assertions) {
            "CORRECTNESS assert"
        } else {
            "TIMING limit missed"
        };
        assert!(
            payload.contains(want),
            "want a panic with {want:?}, got {payload:?}"
        );
    }
}
