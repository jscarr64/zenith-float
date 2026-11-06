//! The operand streams of one test thread must not move when another test thread reseeds.
//! `rebuild_from_seed` replays an operand from its seed; under a multi-threaded run a
//! concurrent `reset_test_rng` / `rebuild_from_seed` on another thread used to restart this
//! thread's crate stream from the other thread's seed, so the rebuilt operand differed from
//! the first one (the cause of the parallel-run failures of `rem_is_exact`,
//! `mpfr_compare_ops(_large)` and `subnormal_band_every_mode`).

use crate::mpfr::common::{rebuild_from_seed, reset_test_rng};
use std::sync::atomic::{AtomicBool, Ordering};
use zenith_float_num::ExactNum;

#[test]
fn rebuild_is_stable_under_concurrent_reseed() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let seed = 0x0123_4567_89AB_CDEF;
    let want = rebuild_from_seed(seed, || ExactNum::random_normal(256, -1000, 1000));
    let stop = AtomicBool::new(false);
    let mut mismatches = 0usize;
    std::thread::scope(|s| {
        s.spawn(|| {
            let mut i = 0u64;
            while !stop.load(Ordering::Relaxed) {
                zenith_float_num::reseed_random(i);
                i = i.wrapping_add(1);
            }
        });
        for _ in 0..20_000 {
            if rebuild_from_seed(seed, || ExactNum::random_normal(256, -1000, 1000)) != want {
                mismatches += 1;
            }
        }
        stop.store(true, Ordering::Relaxed);
    });
    assert_eq!(
        mismatches, 0,
        "rebuilt operand changed under a concurrent reseed"
    );
    Ok(())
}
