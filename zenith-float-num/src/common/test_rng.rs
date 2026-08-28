//! Seeded RNG for tests and for [`crate::ExactNum::random_normal`].
//!
//! Unit tests always use a deterministic stream. Integration tests and the
//! `random` feature use OS entropy unless [`reseed_random`] is called or
//! `ZENITH_TEST_SEED` is set. Replay a failure with
//! `ZENITH_TEST_SEED=<printed value> cargo test <name> -- --test-threads=1`.

use core::cell::RefCell;
use core::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Default seed when `ZENITH_TEST_SEED` is unset.
pub const DEFAULT_RANDOM_SEED: u64 = 0x5EED_CAFE_BADC_0D00;

static SEED: AtomicU64 = AtomicU64::new(DEFAULT_RANDOM_SEED);
static EPOCH: AtomicU64 = AtomicU64::new(1);
static SEEDED: AtomicBool = AtomicBool::new(false);
static INIT: AtomicBool = AtomicBool::new(false);

thread_local! {
    static RNG: RefCell<(u64, StdRng)> = RefCell::new((
        0,
        StdRng::seed_from_u64(DEFAULT_RANDOM_SEED),
    ));
}

fn env_seed() -> Option<u64> {
    #[cfg(feature = "std")]
    {
        std::env::var("ZENITH_TEST_SEED")
            .ok()
            .and_then(|s| s.parse().ok())
    }
    #[cfg(not(feature = "std"))]
    {
        None
    }
}

fn ensure_init() {
    if INIT.load(Ordering::Acquire) {
        return;
    }
    if cfg!(test) {
        reseed_random(env_seed().unwrap_or(DEFAULT_RANDOM_SEED));
        return;
    }
    if let Some(s) = env_seed() {
        reseed_random(s);
        return;
    }
    INIT.store(true, Ordering::Release);
}

fn install_panic_hook() {
    #[cfg(all(test, feature = "std"))]
    {
        use std::sync::Once;
        static HOOK: Once = Once::new();
        HOOK.call_once(|| {
            let prev = std::panic::take_hook();
            std::panic::set_hook(Box::new(move |info| {
                let s = random_seed();
                eprintln!("zenith-float test RNG seed: {s} (replay: ZENITH_TEST_SEED={s})");
                prev(info);
            }));
        });
    }
}

/// Seed currently used by the deterministic generator.
pub fn random_seed() -> u64 {
    ensure_init();
    // Keep `random` live in every cfg that compiles this module (lib, tests, rust-analyzer).
    let _keep: fn() -> u32 = random;
    let _ = _keep;
    SEED.load(Ordering::Relaxed)
}

/// Reseed the generator used by tests and `random_normal`.
pub fn reseed_random(seed: u64) {
    SEED.store(seed, Ordering::Relaxed);
    SEEDED.store(true, Ordering::Relaxed);
    EPOCH.fetch_add(1, Ordering::Relaxed);
    INIT.store(true, Ordering::Relaxed);
    install_panic_hook();
}

fn draw<T>() -> T
where
    rand::distributions::Standard: rand::distributions::Distribution<T>,
{
    ensure_init();
    if cfg!(test) || SEEDED.load(Ordering::Relaxed) {
        let epoch = EPOCH.load(Ordering::Relaxed);
        let seed = SEED.load(Ordering::Relaxed);
        RNG.with(|cell| {
            let mut g = cell.borrow_mut();
            if g.0 != epoch {
                *g = (epoch, StdRng::seed_from_u64(seed));
            }
            g.1.gen()
        })
    } else {
        rand::random()
    }
}

/// Draw a random value from the crate RNG used by tests and `random_normal`.
///
/// After [`reseed_random`], or when `ZENITH_TEST_SEED` is set, the stream is deterministic.
/// Unit tests default to [`DEFAULT_RANDOM_SEED`]. Replay with
/// `ZENITH_TEST_SEED=<n> cargo test <name> -- --test-threads=1`.
pub(crate) fn random<T>() -> T
where
    rand::distributions::Standard: rand::distributions::Distribution<T>,
{
    draw()
}

/// Public name for [`random`] when the `random` feature is enabled.
#[cfg(feature = "random")]
pub fn seeded_random<T>() -> T
where
    rand::distributions::Standard: rand::distributions::Distribution<T>,
{
    random()
}

#[cfg(test)]
#[test]
fn test_rng_draws() {
    let _ = random::<u32>();
    let _ = random::<u64>();
}
