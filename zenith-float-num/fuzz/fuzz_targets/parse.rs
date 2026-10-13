#![no_main]

use libfuzzer_sys::fuzz_target;
use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

fuzz_target!(|data: &[u8]| {
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };
    let Ok(mut consts) = Consts::new() else {
        return;
    };
    for mode in [
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
        RoundingMode::None,
    ] {
        let n = ExactNum::parse(text, Radix::Dec, 53, mode, &mut consts);
        let _ = (n.is_nan(), n.is_inf(), n.is_zero(), n.sign(), n.err());
    }
});
