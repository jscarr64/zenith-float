//! Singular edge values that must be infinities rather than NaN.

use zenith_float_num::{Consts, ExactNum, RoundingMode};

const RM: RoundingMode = RoundingMode::ToEven;

#[test]
fn elliptic_f_at_one_one_is_infinite() {
    // F(±1 | 1) = ±∞ (1.0.4 returned NaN(InvalidArgument)); elliptic_k(1) = +∞ already.
    let mut cc = Consts::new().unwrap();
    for p in [64, 128, 512] {
        let one = ExactNum::from_u8(1, p);
        assert!(one.elliptic_f(&one, p, RM, &mut cc).is_inf_pos());
        assert!(one.neg().elliptic_f(&one, p, RM, &mut cc).is_inf_neg());
        assert!(one.elliptic_k(p, RM, &mut cc).is_inf_pos());
        // Just inside the domain stays finite.
        let m = ExactNum::parse("0.999", zenith_float_num::Radix::Dec, p, RM, &mut cc);
        assert!(!one.elliptic_f(&m, p, RM, &mut cc).is_inf());
    }
}
