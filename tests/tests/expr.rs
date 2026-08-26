use zenith_float::ExactNum;
use zenith_float::Consts;
use zenith_float::RoundingMode;
use zenith_float_macro::expr;

fn main() {
    let rm = RoundingMode::None;
    let mut cc = Consts::new().unwrap();
    let _res: ExactNum = expr!(-6 * atan(1.0 / sqrt(3)), (256, rm, &mut cc));
}
