// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
use zenith_float::ExactNum;
use zenith_float::Consts;
use zenith_float::RoundingMode;
use zenith_float_macro::expr;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rm = RoundingMode::None;
    let mut cc = Consts::new()?;
    let _res: ExactNum = expr!(-6 * atan(1.0 / sqrt(3)), (256, rm, &mut cc));
    Ok(())
}
