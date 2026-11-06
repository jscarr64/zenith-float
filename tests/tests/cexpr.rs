use zenith_float::Consts;
use zenith_float::RoundingMode;
use zenith_float::cexpr;
use zenith_float::ExactComplex;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rm = RoundingMode::None;
    let mut cc = Consts::new()?;
    let _z: ExactComplex = cexpr!(I * I + exp(0), (256, rm, &mut cc));
    Ok(())
}
