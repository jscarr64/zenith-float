//! Oracle checkers: a round-up result must be the smallest `Mag` at or above the exact value and a round-down result the largest at or below it, decided by exact comparisons (cross-multiplication for quotients, squares for square roots).

use super::dyadic::{Dy, MagV};
use core::cmp::Ordering;

pub const TWO63: u64 = 1 << 63;

pub fn max_mag() -> MagV {
    MagV::Finite {
        exp: i64::MAX,
        mant: u64::MAX,
    }
}

pub fn min_pos_mag() -> MagV {
    MagV::Finite {
        exp: i64::MIN,
        mant: TWO63,
    }
}

pub fn zero_mag() -> MagV {
    MagV::Finite { exp: 0, mant: 0 }
}

pub fn is_canonical(m: &MagV) -> bool {
    match m {
        MagV::Infinity => true,
        MagV::Finite { exp, mant } => *mant >= TWO63 || (*exp == 0 && *mant == 0),
    }
}

pub fn succ(m: &MagV) -> MagV {
    match m {
        MagV::Infinity => MagV::Infinity,
        MagV::Finite { exp, mant } => {
            if *mant == 0 {
                min_pos_mag()
            } else if *mant == u64::MAX {
                if *exp == i64::MAX {
                    MagV::Infinity
                } else {
                    MagV::Finite {
                        exp: exp + 1,
                        mant: TWO63,
                    }
                }
            } else {
                MagV::Finite {
                    exp: *exp,
                    mant: *mant + 1,
                }
            }
        }
    }
}

pub fn pred(m: &MagV) -> Option<MagV> {
    match m {
        MagV::Infinity => None,
        MagV::Finite { exp, mant } => {
            if *mant == 0 {
                None
            } else if *mant == TWO63 {
                if *exp == i64::MIN {
                    None
                } else {
                    Some(MagV::Finite {
                        exp: exp - 1,
                        mant: u64::MAX,
                    })
                }
            } else {
                Some(MagV::Finite {
                    exp: *exp,
                    mant: *mant - 1,
                })
            }
        }
    }
}

pub fn cmp_dy(m: &MagV, exact: &Dy) -> Option<Ordering> {
    Dy::from_magv(m)?.cmp(exact)
}

pub fn cmp_quotient(m: &MagV, num: &Dy, den: &Dy) -> Option<Ordering> {
    let cand = Dy::from_magv(m)?;
    let prod = cand.mul(den);
    prod.cmp(num)
}

pub fn cmp_sqrt(m: &MagV, x: &Dy) -> Option<Ordering> {
    let cand = Dy::from_magv(m)?;
    let sq = cand.mul(&cand);
    sq.cmp(x)
}

pub fn cmp_minus(m: &MagV, a: &Dy, b: &Dy) -> Option<Ordering> {
    let cand = Dy::from_magv(m)?;
    let sum = cand.add(b)?;
    sum.cmp(a)
}

pub fn check_up(
    r: &MagV,
    exact_is_zero: bool,
    exact_cmp: &dyn Fn(&MagV) -> Option<Ordering>,
) -> Result<(), String> {
    if !is_canonical(r) {
        return Err(format!("non-canonical {:?}", r));
    }
    if exact_is_zero {
        if *r == zero_mag() {
            return Ok(());
        } else {
            return Err(format!("must be zero {:?}", r));
        }
    }
    match r {
        MagV::Infinity => match exact_cmp(&max_mag()) {
            Some(Ordering::Less) => Ok(()),
            Some(_) => Err(format!(
                "Infinity returned although a finite upper bound exists {:?}",
                r
            )),
            None => Err(format!("model cannot align {:?}", r)),
        },
        MagV::Finite { exp: _, mant } => {
            if *mant == 0 {
                return Err(format!("zero returned for a positive exact value {:?}", r));
            }
            match exact_cmp(r) {
                Some(Ordering::Less) => return Err(format!("not an upper bound {:?}", r)),
                None => return Err(format!("model cannot align {:?}", r)),
                Some(_) => {}
            }
            if let Some(p) = pred(r) {
                match exact_cmp(&p) {
                    Some(Ordering::Less) => Ok(()),
                    Some(_) => Err(format!("not the smallest upper bound {:?}", r)),
                    None => Err(format!("model cannot align {:?}", r)),
                }
            } else {
                Ok(())
            }
        }
    }
}

pub fn check_down(
    r: &MagV,
    exact_is_zero: bool,
    exact_cmp: &dyn Fn(&MagV) -> Option<Ordering>,
) -> Result<(), String> {
    if !is_canonical(r) {
        return Err(format!("non-canonical {:?}", r));
    }
    if exact_is_zero {
        if *r == zero_mag() {
            return Ok(());
        } else {
            return Err(format!("must be zero {:?}", r));
        }
    }
    match r {
        MagV::Infinity => Err(format!("Infinity is not a lower bound {:?}", r)),
        MagV::Finite { exp: _, mant } => {
            if *mant == 0 {
                match exact_cmp(&min_pos_mag()) {
                    Some(Ordering::Greater) => return Ok(()),
                    Some(_) => {
                        return Err(format!(
                            "zero returned for a representable positive value {:?}",
                            r
                        ))
                    }
                    None => return Err(format!("model cannot align {:?}", r)),
                }
            }
            match exact_cmp(r) {
                Some(Ordering::Greater) => return Err(format!("not a lower bound {:?}", r)),
                None => return Err(format!("model cannot align {:?}", r)),
                Some(_) => {}
            }
            let s = succ(r);
            match s {
                MagV::Infinity => Ok(()),
                MagV::Finite { .. } => match exact_cmp(&s) {
                    Some(Ordering::Greater) => Ok(()),
                    Some(_) => Err(format!("not the largest lower bound {:?}", r)),
                    None => Err(format!("model cannot align {:?}", r)),
                },
            }
        }
    }
}
