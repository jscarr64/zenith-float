//! Which error kind a result carries when both operands are NaN with different kinds.
//!
//! The kinds are `A` = `InvalidArgument`, `B` = `DivisionByZero`, `N` = no kind; a case name like `AB` is
//! (first operand, second operand). For `add`, `sub`, the six `+`/`-` operator forms, `add_ref`/`sub_ref`
//! and the complex versions the FIRST operand's kind wins, as it did before Rule R (the `&a + b` and
//! `&a - b` operators are implemented as `b + a` and keep that rule explicitly). Every non-`asinh`
//! expectation below was recorded on the pre-Rule-R code (b7606ed) and is unchanged.
//!
//! One public function does surface a different kind: `ExactComplex::asinh` adds `rad + z` where it used
//! to add `z + rad` (`z` is borrowed there, so the owned `rad` has to be the first operand), so with NaN
//! parts of different kinds its imaginary part now reports the kind of `rad` (which is the kind of the
//! real part of `z`) instead of the kind of `z.im`. That is pinned here as the new behaviour. Every
//! other reordered library site adds a constant (`n`, `1`, `pi`, ...) or a value derived from the other
//! operand, so it cannot meet two NaNs of different kinds; `ln_gamma` (which contains one) is pinned
//! below as unchanged.

use zenith_float_num::{Consts, Error, ExactComplex, ExactNum, RoundingMode};

const P: usize = 128;
const RM: RoundingMode = RoundingMode::ToEven;

fn nan_a() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
}
fn nan_b() -> ExactNum {
    ExactNum::nan(Some(Error::DivisionByZero))
}
fn nan_n() -> ExactNum {
    ExactNum::nan(None)
}
fn kind(x: &ExactNum) -> String {
    if x.is_nan() {
        format!("{:?}", x.err())
    } else {
        "not-nan".to_string()
    }
}
fn ckind(z: &ExactComplex) -> String {
    format!("re={} im={}", kind(z.re()), kind(z.im()))
}
fn cx(re: ExactNum, im: ExactNum) -> ExactComplex {
    ExactComplex::new(re, im)
}

fn observations() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut put = |k: String, v: String| out.push((k, v));
    let kinds: [(&str, fn() -> ExactNum); 3] = [("A", nan_a), ("B", nan_b), ("N", nan_n)];
    for (na, fa) in kinds {
        for (nb, fb) in kinds {
            let t = format!("{na}{nb}");
            put(format!("op &a+b {t}"), kind(&(&fa() + fb())));
            put(format!("op a+&b {t}"), kind(&(fa() + &fb())));
            put(format!("op a+b {t}"), kind(&(fa() + fb())));
            put(format!("op &a-b {t}"), kind(&(&fa() - fb())));
            put(format!("op a-&b {t}"), kind(&(fa() - &fb())));
            put(format!("op a-b {t}"), kind(&(fa() - fb())));
            put(format!("add {t}"), kind(&fa().add(&fb(), P, RM)));
            put(format!("sub {t}"), kind(&fa().sub(&fb(), P, RM)));
            put(
                format!("add_ref {t}"),
                kind(&ExactNum::add_ref(&fa(), &fb(), P, RM)),
            );
            put(
                format!("sub_ref {t}"),
                kind(&ExactNum::sub_ref(&fa(), &fb(), P, RM)),
            );
            let z1 = || cx(fa(), fb());
            let z2 = || cx(fb(), fa());
            put(format!("cadd {t}"), ckind(&z1().add(&z2(), P, RM)));
            put(format!("csub {t}"), ckind(&z1().sub(&z2(), P, RM)));
            put(
                format!("cadd_ref {t}"),
                ckind(&ExactComplex::add_ref(&z1(), &z2(), P, RM)),
            );
            put(
                format!("csub_ref {t}"),
                ckind(&ExactComplex::sub_ref(&z1(), &z2(), P, RM)),
            );
        }
    }
    // NaN against an ordinary value and against Inf: the NaN's kind, in either order.
    for (na, fa) in kinds {
        put(format!("val+{na}"), kind(&(ExactNum::from_u8(3, P) + fa())));
        put(format!("{na}+val"), kind(&(fa() + ExactNum::from_u8(3, P))));
        put(
            format!("&val+{na}"),
            kind(&(&ExactNum::from_u8(3, P) + fa())),
        );
        put(
            format!("&val-{na}"),
            kind(&(&ExactNum::from_u8(3, P) - fa())),
        );
        put(
            format!("inf+{na}"),
            kind(&(zenith_float_num::INF_POS + fa())),
        );
        put(
            format!("{na}-inf"),
            kind(&(fa() - zenith_float_num::INF_POS)),
        );
    }
    // Public functions whose internals were reordered, fed NaN parts of different kinds.
    let mut cc = Consts::new().unwrap();
    for (na, fa) in kinds {
        for (nb, fb) in kinds {
            let t = format!("{na}{nb}");
            let z = || cx(fa(), fb());
            put(format!("asinh {t}"), ckind(&z().asinh(P, RM, &mut cc)));
            put(
                format!("ln_gamma {t}"),
                ckind(&z().ln_gamma(P, RM, &mut cc)),
            );
        }
    }
    out
}

const EXPECTED: &str = "\
op &a+b AA => Some(InvalidArgument)
op a+&b AA => Some(InvalidArgument)
op a+b AA => Some(InvalidArgument)
op &a-b AA => Some(InvalidArgument)
op a-&b AA => Some(InvalidArgument)
op a-b AA => Some(InvalidArgument)
add AA => Some(InvalidArgument)
sub AA => Some(InvalidArgument)
add_ref AA => Some(InvalidArgument)
sub_ref AA => Some(InvalidArgument)
cadd AA => re=Some(InvalidArgument) im=Some(InvalidArgument)
csub AA => re=Some(InvalidArgument) im=Some(InvalidArgument)
cadd_ref AA => re=Some(InvalidArgument) im=Some(InvalidArgument)
csub_ref AA => re=Some(InvalidArgument) im=Some(InvalidArgument)
op &a+b AB => Some(InvalidArgument)
op a+&b AB => Some(InvalidArgument)
op a+b AB => Some(InvalidArgument)
op &a-b AB => Some(InvalidArgument)
op a-&b AB => Some(InvalidArgument)
op a-b AB => Some(InvalidArgument)
add AB => Some(InvalidArgument)
sub AB => Some(InvalidArgument)
add_ref AB => Some(InvalidArgument)
sub_ref AB => Some(InvalidArgument)
cadd AB => re=Some(InvalidArgument) im=Some(DivisionByZero)
csub AB => re=Some(InvalidArgument) im=Some(DivisionByZero)
cadd_ref AB => re=Some(InvalidArgument) im=Some(DivisionByZero)
csub_ref AB => re=Some(InvalidArgument) im=Some(DivisionByZero)
op &a+b AN => Some(InvalidArgument)
op a+&b AN => Some(InvalidArgument)
op a+b AN => Some(InvalidArgument)
op &a-b AN => Some(InvalidArgument)
op a-&b AN => Some(InvalidArgument)
op a-b AN => Some(InvalidArgument)
add AN => Some(InvalidArgument)
sub AN => Some(InvalidArgument)
add_ref AN => Some(InvalidArgument)
sub_ref AN => Some(InvalidArgument)
cadd AN => re=Some(InvalidArgument) im=None
csub AN => re=Some(InvalidArgument) im=None
cadd_ref AN => re=Some(InvalidArgument) im=None
csub_ref AN => re=Some(InvalidArgument) im=None
op &a+b BA => Some(DivisionByZero)
op a+&b BA => Some(DivisionByZero)
op a+b BA => Some(DivisionByZero)
op &a-b BA => Some(DivisionByZero)
op a-&b BA => Some(DivisionByZero)
op a-b BA => Some(DivisionByZero)
add BA => Some(DivisionByZero)
sub BA => Some(DivisionByZero)
add_ref BA => Some(DivisionByZero)
sub_ref BA => Some(DivisionByZero)
cadd BA => re=Some(DivisionByZero) im=Some(InvalidArgument)
csub BA => re=Some(DivisionByZero) im=Some(InvalidArgument)
cadd_ref BA => re=Some(DivisionByZero) im=Some(InvalidArgument)
csub_ref BA => re=Some(DivisionByZero) im=Some(InvalidArgument)
op &a+b BB => Some(DivisionByZero)
op a+&b BB => Some(DivisionByZero)
op a+b BB => Some(DivisionByZero)
op &a-b BB => Some(DivisionByZero)
op a-&b BB => Some(DivisionByZero)
op a-b BB => Some(DivisionByZero)
add BB => Some(DivisionByZero)
sub BB => Some(DivisionByZero)
add_ref BB => Some(DivisionByZero)
sub_ref BB => Some(DivisionByZero)
cadd BB => re=Some(DivisionByZero) im=Some(DivisionByZero)
csub BB => re=Some(DivisionByZero) im=Some(DivisionByZero)
cadd_ref BB => re=Some(DivisionByZero) im=Some(DivisionByZero)
csub_ref BB => re=Some(DivisionByZero) im=Some(DivisionByZero)
op &a+b BN => Some(DivisionByZero)
op a+&b BN => Some(DivisionByZero)
op a+b BN => Some(DivisionByZero)
op &a-b BN => Some(DivisionByZero)
op a-&b BN => Some(DivisionByZero)
op a-b BN => Some(DivisionByZero)
add BN => Some(DivisionByZero)
sub BN => Some(DivisionByZero)
add_ref BN => Some(DivisionByZero)
sub_ref BN => Some(DivisionByZero)
cadd BN => re=Some(DivisionByZero) im=None
csub BN => re=Some(DivisionByZero) im=None
cadd_ref BN => re=Some(DivisionByZero) im=None
csub_ref BN => re=Some(DivisionByZero) im=None
op &a+b NA => None
op a+&b NA => None
op a+b NA => None
op &a-b NA => None
op a-&b NA => None
op a-b NA => None
add NA => None
sub NA => None
add_ref NA => None
sub_ref NA => None
cadd NA => re=None im=Some(InvalidArgument)
csub NA => re=None im=Some(InvalidArgument)
cadd_ref NA => re=None im=Some(InvalidArgument)
csub_ref NA => re=None im=Some(InvalidArgument)
op &a+b NB => None
op a+&b NB => None
op a+b NB => None
op &a-b NB => None
op a-&b NB => None
op a-b NB => None
add NB => None
sub NB => None
add_ref NB => None
sub_ref NB => None
cadd NB => re=None im=Some(DivisionByZero)
csub NB => re=None im=Some(DivisionByZero)
cadd_ref NB => re=None im=Some(DivisionByZero)
csub_ref NB => re=None im=Some(DivisionByZero)
op &a+b NN => None
op a+&b NN => None
op a+b NN => None
op &a-b NN => None
op a-&b NN => None
op a-b NN => None
add NN => None
sub NN => None
add_ref NN => None
sub_ref NN => None
cadd NN => re=None im=None
csub NN => re=None im=None
cadd_ref NN => re=None im=None
csub_ref NN => re=None im=None
val+A => Some(InvalidArgument)
A+val => Some(InvalidArgument)
&val+A => Some(InvalidArgument)
&val-A => Some(InvalidArgument)
inf+A => Some(InvalidArgument)
A-inf => Some(InvalidArgument)
val+B => Some(DivisionByZero)
B+val => Some(DivisionByZero)
&val+B => Some(DivisionByZero)
&val-B => Some(DivisionByZero)
inf+B => Some(DivisionByZero)
B-inf => Some(DivisionByZero)
val+N => None
N+val => None
&val+N => None
&val-N => None
inf+N => None
N-inf => None
asinh AA => re=Some(InvalidArgument) im=Some(InvalidArgument)
ln_gamma AA => re=Some(InvalidArgument) im=Some(InvalidArgument)
asinh AB => re=Some(InvalidArgument) im=Some(InvalidArgument)
ln_gamma AB => re=Some(InvalidArgument) im=Some(DivisionByZero)
asinh AN => re=Some(InvalidArgument) im=Some(InvalidArgument)
ln_gamma AN => re=Some(InvalidArgument) im=None
asinh BA => re=Some(DivisionByZero) im=Some(DivisionByZero)
ln_gamma BA => re=Some(DivisionByZero) im=Some(InvalidArgument)
asinh BB => re=Some(DivisionByZero) im=Some(DivisionByZero)
ln_gamma BB => re=Some(DivisionByZero) im=Some(DivisionByZero)
asinh BN => re=Some(DivisionByZero) im=Some(DivisionByZero)
ln_gamma BN => re=Some(DivisionByZero) im=None
asinh NA => re=None im=None
ln_gamma NA => re=None im=Some(InvalidArgument)
asinh NB => re=None im=None
ln_gamma NB => re=None im=Some(DivisionByZero)
asinh NN => re=None im=None
ln_gamma NN => re=None im=None
";

#[test]
fn both_nan_error_kind_is_pinned() {
    let got: Vec<String> = observations()
        .into_iter()
        .map(|(k, v)| format!("{k} => {v}"))
        .collect();
    let want: Vec<&str> = EXPECTED.lines().collect();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(&want) {
        assert_eq!(g, w);
    }
}
