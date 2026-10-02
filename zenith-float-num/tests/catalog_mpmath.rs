//! Catalog specials vs mpmath 1.4.1 at 50 digits / 2p+64 bits.

use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

const RM: RoundingMode = RoundingMode::ToEven;

fn parse(s: &str, p: usize, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(s, Radix::Dec, p, RM, cc)
}

fn bits(got: &ExactNum, want: &ExactNum, p: usize) -> i32 {
    if got.is_nan() {
        return i32::MIN;
    }
    let d = got.sub(want, 2 * p, RoundingMode::None);
    if d.is_zero() {
        return i32::MAX;
    }
    want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)
}

fn check(label: &str, got: ExactNum, want: &str, p: usize, cc: &mut Consts, failures: &mut Vec<String>) {
    let w = parse(want, 4 * p, cc);
    let bts = bits(&got, &w, p);
    if bts < (p as i32) - 8 {
        failures.push(format!("{label} p={p}: {bts} bits {:?}", got.err()));
    }
}

#[test]
fn catalog_mpmath() {
    let mut cc = Consts::new().unwrap();
    let mut failures = Vec::new();
    // 50-digit golds are ~166 bits, so 128-bit dest is the honest check.
    for &p in &[128usize] {
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RM);
        let x2 = ExactNum::from_u8(2, p);
        let x15 = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, RM);

        check(
            "Gi(1/2)",
            half.scorer_gi(p, RM, &mut cc),
            "0.24472104327655819769105392735646002370770999517103",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "Gi(2)",
            x2.scorer_gi(p, RM, &mut cc),
            "0.16895356565401036277372360565655438323513630836014",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "Gi(-2)",
            x2.neg().scorer_gi(p, RM, &mut cc),
            "-0.55325158419788968902485190507151282381293712677833",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "Hi(1/2)",
            half.scorer_hi(p, RM, &mut cc),
            "0.60955599982659729560899487143878315714907740533371",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "struve_h(0,1)",
            ExactNum::from_u8(1, p).struve_h(&ExactNum::new(p), p, RM, &mut cc),
            "0.56865662704828795098642288632235327430267782725556",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "cl2(pi/2)",
            cc.pi(p, RM).ldexp(-1, p, RM).clausen_cl2(p, RM, &mut cc),
            "0.91596559417721901505460351493238411077414937428167",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "barnes_g(3.5)",
            ExactNum::from_u8(7, p).div(&ExactNum::from_u8(2, p), p, RM).barnes_g(p, RM, &mut cc),
            "1.2596482574951921440863079510960698255629857900463",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "polygamma(1,1.5)",
            x15.polygamma(1, p, RM, &mut cc),
            "0.93480220054467930941724549993807556765684970362040",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "hurwitz(3,1/4)",
            ExactNum::from_u8(3, p).hurwitz_zeta(
                &ExactNum::from_u8(1, p).div(&ExactNum::from_u8(4, p), p, RM),
                p,
                RM,
                &mut cc,
            ),
            "64.663869968768460166668983589421994943644904751419",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "lambert_w0(1.5)",
            x15.lambert_w0(p, RM, &mut cc),
            "0.72586135776622625704868939927630687970616284950429",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "1F1(1;2;1/2)",
            half.hypergeom_1f1(&ExactNum::from_u8(1, p), &ExactNum::from_u8(2, p), p, RM, &mut cc),
            "1.2974425414002562936973015756283271433075522014203",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "ber(1.5)",
            x15.kelvin_ber(p, RM, &mut cc),
            "0.92107218354625576412221298628004717462121862008445",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "bei(1.5)",
            x15.kelvin_bei(p, RM, &mut cc),
            "0.55756006230308669489422237615617643459061643284434",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "polylog2(0.3)",
            parse("0.3", p, &mut cc).polylog(2, p, RM, &mut cc),
            "0.32612951007547606953003569417499604570558867999792",
            p,
            &mut cc,
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
