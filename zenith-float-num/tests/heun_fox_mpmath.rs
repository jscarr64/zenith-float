//! Meijer G / Fox H / Heun vs mpmath 1.4.1 (50-digit golds) and independent reductions.

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

fn check(
    label: &str,
    got: ExactNum,
    want: &str,
    p: usize,
    cc: &mut Consts,
    failures: &mut Vec<String>,
) {
    let w = parse(want, 4 * p, cc);
    let bts = bits(&got, &w, p);
    if bts < (p as i32) - 8 {
        failures.push(format!("{label} p={p}: {bts} bits {:?}", got.err()));
    }
}

fn dec(s: &str, p: usize, cc: &mut Consts) -> ExactNum {
    parse(s, p, cc)
}

#[test]
fn heun_fox_mpmath() {
    let mut cc = Consts::new().unwrap();
    let mut failures = Vec::new();
    // 50-digit golds are ~166 bits; 128-bit dest is the honest check.
    for &p in &[128usize] {
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let half = one.div(&ExactNum::from_u8(2, p), p, RM);
        let two = ExactNum::from_u8(2, p);
        let three = ExactNum::from_u8(3, p);

        check(
            "G^{1,0}_{0,1}(1/2)=e^{-1/2}",
            half.meijer_g(&[], &[], &[zero.clone()], &[], p, RM, &mut cc),
            "0.60653065971263342360379953499118045344191813548719",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "G^{1,0}_{0,1}(2)=e^{-2}",
            two.meijer_g(&[], &[], &[zero.clone()], &[], p, RM, &mut cc),
            "0.13533528323661269189399949497248440340763154590958",
            p,
            &mut cc,
            &mut failures,
        );
        let a03 = dec("0.3", p, &mut cc);
        check(
            "G^{1,1}_{1,1}(1/2|0.3;0)",
            half.meijer_g(&[a03.clone()], &[], &[zero.clone()], &[], p, RM, &mut cc),
            "0.97730320798596565032892615428324051950104280073779",
            p,
            &mut cc,
            &mut failures,
        );
        // series-2 / |z|>1
        check(
            "G^{1,1}_{1,1}(3|0.3;0)",
            three.meijer_g(&[a03.clone()], &[], &[zero.clone()], &[], p, RM, &mut cc),
            "0.49187099298526723194079416580486746098307358094869",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "G^{2,0}_{0,2}(0.3|0,1/2)",
            dec("0.3", p, &mut cc).meijer_g(
                &[],
                &[],
                &[zero.clone(), half.clone()],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "0.59269213972559885103639876208967425787165788448217",
            p,
            &mut cc,
            &mut failures,
        );
        // awkward half-integer poles
        let qtr = one.div(&ExactNum::from_u8(4, p), p, RM);
        check(
            "G^{2,1}_{1,2}(0.15|1/2; 1/4,-1/4)",
            dec("0.15", p, &mut cc).meijer_g(
                &[half.clone()],
                &[],
                &[qtr.clone(), qtr.neg()],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "8.2338974235222697082569708896191469404099880964535",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "G^{1,1}_{1,2}(0.3|0.4; 0.1; 0.6)",
            dec("0.3", p, &mut cc).meijer_g(
                &[dec("0.4", p, &mut cc)],
                &[],
                &[dec("0.1", p, &mut cc)],
                &[dec("0.6", p, &mut cc)],
                p,
                RM,
                &mut cc,
            ),
            "0.41830601347346321883606696857111054341233521506556",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "G^{1,0}_{1,2}(0.8| ;1 ; 0; 1/2)",
            dec("0.8", p, &mut cc).meijer_g(
                &[],
                &[one.clone()],
                &[zero.clone()],
                &[half.clone()],
                p,
                RM,
                &mut cc,
            ),
            "0.564189583547756286948079451560772585844050629329",
            p,
            &mut cc,
            &mut failures,
        );
        // J_{1/4}(1.2) = G^{1,0}_{0,2}((0.6)^2 | ; ν/2; −ν/2)
        let nu = qtr.clone();
        let x = dec("1.2", p, &mut cc);
        let arg = x.div(&two, p, RM).mul(&x.div(&two, p, RM), p, RM);
        check(
            "G Bessel J_{1/4}(1.2)",
            arg.meijer_g(
                &[],
                &[],
                &[nu.div(&two, p, RM)],
                &[nu.neg().div(&two, p, RM)],
                p,
                RM,
                &mut cc,
            ),
            "0.71291095222750948038059963500822076515883269585735",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "G^{1,1}_{1,1}(5|0.4;0) series-2",
            ExactNum::from_u8(5, p).meijer_g(
                &[dec("0.4", p, &mut cc)],
                &[],
                &[zero.clone()],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "0.50822967193437512115335211822765588504313703655064",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "G^{1,1}_{1,1}(0.2|1.7;0) (negative)",
            dec("0.2", p, &mut cc).meijer_g(
                &[dec("1.7", p, &mut cc)],
                &[],
                &[zero.clone()],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "-4.8554317586620597160439431116555107534868625190458",
            p,
            &mut cc,
            &mut failures,
        );

        check(
            "fox_h A=B=1 → e^{-1/2}",
            half.fox_h(
                &[],
                &[],
                &[(zero.clone(), one.clone())],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "0.60653065971263342360379953499118045344191813548719",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "fox_h (0,2)",
            half.fox_h(
                &[],
                &[],
                &[(zero.clone(), two.clone())],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "0.24653434569761989392297869030758873761791888156821",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "fox_h (0,3) z=0.4",
            dec("0.4", p, &mut cc).fox_h(
                &[],
                &[],
                &[(zero.clone(), three.clone())],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "0.15954670338989146420208840484952103429105746090809",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "fox_h [[(0.3,1)],[(0,1)]]",
            half.fox_h(
                &[(a03.clone(), one.clone())],
                &[],
                &[(zero.clone(), one.clone())],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "0.97730320798596565032892615428324051950104280073779",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "fox_h (0,1/2) rational scale",
            half.fox_h(
                &[],
                &[],
                &[(zero.clone(), half.clone())],
                &[],
                p,
                RM,
                &mut cc,
            ),
            "1.5576015661428097364903405339566412945935445808523",
            p,
            &mut cc,
            &mut failures,
        );

        // HeunG → 2F1 (independent: mpmath hyp2f1)
        let alpha = one.div(&three, p, RM);
        let beta = two.div(&ExactNum::from_u8(5, p), p, RM);
        let gamma = ExactNum::from_u8(4, p).div(&three, p, RM);
        let zz = dec("0.35", p, &mut cc);
        check(
            "HeunG 2F1(1/3,2/5;4/3;0.35)",
            zz.heun_g(
                &dec("2.5", p, &mut cc),
                &dec("2.5", p, &mut cc).mul(&alpha, p, RM).mul(&beta, p, RM),
                &alpha,
                &beta,
                &gamma,
                &alpha.add(&beta, p, RM).sub(&gamma, p, RM).add(&one, p, RM),
                p,
                RM,
                &mut cc,
            ),
            "1.0411554332229463578669469840254100095847130753819",
            p,
            &mut cc,
            &mut failures,
        );
        // awkward accessory, not a 2F1 reduction
        check(
            "HeunG(2,0.3;0.4,0.5,0.8,0.6;0.25)",
            dec("0.25", p, &mut cc).heun_g(
                &two,
                &dec("0.3", p, &mut cc),
                &dec("0.4", p, &mut cc),
                &half,
                &dec("0.8", p, &mut cc),
                &dec("0.6", p, &mut cc),
                p,
                RM,
                &mut cc,
            ),
            "1.0532709966795252110568734749696132676374130054028",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "HeunG(1.5,-0.2;0.3,-0.4,1.1,0.7;0.4)",
            dec("0.4", p, &mut cc).heun_g(
                &dec("1.5", p, &mut cc),
                &dec("-0.2", p, &mut cc),
                &dec("0.3", p, &mut cc),
                &dec("-0.4", p, &mut cc),
                &dec("1.1", p, &mut cc),
                &dec("0.7", p, &mut cc),
                p,
                RM,
                &mut cc,
            ),
            "0.94462241427944807688964058137831711175591623829119",
            p,
            &mut cc,
            &mut failures,
        );

        // HeunC → 1F1 (independent: mpmath hyp1f1)
        check(
            "HeunC 1F1(-0.4;1.3;0.35)",
            dec("0.35", p, &mut cc).heun_c(
                &dec("0.4", p, &mut cc),
                &dec("1.3", p, &mut cc),
                &zero,
                &one.neg(),
                &dec("0.4", p, &mut cc),
                p,
                RM,
                &mut cc,
            ),
            "0.88709776171325453338365695154583842259203290747345",
            p,
            &mut cc,
            &mut failures,
        );
        check(
            "HeunC 1F1(0.7;2.2;0.6) (α=-0.7)",
            dec("0.6", p, &mut cc).heun_c(
                &dec("-0.7", p, &mut cc),
                &dec("2.2", p, &mut cc),
                &zero,
                &one.neg(),
                &dec("-0.7", p, &mut cc),
                p,
                RM,
                &mut cc,
            ),
            "1.2257058951900165508102217837364816100227899890624",
            p,
            &mut cc,
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
