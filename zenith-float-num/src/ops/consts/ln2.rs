// Derived from FLINT 3.6 `arb_const_log2` (Fredrik Johansson, LGPL-3.0-or-later)
// and, for the cache wrapper, astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! ln(2) by the FLINT 3.6 `arb_const_log2` hypergeometric series, cached.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::integer::ExactInt;
use crate::num::ExactNumNumber;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

/// `log2(3888) > 11.92`. Using 11 over-estimates the term count.
fn extra_bits(prec: usize) -> usize {
    if prec >= 2 {
        prec.ilog2() as usize + 1
    } else {
        1
    }
}

fn hypgeom_term_count(prec: usize) -> usize {
    let wp = prec.saturating_add(extra_bits(prec)).saturating_add(64);
    wp.div_ceil(11).saturating_add(16)
}

fn eval_a(n: u64) -> ExactInt {
    // A(n) = 1497 + 1794 n
    ExactInt::from_u64(1497).add(&ExactInt::from_u64(1794).mul_u64(n))
}

fn eval_p(n: u64) -> ExactInt {
    // P(n) = n (2n - 1)  (FLINT "3  0 -1 2")
    let nn = ExactInt::from_u64(n);
    nn.mul(&nn.mul_u64(2).sub(&ExactInt::one()))
}

fn eval_q(n: u64) -> ExactInt {
    // Q(n) = 7776 n^2 + 7776 n + 1080  (FLINT "3  1080 7776 7776")
    let nn = ExactInt::from_u64(n);
    let n2 = nn.mul(&nn);
    n2.mul_u64(7776)
        .add(&nn.mul_u64(7776))
        .add(&ExactInt::from_u64(1080))
}

/// Binary split of the FLINT `arb_const_log2_hypgeom` series on raw `ExactInt`.
/// `B = 1` identically, so it is omitted. Returns `(P, Q, T)`.
fn bsplit(a: usize, b: usize, cont: bool) -> (ExactInt, ExactInt, ExactInt) {
    if b - a == 1 {
        if a == 0 {
            let p = ExactInt::one();
            let q = ExactInt::one();
            let t = eval_a(0); // T = A(0) * P(0) = 1497
            return (p, q, t);
        }
        let n = a as u64;
        let p = eval_p(n);
        let q = eval_q(n);
        let t = eval_a(n).mul(&p);
        return (p, q, t);
    }
    let m = a + (b - a) / 2;
    let (p1, q1, t1) = bsplit(a, m, true);
    let (p2, q2, t2) = bsplit(m, b, true);
    // T = T1 Q2 + P1 T2
    let t = t1.mul(&q2).add(&p1.mul(&t2));
    let q = q1.mul(&q2);
    let p = if cont { p1.mul(&p2) } else { ExactInt::zero() };
    (p, q, t)
}

fn eval_log2(k: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
    // Exact ExactInt bsplit of the hypgeom series, then one division at
    // `wp = k + extra_bits(k) + 2 words`. Truncating intermediate merges at
    // dest precision (FLINT can, because `arb` carries a radius) dropped
    // low bits of the ExactNum midpoint from p=4096 upward: a fresh 4096-bit
    // ToEven agreed with FLINT in only 2495 bits, and cache reuse of an
    // 8192-bit fill produced a different 4096-bit value.
    let wp = k
        .saturating_add(extra_bits(k))
        .saturating_add(2 * WORD_BIT_SIZE);
    let n = hypgeom_term_count(wp);
    let (_p, q, t) = bsplit(0, n, false);
    // ln2 = T / (Q * 2160)
    let q2160 = q.mul_u64(2160);
    let t_num = t.to_exact_num_number(wp, RoundingMode::None)?;
    let q_num = q2160.to_exact_num_number(wp, RoundingMode::None)?;
    let mut val = t_num.div(&q_num, wp, RoundingMode::None)?;
    val.set_precision(k, rm)?;
    Ok(val)
}

/// Holds the currently computed ln(2).
#[derive(Debug)]
pub struct Ln2Cache {
    val: ExactNumNumber,
}

impl Ln2Cache {
    pub fn new() -> Result<Self, Error> {
        Ok(Ln2Cache {
            val: eval_log2(WORD_BIT_SIZE, RoundingMode::ToEven)?,
        })
    }

    /// Bits currently stored in the cache.
    pub(crate) fn cached_bit_len(&self) -> usize {
        self.val.mantissa_max_bit_len()
    }

    /// Return ln(2) with precision `k` (compute if the cache is short).
    pub(crate) fn for_prec(&mut self, k: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(k) + p_inc;

        loop {
            if self.val.mantissa_max_bit_len() >= p_wrk {
                // the cache keeps its full value; the caller gets it rounded to k, written straight from it
                if let Some(ret) = self.val.try_rounded_checked(k, rm, p_wrk)? {
                    return Ok(ret);
                }
            } else {
                self.val = eval_log2(p_wrk, RoundingMode::None)?;
                // the cache keeps its full value; the caller gets it rounded to k, written straight from it
                if let Some(ret) = self.val.try_rounded_checked(k, rm, p_wrk)? {
                    return Ok(ret);
                }
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, k)?;
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{RoundingMode, Sign};

    use super::*;

    #[test]
    #[cfg(target_pointer_width = "32")]
    fn test_ln2_const() -> Result<(), crate::Error> {
        let mut ln2 = Ln2Cache::new()?;
        let c = ln2.for_prec(3200, RoundingMode::ToEven)?;
        //println!("{:?}", c);
        let r = ExactNumNumber::from_raw_parts(
            &[
                3303612541, 717600284, 3199892741, 1462999970, 137959948, 6349977, 1794441986,
                1645437721, 1024185081, 96544981, 3471345337, 920660768, 4091103730, 196810490,
                3381928201, 1109492717, 2355549125, 1586635448, 1797039055, 1638923588, 2503522690,
                646601759, 4224649045, 2372882807, 579489529, 2708640721, 1931271310, 1151346004,
                816810139, 2530535123, 2409011252, 1433450182, 451048430, 1972937109, 3009774909,
                1293299123, 1348780427, 1599124760, 455979803, 126841693, 1818258609, 2922462427,
                2984344477, 2505920889, 389352748, 206047828, 1560181411, 122533377, 1578405506,
                1788641162, 895787831, 627481395, 3520465245, 1276780043, 2475905356, 3435941477,
                3027881267, 3376670514, 3683225829, 390465397, 335590294, 2100175838, 4229888533,
                3998653805, 2414170530, 1628254456, 4217109392, 133483025, 255841734, 3660421061,
                790684578, 1700766087, 942682696, 4125075133, 2640660682, 1926137777, 1985476427,
                628072684, 2973130811, 3119160259, 830224510, 449567249, 575334597, 1050069526,
                292141093, 660028201, 3241681220, 3979259445, 1257904912, 1435849467, 1844161688,
                3887625760, 2343238187, 2316113755, 1922610733, 1089684262, 66254511, 3387143064,
                3520035243, 2977044471,
            ],
            3200,
            Sign::Pos,
            0,
            false,
        )?;
        assert!(c.cmp(&r) == 0);
        Ok(())
    }

    #[test]
    #[cfg(not(target_pointer_width = "32"))]
    fn test_ln2_const() -> Result<(), crate::Error> {
        let mut ln2 = Ln2Cache::new()?;
        let c = ln2.for_prec(3200, RoundingMode::ToEven)?;
        //println!("{:?}", c);
        let r = ExactNumNumber::from_raw_parts(
            &[
                3082069754683924605,
                6283537028398873861,
                27272943683312140,
                7067101201094214402,
                414657537012126457,
                3954207892741588665,
                845294622150838770,
                4765234938047111433,
                6814547362189857733,
                7039123212900017103,
                2777133410944596354,
                10191454057530328917,
                11633523313888349945,
                4944993435491556494,
                10868565595481147547,
                6156621654544259124,
                8473700360670835694,
                5554677440240256317,
                6868188547772629387,
                544780923660251931,
                12551880549572046001,
                10762848267602590621,
                884968683061185836,
                526276848443620003,
                7682155296647843458,
                2695012071269245751,
                5483728532390938973,
                14757256277160841548,
                14502689430025391411,
                1677036114017882341,
                9020186540394984342,
                17174087324730849813,
                6993299640500441506,
                573305231163259792,
                15721388746840462790,
                7304734722601575330,
                17717062790720533064,
                8272698762445801674,
                2697551639276418891,
                13396691306361020475,
                1930876632637913214,
                4510014273271556293,
                2834799538024855589,
                17090789181815791940,
                6166926504001936144,
                16697225500131306648,
                9947632833883994667,
                4680158270178506285,
                14547668686819489455,
                12786308645202655659,
            ],
            3200,
            Sign::Pos,
            0,
            false,
        )?;
        assert!(c.cmp(&r) == 0);
        Ok(())
    }

    #[test]
    fn test_ln2_hypgeom_matches_published_digits() -> Result<(), crate::Error> {
        // OEIS A002162 prefix of ln(2), well past dest p = 3200 bits (~964 decimal digits
        // would be needed; 120 digits here cover p <= 384). The existing 3200-bit
        // from_raw_parts gold remains the high-precision check.
        const DEC: &str = "0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481020917637245";
        let mut cache = Ln2Cache::new()?;
        let mut cc = crate::Consts::new()?;
        for p in [64usize, 128, 256] {
            let got: crate::ExactNum = cache.for_prec(p, RoundingMode::ToEven)?.into();
            let want = crate::ExactNum::parse(
                DEC,
                crate::Radix::Dec,
                p.saturating_add(64),
                RoundingMode::ToEven,
                &mut cc,
            );
            let d = crate::ExactNum::sub_ref(&got, &want, p, RoundingMode::None).abs();
            let bits = if d.is_zero() {
                10_000
            } else {
                (want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)) as i32
            };
            assert!(
                bits >= (p as i32) - 4 || d.is_zero(),
                "ln2 cache p={p} only {bits} bits vs published digits"
            );
        }
        // Public ln(2) must hit the power-of-two path and match the cache.
        let p = 256usize;
        let public = crate::ExactNum::from_u8(2, p).ln(p, RoundingMode::ToEven, &mut cc);
        let cached: crate::ExactNum = cache.for_prec(p, RoundingMode::ToEven)?.into();
        let d = crate::ExactNum::sub_ref(&public, &cached, p, RoundingMode::None).abs();
        let bits = if d.is_zero() {
            10_000
        } else {
            (cached.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)) as i32
        };
        assert!(
            bits >= (p as i32) - 4 || d.is_zero(),
            "public ln(2) vs cache only {bits} bits"
        );
        Ok(())
    }

    #[test]
    fn test_ln2_matches_flint_at_row_working_precisions() -> Result<(), crate::Error> {
        // FLINT 3.6 / python-flint 0.9.0 `arb.const_log2` ToEven midpoints.
        // At 07d6ee8 a fresh p=4096 cache agreed with FLINT in only 2495 of 4096
        // bits; p=8192 in 3909 of 8192. p=2048 and the 3200-bit gold still matched.
        // Filling at 8192 and then reading 4096 produced a different 4096 midpoint.
        const LN2_2048: [u64; 32] = [
            6868188547772629387,
            544780923660251931,
            12551880549572046001,
            10762848267602590621,
            884968683061185836,
            526276848443620003,
            7682155296647843458,
            2695012071269245751,
            5483728532390938973,
            14757256277160841548,
            14502689430025391411,
            1677036114017882341,
            9020186540394984342,
            17174087324730849813,
            6993299640500441506,
            573305231163259792,
            15721388746840462790,
            7304734722601575330,
            17717062790720533064,
            8272698762445801674,
            2697551639276418891,
            13396691306361020475,
            1930876632637913214,
            4510014273271556293,
            2834799538024855589,
            17090789181815791940,
            6166926504001936144,
            16697225500131306648,
            9947632833883994667,
            4680158270178506285,
            14547668686819489455,
            12786308645202655659,
        ];
        {
            let mut ln2 = Ln2Cache::new()?;
            let c = ln2.for_prec(2048, RoundingMode::ToEven)?;
            let r = ExactNumNumber::from_raw_parts(&LN2_2048, 2048, Sign::Pos, 0, false)?;
            assert!(
                c.cmp(&r) == 0,
                "ln2 ToEven p=2048 disagrees with FLINT arb.const_log2"
            );
        }
        const LN2_4096: [u64; 64] = [
            1683981977765195225,
            17328443938231695182,
            10333551285590722846,
            14901130174119031709,
            6876852507673314734,
            7802240711941756352,
            10893311997856554306,
            770416509291317028,
            5131894646154344633,
            7044027212584627901,
            1544185093073130979,
            10743412093895996551,
            15520009427689038954,
            17235008467383040344,
            3082069754683924604,
            6283537028398873861,
            27272943683312140,
            7067101201094214402,
            414657537012126457,
            3954207892741588665,
            845294622150838770,
            4765234938047111433,
            6814547362189857733,
            7039123212900017103,
            2777133410944596354,
            10191454057530328917,
            11633523313888349945,
            4944993435491556494,
            10868565595481147547,
            6156621654544259124,
            8473700360670835694,
            5554677440240256317,
            6868188547772629387,
            544780923660251931,
            12551880549572046001,
            10762848267602590621,
            884968683061185836,
            526276848443620003,
            7682155296647843458,
            2695012071269245751,
            5483728532390938973,
            14757256277160841548,
            14502689430025391411,
            1677036114017882341,
            9020186540394984342,
            17174087324730849813,
            6993299640500441506,
            573305231163259792,
            15721388746840462790,
            7304734722601575330,
            17717062790720533064,
            8272698762445801674,
            2697551639276418891,
            13396691306361020475,
            1930876632637913214,
            4510014273271556293,
            2834799538024855589,
            17090789181815791940,
            6166926504001936144,
            16697225500131306648,
            9947632833883994667,
            4680158270178506285,
            14547668686819489455,
            12786308645202655659,
        ];
        {
            let mut ln2 = Ln2Cache::new()?;
            let c = ln2.for_prec(4096, RoundingMode::ToEven)?;
            let r = ExactNumNumber::from_raw_parts(&LN2_4096, 4096, Sign::Pos, 0, false)?;
            assert!(
                c.cmp(&r) == 0,
                "ln2 ToEven p=4096 disagrees with FLINT arb.const_log2"
            );
        }
        const LN2_8192: [u64; 128] = [
            15380334151373072787,
            4853205709837468559,
            15180201254144669436,
            17191618767403771779,
            13472633007291745206,
            10988803879252100845,
            9982350292508439189,
            16357665146209626244,
            2126312375026411415,
            9596014020852838733,
            6526060948148440112,
            3782886268850496926,
            4616937463712423083,
            16094152992545075298,
            12767136170006353273,
            9276362124443397780,
            12640865662083758843,
            1888533448826022573,
            17947010181381825034,
            14510358373783844601,
            18009594371835691210,
            7908467908886976714,
            7496286798650120630,
            7083903438220296811,
            8351280926043609319,
            13918765735208938155,
            3122446419593954965,
            10365386762090731137,
            17498909157122127143,
            2006951244055174130,
            1919935605167611636,
            15477994348553254657,
            12399739463684142264,
            1348239583426777193,
            11804095110535545172,
            15087814138645664243,
            9540087343598882522,
            13220844354109900748,
            12504881928627991056,
            8062879900679745266,
            9859934660139290250,
            2363745270234817426,
            1508931925449069755,
            5123958331046668426,
            1081491877773245656,
            17372500941954464729,
            3732979913381928924,
            1534707341129152403,
            6042401607091599235,
            15801603969349931038,
            13704111861162178889,
            912248912980489244,
            6146010341962107392,
            18089963405297836044,
            17180678772125609222,
            3913447253483908143,
            9906293889306082463,
            5956911102215554104,
            444478003528123916,
            4514934080037916381,
            7216354388686325020,
            14259934078564861152,
            13885076065611555834,
            9645420332068748557,
            1683981977765195224,
            17328443938231695182,
            10333551285590722846,
            14901130174119031709,
            6876852507673314734,
            7802240711941756352,
            10893311997856554306,
            770416509291317028,
            5131894646154344633,
            7044027212584627901,
            1544185093073130979,
            10743412093895996551,
            15520009427689038954,
            17235008467383040344,
            3082069754683924604,
            6283537028398873861,
            27272943683312140,
            7067101201094214402,
            414657537012126457,
            3954207892741588665,
            845294622150838770,
            4765234938047111433,
            6814547362189857733,
            7039123212900017103,
            2777133410944596354,
            10191454057530328917,
            11633523313888349945,
            4944993435491556494,
            10868565595481147547,
            6156621654544259124,
            8473700360670835694,
            5554677440240256317,
            6868188547772629387,
            544780923660251931,
            12551880549572046001,
            10762848267602590621,
            884968683061185836,
            526276848443620003,
            7682155296647843458,
            2695012071269245751,
            5483728532390938973,
            14757256277160841548,
            14502689430025391411,
            1677036114017882341,
            9020186540394984342,
            17174087324730849813,
            6993299640500441506,
            573305231163259792,
            15721388746840462790,
            7304734722601575330,
            17717062790720533064,
            8272698762445801674,
            2697551639276418891,
            13396691306361020475,
            1930876632637913214,
            4510014273271556293,
            2834799538024855589,
            17090789181815791940,
            6166926504001936144,
            16697225500131306648,
            9947632833883994667,
            4680158270178506285,
            14547668686819489455,
            12786308645202655659,
        ];
        {
            let mut ln2 = Ln2Cache::new()?;
            let c = ln2.for_prec(8192, RoundingMode::ToEven)?;
            let r = ExactNumNumber::from_raw_parts(&LN2_8192, 8192, Sign::Pos, 0, false)?;
            assert!(
                c.cmp(&r) == 0,
                "ln2 ToEven p=8192 disagrees with FLINT arb.const_log2"
            );
        }
        // Cache reuse must not change a lower ToEven request.
        let mut reuse = Ln2Cache::new()?;
        let _ = reuse.for_prec(8192, RoundingMode::ToEven)?;
        let c = reuse.for_prec(4096, RoundingMode::ToEven)?;
        let r = ExactNumNumber::from_raw_parts(&LN2_4096, 4096, Sign::Pos, 0, false)?;
        assert!(
            c.cmp(&r) == 0,
            "ln2 ToEven p=4096 after an 8192 fill disagrees with FLINT"
        );
        Ok(())
    }
}
