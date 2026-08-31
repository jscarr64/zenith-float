//! Deterministic decimal codec: value plus explicit precision (`s@p=N`).

use crate::defs::DEFAULT_P;
use crate::Consts;
use crate::ExactNum;
use crate::Radix;
use crate::RoundingMode;
use alloc::string::String;
use alloc::string::ToString;

const P_MARK: &str = "@p=";

pub(super) fn encode_exact_num(n: &ExactNum) -> String {
    let p = n.precision().unwrap_or(DEFAULT_P);
    let mut s = n.to_string();
    s.push_str(P_MARK);
    s.push_str(&p.to_string());
    s
}

pub(super) fn decode_exact_num(s: &str) -> Result<ExactNum, String> {
    let (body, p) = split_prec(s);
    let mut cc = Consts::new().map_err(|e| alloc::format!("{e:?}"))?;
    let n = ExactNum::parse(body, Radix::Dec, p, RoundingMode::ToEven, &mut cc);
    if n.is_nan() && !body.eq_ignore_ascii_case("nan") && !body.eq_ignore_ascii_case("err") {
        return Err(alloc::format!("invalid exact decimal: {body}"));
    }
    Ok(n)
}

fn split_prec(s: &str) -> (&str, usize) {
    if let Some(i) = s.rfind(P_MARK) {
        let body = &s[..i];
        let p = s[i + P_MARK.len()..].parse::<usize>().unwrap_or(DEFAULT_P);
        (body, p.max(1))
    } else {
        (s, DEFAULT_P)
    }
}

pub(super) fn rm_name(rm: RoundingMode) -> &'static str {
    match rm {
        RoundingMode::None => "None",
        RoundingMode::Up => "Up",
        RoundingMode::Down => "Down",
        RoundingMode::ToZero => "ToZero",
        RoundingMode::FromZero => "FromZero",
        RoundingMode::ToEven => "ToEven",
        RoundingMode::ToOdd => "ToOdd",
    }
}

pub(super) fn parse_rm(s: &str) -> Option<RoundingMode> {
    Some(match s {
        "None" => RoundingMode::None,
        "Up" => RoundingMode::Up,
        "Down" => RoundingMode::Down,
        "ToZero" => RoundingMode::ToZero,
        "FromZero" => RoundingMode::FromZero,
        "ToEven" => RoundingMode::ToEven,
        "ToOdd" => RoundingMode::ToOdd,
        _ => return None,
    })
}
