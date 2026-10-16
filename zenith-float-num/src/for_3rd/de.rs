//! Deserialization of software numbers and arrays.

use super::codec::{decode_exact_num, parse_rm};
use crate::Ball;
use crate::ExactComplex;
use crate::ExactInt;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::ExactRational;
use crate::Ieee32;
use crate::Ieee32Array;
use crate::Ieee64;
use crate::Ieee64Array;
use crate::RoundingMode;
use core::fmt::Formatter;
use core::str::FromStr;
use serde::de::Error;
use serde::de::MapAccess;
use serde::de::Visitor;
use serde::{Deserialize, Deserializer};

struct ExactNumVisitor {}

impl<'de> Deserialize<'de> for ExactNum {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(ExactNumVisitor {})
    }
}

impl<'de> Visitor<'de> for ExactNumVisitor {
    type Value = ExactNum;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "a decimal string with optional @p= precision")
    }

    fn visit_u64<E: Error>(self, v: u64) -> Result<Self::Value, E> {
        Ok(ExactNum::from(v))
    }

    fn visit_i64<E: Error>(self, v: i64) -> Result<Self::Value, E> {
        Ok(ExactNum::from(v))
    }

    fn visit_u128<E: Error>(self, v: u128) -> Result<Self::Value, E> {
        Ok(ExactNum::from(v))
    }

    fn visit_i128<E: Error>(self, v: i128) -> Result<Self::Value, E> {
        Ok(ExactNum::from(v))
    }

    fn visit_str<E: Error>(self, v: &str) -> Result<Self::Value, E> {
        decode_exact_num(v)
            .or_else(|_| ExactNum::from_str(v).map_err(|e| Error::custom(format!("{e:?}"))))
    }

    fn visit_string<E: Error>(self, v: String) -> Result<Self::Value, E> {
        self.visit_str(&v)
    }
}

impl<'de> Deserialize<'de> for ExactComplex {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_struct("ExactComplex", &["re", "im"], ExactComplexVisitor {})
    }
}

struct ExactComplexVisitor {}

impl<'de> Visitor<'de> for ExactComplexVisitor {
    type Value = ExactComplex;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "a map with decimal-string fields re and im")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut re = None;
        let mut im = None;
        while let Some(key) = map.next_key::<&str>()? {
            match key {
                "re" => {
                    if re.is_some() {
                        return Err(Error::duplicate_field("re"));
                    }
                    re = Some(map.next_value()?);
                }
                "im" => {
                    if im.is_some() {
                        return Err(Error::duplicate_field("im"));
                    }
                    im = Some(map.next_value()?);
                }
                _ => {
                    let _: serde::de::IgnoredAny = map.next_value()?;
                }
            }
        }
        let re = re.ok_or_else(|| Error::missing_field("re"))?;
        let im = im.ok_or_else(|| Error::missing_field("im"))?;
        Ok(ExactComplex::new(re, im))
    }
}

impl<'de> Deserialize<'de> for ExactRational {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_struct("ExactRational", &["num", "den"], RatVisitor {})
    }
}

struct RatVisitor {}

impl<'de> Visitor<'de> for RatVisitor {
    type Value = ExactRational;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "a map with num and den decimal@p strings")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut num = None;
        let mut den = None;
        while let Some(key) = map.next_key::<&str>()? {
            match key {
                "num" => {
                    if num.is_some() {
                        return Err(Error::duplicate_field("num"));
                    }
                    let s: String = map.next_value()?;
                    num = Some(decode_exact_num(&s).map_err(Error::custom)?);
                }
                "den" => {
                    if den.is_some() {
                        return Err(Error::duplicate_field("den"));
                    }
                    let s: String = map.next_value()?;
                    den = Some(decode_exact_num(&s).map_err(Error::custom)?);
                }
                _ => {
                    let _: serde::de::IgnoredAny = map.next_value()?;
                }
            }
        }
        let num = num.ok_or_else(|| Error::missing_field("num"))?;
        let den = den.ok_or_else(|| Error::missing_field("den"))?;
        Ok(ExactRational::new(num, den))
    }
}

impl<'de> Deserialize<'de> for ExactInt {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_str(IntVisitor {})
    }
}

struct IntVisitor {}

impl<'de> Visitor<'de> for IntVisitor {
    type Value = ExactInt;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "a decimal integer string")
    }

    fn visit_str<E: Error>(self, v: &str) -> Result<Self::Value, E> {
        ExactInt::from_dec_string(v).ok_or_else(|| Error::custom("invalid ExactInt decimal"))
    }

    fn visit_string<E: Error>(self, v: String) -> Result<Self::Value, E> {
        self.visit_str(&v)
    }
}

impl<'de> Deserialize<'de> for Ball {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_struct("Ball", &["mid", "rad"], BallVisitor {})
    }
}

struct BallVisitor {}

impl<'de> Visitor<'de> for BallVisitor {
    type Value = Ball;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "a map with mid and rad decimal@p strings")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut mid = None;
        let mut rad = None;
        while let Some(key) = map.next_key::<&str>()? {
            match key {
                "mid" => {
                    if mid.is_some() {
                        return Err(Error::duplicate_field("mid"));
                    }
                    let s: String = map.next_value()?;
                    mid = Some(decode_exact_num(&s).map_err(Error::custom)?);
                }
                "rad" => {
                    if rad.is_some() {
                        return Err(Error::duplicate_field("rad"));
                    }
                    let s: String = map.next_value()?;
                    rad = Some(decode_exact_num(&s).map_err(Error::custom)?);
                }
                _ => {
                    let _: serde::de::IgnoredAny = map.next_value()?;
                }
            }
        }
        let mid = mid.ok_or_else(|| Error::missing_field("mid"))?;
        let rad = rad.ok_or_else(|| Error::missing_field("rad"))?;
        Ok(Ball::new(mid, rad))
    }
}

impl<'de> Deserialize<'de> for ExactNumArray {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_struct(
            "ExactNumArray",
            &["shape", "data", "p", "rm"],
            ArrVisitor {},
        )
    }
}

struct ArrVisitor {}

impl<'de> Visitor<'de> for ArrVisitor {
    type Value = ExactNumArray;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "shape, data, p, and rm")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut shape = None;
        let mut data = None;
        let mut p = None;
        let mut rm = RoundingMode::ToEven;
        while let Some(key) = map.next_key::<&str>()? {
            match key {
                "shape" => {
                    if shape.is_some() {
                        return Err(Error::duplicate_field("shape"));
                    }
                    shape = Some(map.next_value::<[usize; 2]>()?);
                }
                "data" => {
                    if data.is_some() {
                        return Err(Error::duplicate_field("data"));
                    }
                    let raw: alloc::vec::Vec<String> = map.next_value()?;
                    let mut vals = alloc::vec::Vec::new();
                    for s in raw {
                        vals.push(decode_exact_num(&s).map_err(Error::custom)?);
                    }
                    data = Some(vals);
                }
                "p" => {
                    if p.is_some() {
                        return Err(Error::duplicate_field("p"));
                    }
                    p = Some(map.next_value::<usize>()?);
                }
                "rm" => {
                    let s: String = map.next_value()?;
                    rm = parse_rm(&s).ok_or_else(|| Error::custom("unknown rounding mode"))?;
                }
                _ => {
                    let _: serde::de::IgnoredAny = map.next_value()?;
                }
            }
        }
        let [rows, cols] = shape.ok_or_else(|| Error::missing_field("shape"))?;
        let mut data = data.ok_or_else(|| Error::missing_field("data"))?;
        let p = p.ok_or_else(|| Error::missing_field("p"))?;
        let n = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::custom("shape overflow"))?;
        if n != data.len() {
            return Err(Error::custom("shape does not match data length"));
        }
        for v in &mut data {
            let _ = v.set_precision(p, rm);
        }
        ExactNumArray::from_shape(p, rows, cols, &data)
            .ok_or_else(|| Error::custom("ExactNumArray::from_shape rejected the buffer"))
    }
}

impl<'de> Deserialize<'de> for Ieee64Array {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_struct("Ieee64Array", &["shape", "data"], Ieee64Visitor {})
    }
}

struct Ieee64Visitor {}

impl<'de> Visitor<'de> for Ieee64Visitor {
    type Value = Ieee64Array;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "shape and u64 IEEE bit patterns")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut shape = None;
        let mut data = None;
        while let Some(key) = map.next_key::<&str>()? {
            match key {
                "shape" => shape = Some(map.next_value::<[usize; 2]>()?),
                "data" => data = Some(map.next_value::<alloc::vec::Vec<u64>>()?),
                _ => {
                    let _: serde::de::IgnoredAny = map.next_value()?;
                }
            }
        }
        let [rows, cols] = shape.ok_or_else(|| Error::missing_field("shape"))?;
        let bits = data.ok_or_else(|| Error::missing_field("data"))?;
        let n = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::custom("shape overflow"))?;
        if n != bits.len() {
            return Err(Error::custom("shape does not match data length"));
        }
        let vals: alloc::vec::Vec<Ieee64> = bits.into_iter().map(Ieee64::from_bits).collect();
        Ieee64Array::from_shape(rows, cols, &vals)
            .ok_or_else(|| Error::custom("Ieee64Array::from_shape rejected the buffer"))
    }
}

impl<'de> Deserialize<'de> for Ieee32Array {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_struct("Ieee32Array", &["shape", "data"], Ieee32Visitor {})
    }
}

struct Ieee32Visitor {}

impl<'de> Visitor<'de> for Ieee32Visitor {
    type Value = Ieee32Array;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "shape and u32 IEEE bit patterns")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut shape = None;
        let mut data = None;
        while let Some(key) = map.next_key::<&str>()? {
            match key {
                "shape" => shape = Some(map.next_value::<[usize; 2]>()?),
                "data" => data = Some(map.next_value::<alloc::vec::Vec<u32>>()?),
                _ => {
                    let _: serde::de::IgnoredAny = map.next_value()?;
                }
            }
        }
        let [rows, cols] = shape.ok_or_else(|| Error::missing_field("shape"))?;
        let bits = data.ok_or_else(|| Error::missing_field("data"))?;
        let n = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::custom("shape overflow"))?;
        if n != bits.len() {
            return Err(Error::custom("shape does not match data length"));
        }
        let vals: alloc::vec::Vec<Ieee32> = bits.into_iter().map(Ieee32::from_bits).collect();
        Ieee32Array::from_shape(rows, cols, &vals)
            .ok_or_else(|| Error::custom("Ieee32Array::from_shape rejected the buffer"))
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        Ball, ExactComplex, ExactInt, ExactNum, ExactNumArray, ExactRational, Ieee64, Ieee64Array,
        INF_NEG, INF_POS, NAN,
    };

    fn rt(n: &ExactNum) -> Result<ExactNum, Box<dyn std::error::Error>> {
        let s = serde_json::to_string(n)?;
        Ok(serde_json::from_str(&s)?)
    }

    #[test]
    fn serde_json_roundtrip() -> Result<(), Box<dyn std::error::Error>> {
        let n = ExactNum::from(42u32);
        let s = serde_json::to_string(&n)?;
        assert!(s.starts_with('"'), "serialized as a JSON string: {s}");
        assert!(s.contains("@p="), "precision suffix: {s}");
        let m: ExactNum = serde_json::from_str(&s)?;
        assert_eq!(n.compare(&m), Some(0));
        assert_eq!(n.precision(), m.precision());
        Ok(())
    }

    #[test]
    fn serde_json_integer_and_string() -> Result<(), Box<dyn std::error::Error>> {
        let from_int: ExactNum = serde_json::from_str("42")?;
        assert_eq!(from_int.compare(&ExactNum::from(42u32)), Some(0));

        let from_neg: ExactNum = serde_json::from_str("-7")?;
        assert_eq!(from_neg.compare(&ExactNum::from(-7i32)), Some(0));

        let from_str: ExactNum = serde_json::from_str("\"1.25e+0\"")?;
        assert!(!from_str.is_nan());
        assert!(from_str.is_positive());
        Ok(())
    }

    #[test]
    fn serde_json_specials() -> Result<(), Box<dyn std::error::Error>> {
        let inf = rt(&INF_POS)?;
        assert!(inf.is_inf_pos());
        let ninf = rt(&INF_NEG)?;
        assert!(ninf.is_inf_neg());
        let nan = rt(&NAN)?;
        assert!(nan.is_nan());
        Ok(())
    }

    #[test]
    fn serde_json_negative_fraction() -> Result<(), Box<dyn std::error::Error>> {
        let n: ExactNum = serde_json::from_str("\"-1.234567890123456789e-2\"")?;
        let m = rt(&n)?;
        assert_eq!(n.compare(&m), Some(0));
        assert!(m.is_negative());
        Ok(())
    }

    #[test]
    fn serde_json_complex() -> Result<(), Box<dyn std::error::Error>> {
        let z = ExactComplex::new(ExactNum::from(3u32), ExactNum::from(-4i32));
        let s = serde_json::to_string(&z)?;
        assert!(s.contains("\"re\""), "{s}");
        assert!(s.contains("\"im\""), "{s}");
        let w: ExactComplex = serde_json::from_str(&s)?;
        assert_eq!(z.re().compare(w.re()), Some(0));
        assert_eq!(z.im().compare(w.im()), Some(0));
        Ok(())
    }

    #[test]
    fn serde_json_reject_garbage() {
        let r: Result<ExactNum, _> = serde_json::from_str("true");
        assert!(r.is_err());
    }

    #[test]
    fn serde_rational_one_third() -> Result<(), Box<dyn std::error::Error>> {
        let r = ExactRational::from_i64(1, 3);
        let s = serde_json::to_string(&r)?;
        assert!(s.contains("@p="), "{s}");
        let q: ExactRational = serde_json::from_str(&s)?;
        assert_eq!(r.num().compare(q.num()), Some(0));
        assert_eq!(r.den().compare(q.den()), Some(0));
        Ok(())
    }

    #[test]
    fn serde_exact_int_and_array() -> Result<(), Box<dyn std::error::Error>> {
        let n = ExactInt::from_i64(-123);
        let s = serde_json::to_string(&n)?;
        let m: ExactInt = serde_json::from_str(&s)?;
        assert_eq!(n, m);

        let p = 256;
        let a = ExactNumArray::from_shape(
            p,
            2,
            2,
            &[
                ExactNum::from_u8(1, p),
                ExactNum::from_u8(2, p),
                ExactNum::from_u8(3, p),
                ExactNum::from_u8(4, p),
            ],
        )
        .ok_or("none value")?;
        let js = serde_json::to_string(&a)?;
        assert!(js.contains("\"p\""), "{js}");
        let b: ExactNumArray = serde_json::from_str(&js)?;
        assert_eq!(b.shape(), (2, 2));
        assert_eq!(b.precision(), p);
        for i in 0..4 {
            assert_eq!(
                a.get(i)
                    .ok_or("none value")?
                    .compare(b.get(i).ok_or("none value")?),
                Some(0)
            );
        }
        let bad = js.replacen("[2,2]", "[3,3]", 1);
        let r: Result<ExactNumArray, _> = serde_json::from_str(&bad);
        assert!(r.is_err());
        Ok(())
    }

    #[test]
    fn serde_ieee64_bits_and_ball() -> Result<(), Box<dyn std::error::Error>> {
        let a =
            Ieee64Array::from_values(&[Ieee64::from_i32(1), Ieee64::from_bits(0x7ff8000000000001)]);
        let js = serde_json::to_string(&a)?;
        let b: Ieee64Array = serde_json::from_str(&js)?;
        assert_eq!(a.as_bits(), b.as_bits());

        let mid = ExactNum::from_u8(1, 128);
        let rad = ExactNum::from_u8(1, 128).ldexp(-8, 128, crate::RoundingMode::ToEven);
        let ball = Ball::new(mid, rad);
        let s = serde_json::to_string(&ball)?;
        let ball2: Ball = serde_json::from_str(&s)?;
        assert_eq!(ball.mid().compare(ball2.mid()), Some(0));
        assert_eq!(ball.rad().compare(ball2.rad()), Some(0));
        Ok(())
    }
}
