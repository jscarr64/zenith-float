//! Deserialization of [`ExactNum`] and [`ExactComplex`].

use core::fmt::Formatter;
use core::str::FromStr;

use crate::ExactComplex;
use crate::ExactNum;
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
        write!(formatter, "a decimal string or an integer")
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
        ExactNum::from_str(v).map_err(|e| Error::custom(format!("{e:?}")))
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

#[cfg(test)]
mod tests {
    use crate::{ExactComplex, ExactNum, INF_NEG, INF_POS, NAN};

    fn rt(n: &ExactNum) -> ExactNum {
        let s = serde_json::to_string(n).expect("serialize");
        serde_json::from_str(&s).expect("deserialize")
    }

    #[test]
    fn serde_json_roundtrip() {
        let n = ExactNum::from(42u32);
        let s = serde_json::to_string(&n).expect("serialize");
        assert!(s.starts_with('"'), "serialized as a JSON string: {s}");
        let m: ExactNum = serde_json::from_str(&s).expect("deserialize");
        assert_eq!(n.cmp(&m), Some(0));
    }

    #[test]
    fn serde_json_integer_and_string() {
        let from_int: ExactNum = serde_json::from_str("42").expect("int");
        assert_eq!(from_int.cmp(&ExactNum::from(42u32)), Some(0));

        let from_neg: ExactNum = serde_json::from_str("-7").expect("neg int");
        assert_eq!(from_neg.cmp(&ExactNum::from(-7i32)), Some(0));

        let from_str: ExactNum =
            serde_json::from_str("\"1.25e+0\"").expect("decimal string");
        assert!(!from_str.is_nan());
        assert!(from_str.is_positive());
    }

    #[test]
    fn serde_json_specials() {
        let inf = rt(&INF_POS);
        assert!(inf.is_inf_pos());
        let ninf = rt(&INF_NEG);
        assert!(ninf.is_inf_neg());
        let nan = rt(&NAN);
        assert!(nan.is_nan());
    }

    #[test]
    fn serde_json_negative_fraction() {
        let n: ExactNum = serde_json::from_str("\"-1.234567890123456789e-2\"").unwrap();
        let m = rt(&n);
        assert_eq!(n.cmp(&m), Some(0));
        assert!(m.is_negative());
    }

    #[test]
    fn serde_json_complex() {
        let z = ExactComplex::new(ExactNum::from(3u32), ExactNum::from(-4i32));
        let s = serde_json::to_string(&z).expect("serialize complex");
        assert!(s.contains("\"re\""), "{s}");
        assert!(s.contains("\"im\""), "{s}");
        let w: ExactComplex = serde_json::from_str(&s).expect("deserialize complex");
        assert_eq!(z.re().cmp(w.re()), Some(0));
        assert_eq!(z.im().cmp(w.im()), Some(0));
    }

    #[test]
    fn serde_json_reject_garbage() {
        let r: Result<ExactNum, _> = serde_json::from_str("true");
        assert!(r.is_err());
    }
}
