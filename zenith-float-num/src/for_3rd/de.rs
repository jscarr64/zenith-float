//! Deserialization of ExactNum.

use core::fmt::Formatter;
use core::str::FromStr;

use crate::num::ExactNumNumber;
use crate::ExactNum;
use serde::de::Error;
use serde::de::Visitor;
use serde::{Deserialize, Deserializer};

pub struct ExactNumVisitor {}

impl<'de> Deserialize<'de> for ExactNum {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(ExactNumVisitor {})
    }
}

impl<'de> Visitor<'de> for ExactNumVisitor {
    type Value = ExactNum;

    fn expecting(&self, formatter: &mut Formatter) -> core::fmt::Result {
        write!(formatter, "except `String`, `Number`, `Bytes`")
    }

    fn visit_u64<E: Error>(self, v: u64) -> Result<Self::Value, E> {
        match ExactNumNumber::from_usize(v as usize) {
            Ok(o) => Ok(o.into()),
            Err(e) => Err(Error::custom(format!("{e:?}"))),
        }
    }

    fn visit_str<E: Error>(self, v: &str) -> Result<Self::Value, E> {
        match ExactNum::from_str(v) {
            Ok(o) => Ok(o),
            Err(e) => Err(Error::custom(format!("{e:?}"))),
        }
    }

    fn visit_string<E: Error>(self, v: String) -> Result<Self::Value, E> {
        self.visit_str(&v)
    }

    // lossless conversion
    // (&[Word], usize, Sign, Exponent)
    // (s * len, s    , 1   , 1       )
    // fn visit_bytes<E: Error>(self, _: &[u8]) -> Result<Self::Value, E> {
    //     todo!()
    // }
}



