//! Deserialization of ExactNum.

use core::fmt::Formatter;
use core::str::FromStr;

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
        write!(formatter, "a decimal string or an integer")
    }

    fn visit_u64<E: Error>(self, v: u64) -> Result<Self::Value, E> {
        Ok(ExactNum::from(v))
    }

    fn visit_i64<E: Error>(self, v: i64) -> Result<Self::Value, E> {
        Ok(ExactNum::from(v))
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
}



