//! Serialization of software numbers and arrays.
//! Decimal strings carry an explicit `@p=` precision suffix.

use super::codec::{encode_exact_num, rm_name};
use crate::Ball;
use crate::ExactComplex;
use crate::ExactInt;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::ExactRational;
use crate::Ieee32Array;
use crate::Ieee64Array;
use crate::RoundingMode;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

impl Serialize for ExactNum {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&encode_exact_num(self))
    }
}

impl Serialize for ExactComplex {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut st = serializer.serialize_struct("ExactComplex", 2)?;
        st.serialize_field("re", self.re())?;
        st.serialize_field("im", self.im())?;
        st.end()
    }
}

impl Serialize for ExactRational {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut st = serializer.serialize_struct("ExactRational", 2)?;
        st.serialize_field("num", &encode_exact_num(self.num()))?;
        st.serialize_field("den", &encode_exact_num(self.den()))?;
        st.end()
    }
}

impl Serialize for ExactInt {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_dec_string())
    }
}

impl Serialize for Ball {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut st = serializer.serialize_struct("Ball", 2)?;
        st.serialize_field("mid", &encode_exact_num(self.mid()))?;
        st.serialize_field("rad", &encode_exact_num(self.rad()))?;
        st.end()
    }
}

impl Serialize for ExactNumArray {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let (rows, cols) = self.shape();
        let data: alloc::vec::Vec<alloc::string::String> =
            self.as_slice().iter().map(encode_exact_num).collect();
        let mut st = serializer.serialize_struct("ExactNumArray", 4)?;
        st.serialize_field("shape", &[rows, cols])?;
        st.serialize_field("data", &data)?;
        st.serialize_field("p", &self.precision())?;
        st.serialize_field("rm", rm_name(RoundingMode::ToEven))?;
        st.end()
    }
}

impl Serialize for Ieee64Array {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let (rows, cols) = self.shape();
        let mut st = serializer.serialize_struct("Ieee64Array", 2)?;
        st.serialize_field("shape", &[rows, cols])?;
        st.serialize_field("data", self.as_bits())?;
        st.end()
    }
}

impl Serialize for Ieee32Array {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let (rows, cols) = self.shape();
        let mut st = serializer.serialize_struct("Ieee32Array", 2)?;
        st.serialize_field("shape", &[rows, cols])?;
        st.serialize_field("data", self.as_bits())?;
        st.end()
    }
}
