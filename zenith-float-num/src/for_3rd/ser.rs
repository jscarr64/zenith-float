//! Serialization of [`ExactNum`] and [`ExactComplex`].
//! Numbers are written as decimal strings (`Display`).

use crate::ExactComplex;
use crate::ExactNum;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

impl Serialize for ExactNum {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
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
