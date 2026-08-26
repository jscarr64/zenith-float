//! Serialization of ExactNum.
//! Serialization to a string uses decimal radix.

use crate::ExactNum;
use serde::{Serialize, Serializer};

impl Serialize for ExactNum {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}


