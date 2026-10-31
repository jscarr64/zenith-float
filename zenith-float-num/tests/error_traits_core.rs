//! ZF-M16: `Error` and `QuadratureError` implement `core::error::Error` with and without the
//! `std` feature. The check is at compile time: this file does not build if either impl is
//! missing, so `cargo test --no-default-features` covers the `no_std` library.

use zenith_float_num::{Error, QuadratureError};

const fn assert_core_error<T: core::error::Error + 'static>() {}

const _: () = {
    assert_core_error::<Error>();
    assert_core_error::<QuadratureError>();
};

#[test]
fn error_types_implement_core_error() {
    assert_core_error::<Error>();
    assert_core_error::<QuadratureError>();
}
