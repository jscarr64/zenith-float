//! Host runner for the `no_std` reference cases (Jeff 2026-11-24 4:40 PM CT item 6: test code is safe Rust only).
//!
//! The cases ([`zenith_float_nostd_tests`], `#![no_std]`, core + alloc) and `zenith-float-num` (built with
//! `default-features = false`) are compiled without `std`; only this binary uses `std`, for its allocator, output and
//! exit status. It used to be a `#![no_main]` binary linking libc by hand (a `GlobalAlloc` over `posix_memalign`,
//! `write`, `_exit`), which needed five raw FFI blocks. The no-`std` property it guarded is checked where it is
//! decided: a `#![no_std]` crate cannot name `std` at compile time, and `scripts/ci_nostd.sh` builds both libraries
//! for the bare-metal `x86_64-unknown-none` and `aarch64-unknown-none` targets, which have no `std` to link.
//!
//! Exit status: the number of failed cases (capped at 255), `0` when every reference case passes.
#![forbid(unsafe_code)]

fn main() {
    println!("zenith-float no_std reference tests (host runner)");
    let (ok, bad) = zenith_float_nostd_tests::run(&mut |line| println!("{line}"));
    println!("no_std result: {ok} passed, {bad} failed");
    std::process::exit(i32::try_from(bad.min(255)).unwrap_or(255));
}
