//! `no_std` runner for hosted Linux targets: links only libc (no Rust `std`).
//!
//! Allocator: `posix_memalign` / `free`. Output: `write(1, ..)`. Exit status: number of failures
//! (capped at 255), `0` when every reference case passes.
#![no_std]
#![no_main]

extern crate alloc;

use core::alloc::{GlobalAlloc, Layout};
use core::ffi::{c_char, c_int, c_void};

#[link(name = "c")]
extern "C" {
    fn posix_memalign(memptr: *mut *mut c_void, alignment: usize, size: usize) -> c_int;
    fn free(ptr: *mut c_void);
    fn write(fd: c_int, buf: *const c_void, count: usize) -> isize;
    fn _exit(status: c_int) -> !;
}

// The prebuilt `liballoc` for hosted targets is compiled with unwinding tables and references
// these two symbols even though this binary is built with `panic = "abort"` and never unwinds:
// `_Unwind_Resume` comes from libgcc_s; the personality routine is never called.
#[link(name = "gcc_s")]
extern "C" {}

#[no_mangle]
extern "C" fn rust_eh_personality() {}

struct LibcAlloc;

unsafe impl GlobalAlloc for LibcAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let align = layout.align().max(core::mem::size_of::<usize>());
        let mut p: *mut c_void = core::ptr::null_mut();
        if posix_memalign(&mut p, align, layout.size().max(1)) != 0 {
            return core::ptr::null_mut();
        }
        p.cast()
    }
    unsafe fn dealloc(&self, ptr: *mut u8, _layout: Layout) {
        free(ptr.cast());
    }
}

#[global_allocator]
static ALLOC: LibcAlloc = LibcAlloc;

fn out(s: &str) {
    let mut b = s.as_bytes();
    while !b.is_empty() {
        let n = unsafe { write(1, b.as_ptr().cast(), b.len()) };
        if n <= 0 {
            return;
        }
        b = &b[n as usize..];
    }
}

#[panic_handler]
fn panic(info: &core::panic::PanicInfo) -> ! {
    use core::fmt::Write;
    struct W;
    impl core::fmt::Write for W {
        fn write_str(&mut self, s: &str) -> core::fmt::Result {
            out(s);
            Ok(())
        }
    }
    let _ = writeln!(W, "PANIC {info}");
    unsafe { _exit(101) }
}

#[no_mangle]
pub extern "C" fn main(_argc: c_int, _argv: *const *const c_char) -> c_int {
    out("zenith-float no_std reference tests (host, libc only)\n");
    let (ok, bad) = zenith_float_nostd_tests::run(&mut |line| {
        out(line);
        out("\n");
    });
    let summary = alloc::format!("no_std result: {ok} passed, {bad} failed\n");
    out(&summary);
    bad.min(255) as c_int
}
