//! Bare-metal runner: `thumbv7em-none-eabihf` (Cortex-M4F) under
//! `qemu-system-arm -machine mps2-an386 -semihosting`. Output and exit status via semihosting.
#![no_std]
#![no_main]

extern crate alloc;

use core::mem::MaybeUninit;
use cortex_m_rt::entry;
use cortex_m_semihosting::{debug, hprintln};
use embedded_alloc::LlffHeap as Heap;

#[global_allocator]
static HEAP: Heap = Heap::empty();

/// Heap size (bytes). The MPS2-AN386 model has 4 MiB of RAM at 0x2000_0000.
const HEAP_SIZE: usize = 3 * 1024 * 1024;

#[panic_handler]
fn panic(info: &core::panic::PanicInfo) -> ! {
    hprintln!("PANIC {}", info);
    debug::exit(debug::EXIT_FAILURE);
    loop {}
}

#[entry]
fn main() -> ! {
    static mut HEAP_MEM: [MaybeUninit<u8>; HEAP_SIZE] = [MaybeUninit::uninit(); HEAP_SIZE];
    // `#[entry]` turns this `static mut` into a `&'static mut` handed out exactly once.
    unsafe { HEAP.init(HEAP_MEM.as_mut_ptr() as usize, HEAP_SIZE) }
    hprintln!("zenith-float no_std reference tests (thumbv7em-none-eabihf, qemu mps2-an386)");
    let (ok, bad) = zenith_float_nostd_tests::run(&mut |line| hprintln!("{}", line));
    hprintln!("no_std result: {} passed, {} failed", ok, bad);
    debug::exit(if bad == 0 { debug::EXIT_SUCCESS } else { debug::EXIT_FAILURE });
    // Only reached if the host ignores the semihosting exit.
    loop {
        cortex_m::asm::wfi();
    }
}
