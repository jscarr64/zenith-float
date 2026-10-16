//! Puts `memory.x` on the linker search path for the bare-metal runner (cortex-m-rt `link.x`).
use std::{env, fs, path::PathBuf};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = PathBuf::from(env::var_os("OUT_DIR").ok_or("OUT_DIR is not set")?);
    fs::copy("memory.x", out.join("memory.x"))?;
    println!("cargo:rustc-link-search={}", out.display());
    println!("cargo:rerun-if-changed=memory.x");
    println!("cargo:rerun-if-changed=build.rs");
    Ok(())
}
