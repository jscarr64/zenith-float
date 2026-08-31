# Hex limb golds

`src/bin/hex_golds.rs` prints `ExactNum::to_bytes()` (and complex re/im) as hex for one representative of each special at `p=64`, `ToEven`.

The interchange format uses big-endian `u32` limbs, so the bytes are identical on x86_64 (`WORD_BIT_SIZE=64`), i686 musl (`WORD_BIT_SIZE=32`), wasm32-wasip1, and aarch64 musl.

Locked file: `golds/hex/reference.txt`.

```bash
./scripts/ci_hex_32bit.sh   # i686-unknown-linux-musl, static
./scripts/ci_hex_wasm.sh    # wasm32-wasip1 + wasmtime
./scripts/ci_hex_arm.sh     # aarch64-unknown-linux-musl + qemu-user
```
