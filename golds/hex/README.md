# Hex limb golds

`src/bin/hex_golds.rs` prints `ExactNum::to_bytes()` (and complex re/im) as hex for one representative of each special at `p=64`, `ToEven`.

The interchange format uses big-endian `u32` limbs, so the bytes are identical on x86_64 and aarch64 musl (64-bit targets only, P6).

Locked file: `golds/hex/reference.txt`.

```bash
./scripts/ci_hex_arm.sh     # aarch64-unknown-linux-musl + qemu-user
```
