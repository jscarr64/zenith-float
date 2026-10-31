# ZF-M32 decision (Paul, 2026-10-30 6:37 PM CT, under P1/P6; Spence recommendation 6:36 PM)
- MAX_PRECISION = 2^32-1 bits on every target (zf defs.rs and zb precision.rs read the same constant).
- One up-front size check before any allocation: a single buffer larger than 2^31 bytes (2 GiB) is refused with Error::MemoryAllocation, identically on every OS/arch, never left to the allocator/overcommit.
  (Max mantissa 2^32 bits = 512 MiB; 4p working buffer = 2 GiB fits exactly.)
- Residual, documented in LIBRARY.md: a machine with less free RAM than the work needs can still be OOM-killed by the OS; that is physical, not a math answer.
- M25/M30 cost-model overflow stays live (u128 or clamp).
- Red rows: p=2^32-1 accepted+correct, p=2^32 InvalidArgument, budget-1 and budget+1 byte boundary, all on host + aarch64.
