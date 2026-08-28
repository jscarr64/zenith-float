//! Transcendental-only fixtures (compiled only from `transcendentals.rs`).

/// Batch size for transcendental kernels (each call is substantially more work).
pub const TRANSCENDENTAL_BATCH: usize = 64;

/// Positive operands for logarithms (domain is (0, +Inf)).
pub const LN_FIXTURES: &[&str] = &[
    "1.234567890123456789",
    "2.718281828459045235360287471",
    "10.0",
    "0.001234567890123456789",
    "1234.567890123456789",
    "3.141592653589793238462643383279",
    "1.000000000000000001",
    "0.999999999999999999",
    "987.654321098765432109",
    "42.0",
];

/// Arguments where `exp` stays in range without premature overflow at modest precision.
pub const EXP_FIXTURES: &[&str] = &[
    "0.01",
    "0.1",
    "0.5",
    "1.0",
    "2.3",
    "-0.5",
    "-1.0",
    "0.693147180559945309417232121",
    "2.302585092994045684017991455",
    "3.141592653589793238462643383279",
];

/// Radian arguments for sin / cos / tan (typical after range reduction).
pub const TRIG_FIXTURES: &[&str] = &[
    "0.1",
    "0.5",
    "1.0",
    "1.234567890123456789",
    "2.718281828459045235360287471",
    "3.141592653589793238462643383279",
    "-0.75",
    "-2.3",
    "0.999999999999999999",
    "1.570796326794896619231321691",
];
