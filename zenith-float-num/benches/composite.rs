//! Composite formula benchmarks — end-to-end workloads resembling Accumath evaluation.

mod shared;

use shared::{init_cc, sink, PRECISIONS};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use zenith_float_num::{Consts, ExactNum, RoundingMode};

/// π = 6 · arctan(1 / √3) — the README example, measured end-to-end.
fn bench_pi_via_atan(c: &mut Criterion) {
    let mut group = c.benchmark_group("composite/pi_via_atan");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                let three = ExactNum::from_word(3, p);
                let root3 = three.sqrt(p, rm);
                let inv = ExactNum::from_word(1, p).div(&root3, p, rm);
                let atan = inv.atan(p, rm, &mut cc);
                let pi = ExactNum::from_word(6, p).mul(&atan, p, rm);
                sink(pi);
            });
        });
    }
    group.finish();
}

/// e^x · ln(x) at a fixed operand — common in log-domain transforms.
fn bench_exp_ln_product(c: &mut Criterion) {
    let mut group = c.benchmark_group("composite/exp_ln_product");
    let rm = RoundingMode::ToEven;

    for &p in &[128, 1024, 4096] {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            let x = ExactNum::parse(
                "2.718281828459045235360287471",
                zenith_float_num::Radix::Dec,
                p,
                rm,
                &mut cc,
            );
            b.iter(|| {
                let ln_x = x.ln(p, rm, &mut cc);
                let exp_x = x.exp(p, rm, &mut cc);
                sink(exp_x.mul(&ln_x, p, rm));
            });
        });
    }
    group.finish();
}

/// Distance √(x² + y²) via hypot vs explicit sqrt — validates the fused path.
fn bench_pythagorean(c: &mut Criterion) {
    let mut group = c.benchmark_group("composite/pythagorean");
    let rm = RoundingMode::ToEven;

    for &p in &[128, 1024, 4096, 32768] {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = init_cc();
            let x = ExactNum::parse("3.141592653589793", zenith_float_num::Radix::Dec, p, rm, &mut cc);
            let y = ExactNum::parse("2.718281828459045", zenith_float_num::Radix::Dec, p, rm, &mut cc);
            b.iter(|| {
                sink(x.hypot(&y, p, rm));
            });
        });
    }
    group.finish();
}

/// sin²(x) + cos²(x) — paired trig with cancellation (identity should ≈ 1).
fn bench_trig_identity(c: &mut Criterion) {
    let mut group = c.benchmark_group("composite/trig_identity");
    let rm = RoundingMode::ToEven;

    for &p in &[128, 1024, 4096] {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            let x = ExactNum::parse(
                "1.234567890123456789",
                zenith_float_num::Radix::Dec,
                p,
                rm,
                &mut cc,
            );
            b.iter(|| {
                let s = x.sin(p, rm, &mut cc);
                let c = x.cos(p, rm, &mut cc);
                let s2 = s.mul(&s, p, RoundingMode::None);
                let c2 = c.mul(&c, p, RoundingMode::None);
                sink(s2.add(&c2, p, rm));
            });
        });
    }
    group.finish();
}

/// Decimal format of computed values — typical I/O path for Accumath results.
fn bench_decimal_format(c: &mut Criterion) {
    let mut group = c.benchmark_group("composite/decimal_format");
    let rm = RoundingMode::ToEven;

    for &p in &[128, 1024, 4096] {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            let x = ExactNum::parse(
                "3.141592653589793238462643383279",
                zenith_float_num::Radix::Dec,
                p,
                rm,
                &mut cc,
            );
            b.iter(|| {
                sink(x.format(zenith_float_num::Radix::Dec, rm, &mut cc).unwrap());
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_pi_via_atan,
    bench_exp_ln_product,
    bench_pythagorean,
    bench_trig_identity,
    bench_decimal_format,
);
criterion_main!(benches);
