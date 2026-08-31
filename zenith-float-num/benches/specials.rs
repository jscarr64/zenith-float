//! Special-function Criterion benches at 64 / 128 / 256 / 512 / 1024 bits.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;
use zenith_float_num::{Consts, ExactNum, RoundingMode};

const SPECIAL_PRECISIONS: &[usize] = &[64, 128, 256, 512, 1024];

fn sink<T>(v: T) {
    black_box(v);
}

fn bench_erf(c: &mut Criterion) {
    let mut group = c.benchmark_group("specials/erf");
    let rm = RoundingMode::ToEven;
    for &p in SPECIAL_PRECISIONS {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            let x = ExactNum::from_u8(1, p);
            b.iter(|| sink(x.erf(p, rm, &mut cc)));
        });
    }
    group.finish();
}

fn bench_gamma(c: &mut Criterion) {
    let mut group = c.benchmark_group("specials/gamma");
    let rm = RoundingMode::ToEven;
    for &p in SPECIAL_PRECISIONS {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            let x = ExactNum::from_u8(5, p);
            b.iter(|| sink(x.gamma(p, rm, &mut cc)));
        });
    }
    group.finish();
}

fn bench_bessel_j(c: &mut Criterion) {
    let mut group = c.benchmark_group("specials/bessel_j");
    let rm = RoundingMode::ToEven;
    for &p in SPECIAL_PRECISIONS {
        group.bench_with_input(BenchmarkId::from_parameter(p), &p, |b, &p| {
            let mut cc = Consts::new().expect("constants cache");
            let x = ExactNum::from_u8(1, p);
            b.iter(|| sink(x.bessel_j(0, p, rm, &mut cc)));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_erf, bench_gamma, bench_bessel_j);
criterion_main!(benches);
