//! Linear-algebra and FFT Criterion benches (software limbs / software IEEE).

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;
use zenith_float_num::{Consts, ExactNum, ExactNumArray, Ieee64, Ieee64Array, RoundingMode};

fn sink<T>(v: T) {
    black_box(v);
}

fn diag_dominant(n: usize, p: usize) -> Option<ExactNumArray> {
    let mut vals = Vec::with_capacity(n.saturating_mul(n));
    for i in 0..n {
        for j in 0..n {
            vals.push(if i == j {
                ExactNum::from_u32((n as u32).saturating_add(3), p)
            } else {
                ExactNum::from_u8(1, p)
            });
        }
    }
    ExactNumArray::from_shape(p, n, n, &vals)
}

fn bench_matmul(c: &mut Criterion) {
    let mut group = c.benchmark_group("linalg/matmul");
    let p = 64;
    let one = ExactNum::from_u8(1, p);
    for &n in &[10usize, 100] {
        let Some(a) = ExactNumArray::filled_2d(p, n, n, &one) else {
            continue;
        };
        group.bench_with_input(BenchmarkId::from_parameter(n), &a, |b, a| {
            b.iter(|| sink(a.matmul(a)));
        });
    }
    let one64 = Ieee64::from_i32(1);
    let Some(big) = Ieee64Array::filled_2d(1000, 1000, one64) else {
        group.finish();
        return;
    };
    group.sample_size(10);
    group.bench_with_input(BenchmarkId::from_parameter(1000), &big, |b, a| {
        b.iter(|| sink(a.matmul(a)));
    });
    group.finish();
}

fn bench_fft(c: &mut Criterion) {
    let mut group = c.benchmark_group("linalg/fft");
    let p = 64;
    let rm = RoundingMode::ToEven;
    let one = ExactNum::from_u8(1, p);
    for &n in &[256usize, 1024, 4096] {
        let a = ExactNumArray::filled(p, n, &one);
        group.bench_with_input(BenchmarkId::from_parameter(n), &a, |b, a| {
            let Ok(mut cc) = Consts::new() else {
                return;
            };
            b.iter(|| sink(a.fft(p, rm, &mut cc)));
        });
    }
    group.finish();
}

fn bench_lu(c: &mut Criterion) {
    let mut group = c.benchmark_group("linalg/lu");
    let p = 64;
    let rm = RoundingMode::ToEven;
    for &n in &[50usize, 200] {
        let Some(a) = diag_dominant(n, p) else {
            continue;
        };
        group.bench_with_input(BenchmarkId::from_parameter(n), &a, |b, a| {
            b.iter(|| sink(a.lu_decomp(p, rm)));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_matmul, bench_fft, bench_lu);
criterion_main!(benches);
