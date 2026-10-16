//! Arithmetic benchmarks: add, mul, div, reciprocal, remainder, sqrt.

mod common;
mod shared;

use common::{pair_fixtures, ARITH_BATCH, DECIMAL_FIXTURES, FFT_MUL_BATCH, FFT_MUL_PRECISIONS};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use shared::{cycle_batch, init_cc, sink, PRECISIONS};
use zenith_float_num::{ExactNum, Radix, RoundingMode};

fn bench_add(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/add");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let pairs = pair_fixtures(p, &mut cc, ARITH_BATCH);
        group.throughput(Throughput::Elements(ARITH_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_val) in pairs {
                    sink(ExactNum::add_ref(a, b_val, p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/mul");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let pairs = pair_fixtures(p, &mut cc, ARITH_BATCH);
        group.throughput(Throughput::Elements(ARITH_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_val) in pairs {
                    sink(ExactNum::mul_ref(a, b_val, p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_mul_fft_scale(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/mul_fft");
    group.sample_size(10);
    let rm = RoundingMode::ToEven;

    for &p in FFT_MUL_PRECISIONS {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let pairs = pair_fixtures(p, &mut cc, FFT_MUL_BATCH);
        group.throughput(Throughput::Elements(FFT_MUL_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_val) in pairs {
                    sink(ExactNum::mul_ref(a, b_val, p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_div(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/div");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let pairs = pair_fixtures(p, &mut cc, ARITH_BATCH);
        group.throughput(Throughput::Elements(ARITH_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_val) in pairs {
                    sink(ExactNum::div_ref(a, b_val, p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_reciprocal(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/reciprocal");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let values = cycle_batch(p, &mut cc, DECIMAL_FIXTURES, ARITH_BATCH);
        group.throughput(Throughput::Elements(ARITH_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &values, |b, values| {
            b.iter(|| {
                for x in values {
                    sink(x.reciprocal(p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_recip_mul_vs_div(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/recip_mul_vs_div");
    let rm = RoundingMode::ToEven;
    // Division at very high precision is a common hot path.
    let p = 14400;

    let Some(mut cc) = init_cc() else {
        return;
    };
    let values = cycle_batch(p, &mut cc, DECIMAL_FIXTURES, 64);
    let divisor = ExactNum::parse(
        DECIMAL_FIXTURES[0],
        Radix::Dec,
        p,
        RoundingMode::ToEven,
        &mut cc,
    );

    group.throughput(Throughput::Elements(63));
    group.bench_function("recip_then_mul", |b| {
        b.iter(|| {
            let inv = divisor.reciprocal(p, RoundingMode::None);
            for x in values.iter().skip(1) {
                sink(ExactNum::mul_ref(&inv, x, p, RoundingMode::None));
            }
        });
    });
    group.bench_function("direct_div", |b| {
        b.iter(|| {
            for x in values.iter().skip(1) {
                sink(ExactNum::div_ref(x, &divisor, p, rm));
            }
        });
    });
    group.finish();
}

fn bench_rem(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/rem");
    let batch = 64;

    for &p in &[1024, 4096, 32768] {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let pairs = pair_fixtures(p, &mut cc, batch);
        group.throughput(Throughput::Elements(batch as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_val) in pairs {
                    sink(a.rem(b_val));
                }
            });
        });
    }
    group.finish();
}

fn bench_sqrt(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/sqrt");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let values = cycle_batch(p, &mut cc, DECIMAL_FIXTURES, ARITH_BATCH);
        group.throughput(Throughput::Elements(ARITH_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &values, |b, values| {
            b.iter(|| {
                for x in values {
                    sink(x.sqrt(p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_div_one(c: &mut Criterion) {
    let mut group = c.benchmark_group("arithmetic/div_one");
    let rm = RoundingMode::ToEven;
    let batch = 64;

    for &p in &[4096, 14400, 32768] {
        let Some(mut cc) = init_cc() else {
            return;
        };
        let values = cycle_batch(p, &mut cc, DECIMAL_FIXTURES, batch);
        group.throughput(Throughput::Elements(batch as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &values, |b, values| {
            let one = ExactNum::from_word(1, p);
            b.iter(|| {
                for x in values {
                    sink(ExactNum::div_ref(&one, x, p, rm));
                }
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_add,
    bench_mul,
    bench_mul_fft_scale,
    bench_div,
    bench_reciprocal,
    bench_recip_mul_vs_div,
    bench_rem,
    bench_sqrt,
    bench_div_one,
);
criterion_main!(benches);
