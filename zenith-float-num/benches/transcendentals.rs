//! Transcendental function benchmarks with fixed, domain-valid operands.

mod shared;
mod trans_common;

use shared::{cycle_batch, init_cc, parse_list, sink, PRECISIONS};
use trans_common::{
    EXP_FIXTURES, LN_FIXTURES, TRIG_FIXTURES, TRANSCENDENTAL_BATCH,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use zenith_float_num::{Consts, ExactNum, RoundingMode};

fn bench_ln(c: &mut Criterion) {
    let mut group = c.benchmark_group("transcendentals/ln");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let mut cc = init_cc();
        let values = cycle_batch(&parse_list(p, &mut cc, LN_FIXTURES), TRANSCENDENTAL_BATCH);
        group.throughput(Throughput::Elements(TRANSCENDENTAL_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &values, |b, values| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                for x in values {
                    sink(x.ln(p, rm, &mut cc));
                }
            });
        });
    }
    group.finish();
}

fn bench_exp(c: &mut Criterion) {
    let mut group = c.benchmark_group("transcendentals/exp");
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let mut cc = init_cc();
        let values = cycle_batch(&parse_list(p, &mut cc, EXP_FIXTURES), TRANSCENDENTAL_BATCH);
        group.throughput(Throughput::Elements(TRANSCENDENTAL_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &values, |b, values| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                for x in values {
                    sink(x.exp(p, rm, &mut cc));
                }
            });
        });
    }
    group.finish();
}

fn bench_sin(c: &mut Criterion) {
    bench_trig_unary(c, "sin", ExactNum::sin);
}

fn bench_cos(c: &mut Criterion) {
    bench_trig_unary(c, "cos", ExactNum::cos);
}

fn bench_atan(c: &mut Criterion) {
    bench_trig_unary(c, "atan", ExactNum::atan);
}

fn bench_sinh(c: &mut Criterion) {
    bench_trig_unary(c, "sinh", ExactNum::sinh);
}

fn bench_cosh(c: &mut Criterion) {
    bench_trig_unary(c, "cosh", ExactNum::cosh);
}

fn bench_tanh(c: &mut Criterion) {
    bench_trig_unary(c, "tanh", ExactNum::tanh);
}

fn bench_trig_unary<F>(c: &mut Criterion, name: &str, op: F)
where
    F: Copy + Fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let mut group = c.benchmark_group(format!("transcendentals/{name}"));
    let rm = RoundingMode::ToEven;

    for &p in PRECISIONS {
        let mut cc = init_cc();
        let values = cycle_batch(&parse_list(p, &mut cc, TRIG_FIXTURES), TRANSCENDENTAL_BATCH);
        group.throughput(Throughput::Elements(TRANSCENDENTAL_BATCH as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &values, |b, values| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                for x in values {
                    sink(op(x, p, rm, &mut cc));
                }
            });
        });
    }
    group.finish();
}

fn bench_pow(c: &mut Criterion) {
    let mut group = c.benchmark_group("transcendentals/pow");
    let rm = RoundingMode::ToEven;
    let batch = TRANSCENDENTAL_BATCH;

    // y^x with fixed bases and exponents from real constants.
    let bases: &[&str] = &["2.718281828459045", "3.141592653589793", "1.414213562373095"];
    let exponents: &[&str] = &["0.5", "1.5", "2.0", "3.0", "0.25"];

    for &p in &[128, 1024, 4096] {
        let mut cc = init_cc();
        let mut pairs = Vec::with_capacity(batch);
        let parsed_bases = parse_list(p, &mut cc, bases);
        let parsed_exps = parse_list(p, &mut cc, exponents);
        for i in 0..batch {
            pairs.push((
                parsed_bases[i % parsed_bases.len()].clone(),
                parsed_exps[i % parsed_exps.len()].clone(),
            ));
        }

        group.throughput(Throughput::Elements(batch as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                for (base, exp) in pairs {
                    sink(base.pow(exp, p, rm, &mut cc));
                }
            });
        });
    }
    group.finish();
}

fn bench_hypot(c: &mut Criterion) {
    let mut group = c.benchmark_group("transcendentals/hypot");
    let rm = RoundingMode::ToEven;
    let batch = TRANSCENDENTAL_BATCH;

    let legs: &[&str] = &["3.0", "4.0", "5.0", "12.0", "1.234567890123456789"];

    for &p in PRECISIONS {
        let mut cc = init_cc();
        let parsed = parse_list(p, &mut cc, legs);
        let n = parsed.len();
        let pairs: Vec<_> = (0..batch)
            .map(|i| {
                (
                    parsed[i % n].clone(),
                    parsed[(i.wrapping_mul(5).wrapping_add(2)) % n].clone(),
                )
            })
            .collect();

        group.throughput(Throughput::Elements(batch as u64));
        group.bench_with_input(BenchmarkId::from_parameter(p), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_val) in pairs {
                    sink(a.hypot(b_val, p, rm));
                }
            });
        });
    }
    group.finish();
}

fn bench_log2_log10(c: &mut Criterion) {
    let mut group = c.benchmark_group("transcendentals/log2_log10");
    let rm = RoundingMode::ToEven;

    for &p in &[128, 1024, 4096] {
        let mut cc = init_cc();
        let values = cycle_batch(&parse_list(p, &mut cc, LN_FIXTURES), TRANSCENDENTAL_BATCH);
        group.throughput(Throughput::Elements(TRANSCENDENTAL_BATCH as u64));
        group.bench_with_input(BenchmarkId::new("log2", p), &values, |b, values| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                for x in values {
                    sink(x.log2(p, rm, &mut cc));
                }
            });
        });
        group.bench_with_input(BenchmarkId::new("log10", p), &values, |b, values| {
            let mut cc = Consts::new().expect("constants cache");
            b.iter(|| {
                for x in values {
                    sink(x.log10(p, rm, &mut cc));
                }
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_ln,
    bench_exp,
    bench_sin,
    bench_cos,
    bench_atan,
    bench_sinh,
    bench_cosh,
    bench_tanh,
    bench_pow,
    bench_hypot,
    bench_log2_log10,
);
criterion_main!(benches);
