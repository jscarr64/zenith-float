//! Cross-library comparison runner (bigfloat-bench compatible workloads, TSV output).

mod astro;
mod backend;
#[cfg(feature = "dashu")]
mod dashu;
mod task;
mod zenith;

use clap::Parser;

use crate::backend::BenchFloat;
use crate::task::{workload, TASKS};

#[derive(Parser, Debug)]
#[command(
    name = "zenith-float-compare",
    about = "Apples-to-apples big-float benchmarks vs astro-float and dashu-float"
)]
struct Args {
    /// Libraries to benchmark.
    #[arg(long = "lib", value_name = "NAME", required = true)]
    libraries: Vec<String>,

    /// Tasks (add, mul, ln, sin, …).
    #[arg(long = "task", value_name = "TASK", required = true)]
    tasks: Vec<String>,

    /// Binary precision in bits (repeat for multiple tiers).
    #[arg(long = "precision", value_name = "BITS", default_values_t = vec![132, 1000, 10000])]
    precisions: Vec<usize>,

    /// Best-of-N runs (bigfloat-bench default is 5).
    #[arg(short = 'n', default_value_t = 5)]
    runs: usize,

    /// Emit TSV to stdout (default) or a human table to stderr summary only.
    #[arg(long)]
    table: bool,
}

fn main() {
    let args = Args::parse();
    for task in &args.tasks {
        if !TASKS.contains(&task.as_str()) {
            eprintln!("error: unknown task {task:?} (expected one of: {TASKS:?})");
            std::process::exit(1);
        }
    }
    if !args.table {
        println!("library\tprecision_bits\ttask\tbatch_size\tbest_time_us");
    }

    for &precision in &args.precisions {
        for task in &args.tasks {
            let (batch, exp_from, exp_to, sign_positive) = workload(task);
            for lib in &args.libraries {
                let best_us = match lib.as_str() {
                    "zenith" => bench_lib::<zenith::ZenithFloat>(precision, task, batch, exp_from, exp_to, sign_positive, args.runs),
                    #[cfg(feature = "astro")]
                    "astro" | "astro-float" => bench_lib::<astro::AstroFloat>(precision, task, batch, exp_from, exp_to, sign_positive, args.runs),
                    #[cfg(feature = "dashu")]
                    "dashu" | "dashu-float" => bench_lib::<dashu::DashuFloat>(precision, task, batch, exp_from, exp_to, sign_positive, args.runs),
                    other => {
                        eprintln!("error: unknown library {other:?} (enable features: astro, dashu)");
                        std::process::exit(1);
                    }
                };
                if args.table {
                    eprint!("{lib:>15} p={precision:<6} {task:>6} ");
                    eprintln!("{best_us} us (batch {batch})");
                } else {
                    println!("{lib}\t{precision}\t{task}\t{batch}\t{best_us}");
                }
            }
        }
    }
}

fn bench_lib<F: BenchFloat>(
    precision: usize,
    task: &str,
    batch: usize,
    exp_from: i32,
    exp_to: i32,
    sign_positive: bool,
    runs: usize,
) -> u64 {
    let state = F::state(precision);
    let values = F::rand_normal(batch, exp_from, exp_to, &state, sign_positive);
    let mut samples = Vec::with_capacity(runs);
    for _ in 0..runs {
        samples.push(measure_batch::<F>(&state, task, &values));
    }
    samples.sort_unstable();
    samples[0]
}

/// Adaptive timing loop from bigfloat-bench: accumulate until ≥1 ms, then report µs/run.
fn measure_batch<F: BenchFloat>(state: &F::State, task: &str, values: &[F]) -> u64 {
    let mut full_us = 0u128;
    let mut iter = 1usize;
    let mut niter = 0usize;
    while full_us < 1_000 && iter < 16 {
        niter += iter;
        for _ in 0..iter {
            let (_, d) = F::run_task(state, task, values);
            full_us += d.as_micros();
        }
        iter *= 2;
    }
    (full_us / niter as u128) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn workload_sizes_match_bigfloat_bench() {
        assert_eq!(workload("mul"), (1_000_000, -10, 10, false));
        assert_eq!(workload("ln"), (10_000, -10, 10, true));
    }
}
