//! dashu-float backend (`FBig`, binary significand).

#![cfg(feature = "dashu")]

use std::time::Instant;

use dashu_float::ops::CubicRoot;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::{IBig, Sign, UBig};
use rand::random;

use crate::backend::{convert_exp, BenchFloat, BenchState};
use crate::task::{task_for_one_arg, task_for_two_args};

type FBin = FBig<HalfEven, 2>;

#[derive(Clone, Copy)]
pub struct DashuState {
    pub precision: usize,
}

impl BenchState for DashuState {}

#[derive(Clone)]
pub struct DashuFloat {
    inner: FBin,
}

impl DashuFloat {
    fn new(inner: FBin) -> Self {
        Self { inner }
    }

    fn inner(&self) -> &FBin {
        &self.inner
    }
}

fn word_count(precision: usize) -> usize {
    precision.div_ceil(64).max(1)
}

impl BenchFloat for DashuFloat {
    type State = DashuState;

    fn state(precision: usize) -> Self::State {
        DashuState { precision }
    }

    fn rand_normal(
        count: usize,
        exp_from: i32,
        exp_to: i32,
        state: &Self::State,
        sign_positive: bool,
    ) -> Vec<Self> {
        let words = word_count(state.precision);
        let exp_from = convert_exp(exp_from);
        let exp_to = convert_exp(exp_to);
        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            let mut mantissa = vec![0u64; words];
            for w in &mut mantissa {
                *w = random();
            }
            let last = mantissa.len() - 1;
            if mantissa[last] == 0 {
                mantissa[last] = u64::MAX;
            }
            while mantissa[last] <= (u64::MAX >> 1) {
                mantissa[last] <<= 1;
            }
            let exp_range = exp_to - exp_from;
            let sign = if sign_positive || random::<u8>() & 1 == 0 {
                Sign::Positive
            } else {
                Sign::Negative
            };
            let exp = (if exp_range != 0 { random::<i32>().abs() % exp_range } else { 0 })
                - words as i32 * 64
                + exp_from;
            let exp = exp as isize * 3_321_928_095 / 1_000_000_000;
            let m = UBig::from_words(&mantissa);
            let i = IBig::from_parts(sign, m);
            out.push(Self::new(FBin::from_parts(i, exp)));
        }
        out
    }

    fn run_task(_state: &Self::State, task: &str, values: &[Self]) -> (Self, std::time::Duration) {
        let start = Instant::now();
        let result = match task {
            "add" => task_for_two_args(values, |a, b| Self::new(a.inner() + b.inner())),
            "sub" => task_for_two_args(values, |a, b| Self::new(a.inner() - b.inner())),
            "mul" => task_for_two_args(values, |a, b| Self::new(a.inner() * b.inner())),
            "div" => task_for_two_args(values, |a, b| Self::new(a.inner() / b.inner())),
            "sqrt" => task_for_one_arg(values, |_a, b| Self::new(b.inner().sqrt())),
            "cbrt" => task_for_one_arg(values, |_a, b| Self::new(b.inner().cbrt())),
            "ln" => task_for_one_arg(values, |_a, b| Self::new(FBin::ln(b.inner()))),
            "exp" => task_for_one_arg(values, |_a, b| Self::new(FBin::exp(b.inner()))),
            "pow" => task_for_two_args(values, |a, b| Self::new(FBin::powf(a.inner(), b.inner()))),
            "sin" => task_for_one_arg(values, |_a, b| Self::new(FBin::sin(b.inner()))),
            "asin" => task_for_one_arg(values, |_a, b| Self::new(FBin::asin(b.inner()))),
            "cos" => task_for_one_arg(values, |_a, b| Self::new(FBin::cos(b.inner()))),
            "acos" => task_for_one_arg(values, |_a, b| Self::new(FBin::acos(b.inner()))),
            "tan" => task_for_one_arg(values, |_a, b| Self::new(FBin::tan(b.inner()))),
            "atan" => task_for_one_arg(values, |_a, b| Self::new(FBin::atan(b.inner()))),
            "sinh" => task_for_one_arg(values, |_a, b| Self::new(FBin::sinh(b.inner()))),
            "asinh" => task_for_one_arg(values, |_a, b| Self::new(FBin::asinh(b.inner()))),
            "cosh" => task_for_one_arg(values, |_a, b| Self::new(FBin::cosh(b.inner()))),
            "acosh" => task_for_one_arg(values, |_a, b| Self::new(FBin::acosh(b.inner()))),
            "tanh" => task_for_one_arg(values, |_a, b| Self::new(FBin::tanh(b.inner()))),
            "atanh" => task_for_one_arg(values, |_a, b| Self::new(FBin::atanh(b.inner()))),
            _ => panic!("unknown task: {task}"),
        };
        (result, start.elapsed())
    }
}
