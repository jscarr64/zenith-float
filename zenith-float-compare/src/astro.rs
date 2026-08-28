//! astro-float backend (BigFloat).

use std::cell::RefCell;
use std::rc::Rc;
use std::time::Instant;

use astro_float::{BigFloat, Consts, Exponent, RoundingMode, Sign};

use crate::backend::{convert_exp, BenchFloat, BenchState};
use crate::task::{task_for_one_arg, task_for_two_args};

#[derive(Clone)]
pub struct AstroState {
    pub precision: usize,
    pub cc: Rc<RefCell<Consts>>,
}

impl BenchState for AstroState {}

#[derive(Clone)]
pub struct AstroFloat {
    inner: BigFloat,
    cc: Rc<RefCell<Consts>>,
    precision: usize,
}

impl AstroFloat {
    fn new(inner: BigFloat, cc: Rc<RefCell<Consts>>, precision: usize) -> Self {
        Self {
            inner,
            cc,
            precision,
        }
    }

    fn p(&self) -> usize {
        self.precision
    }

    fn rm(&self) -> RoundingMode {
        RoundingMode::ToEven
    }

    fn cc_mut(&self) -> std::cell::RefMut<'_, Consts> {
        self.cc.borrow_mut()
    }

    fn inner(&self) -> &BigFloat {
        &self.inner
    }
}

impl BenchFloat for AstroFloat {
    type State = AstroState;

    fn state(precision: usize) -> Self::State {
        AstroState {
            precision,
            cc: Rc::new(RefCell::new(Consts::new().expect("constants cache"))),
        }
    }

    fn rand_normal(
        count: usize,
        exp_from: i32,
        exp_to: i32,
        state: &Self::State,
        sign_positive: bool,
    ) -> Vec<Self> {
        let exp_from = convert_exp(exp_from) as Exponent;
        let exp_to = convert_exp(exp_to) as Exponent;
        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            let mut inner = BigFloat::random_normal(state.precision, exp_from, exp_to);
            if sign_positive {
                inner.set_sign(Sign::Pos);
            }
            out.push(Self::new(inner, state.cc.clone(), state.precision));
        }
        out
    }

    fn run_task(state: &Self::State, task: &str, values: &[Self]) -> (Self, std::time::Duration) {
        let start = Instant::now();
        let p = state.precision;
        let cc = state.cc.clone();
        let result = match task {
            "add" => task_for_two_args(values, |a, b| {
                Self::new(a.inner().add(b.inner(), a.p(), a.rm()), cc.clone(), p)
            }),
            "sub" => task_for_two_args(values, |a, b| {
                Self::new(a.inner().sub(b.inner(), a.p(), a.rm()), cc.clone(), p)
            }),
            "mul" => task_for_two_args(values, |a, b| {
                Self::new(a.inner().mul(b.inner(), a.p(), a.rm()), cc.clone(), p)
            }),
            "div" => task_for_two_args(values, |a, b| {
                Self::new(a.inner().div(b.inner(), a.p(), a.rm()), cc.clone(), p)
            }),
            "sqrt" => task_for_one_arg(values, |a, b| {
                Self::new(b.inner().sqrt(a.p(), a.rm()), cc.clone(), p)
            }),
            "cbrt" => task_for_one_arg(values, |a, b| {
                Self::new(b.inner().cbrt(a.p(), a.rm()), cc.clone(), p)
            }),
            "ln" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().ln(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "exp" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().exp(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "pow" => task_for_two_args(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(
                    a.inner().pow(b.inner(), a.p(), a.rm(), &mut c),
                    cc.clone(),
                    p,
                )
            }),
            "sin" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().sin(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "asin" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().asin(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "cos" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().cos(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "acos" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().acos(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "tan" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().tan(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "atan" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().atan(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "sinh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().sinh(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "asinh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().asinh(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "cosh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().cosh(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "acosh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().acosh(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "tanh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().tanh(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            "atanh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().atanh(a.p(), a.rm(), &mut c), cc.clone(), p)
            }),
            _ => panic!("unknown task: {task}"),
        };
        (result, start.elapsed())
    }
}
