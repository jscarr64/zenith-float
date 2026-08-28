//! zenith-float backend (ExactNum).

use std::cell::RefCell;
use std::rc::Rc;
use std::time::Instant;

use zenith_float::{Consts, ExactNum, Exponent, RoundingMode, Sign};

use crate::backend::{convert_exp, BenchFloat, BenchState};
use crate::task::{task_for_one_arg, task_for_two_args};

#[derive(Clone)]
pub struct ZenithState {
    pub precision: usize,
    pub cc: Rc<RefCell<Consts>>,
}

impl BenchState for ZenithState {}

#[derive(Clone)]
pub struct ZenithFloat {
    inner: ExactNum,
    cc: Rc<RefCell<Consts>>,
    precision: usize,
}

impl ZenithFloat {
    fn new(inner: ExactNum, cc: Rc<RefCell<Consts>>, precision: usize) -> Self {
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
}

impl BenchFloat for ZenithFloat {
    type State = ZenithState;

    fn state(precision: usize) -> Self::State {
        ZenithState {
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
            let mut v = ExactNum::random_normal(state.precision, exp_from, exp_to);
            if sign_positive {
                v.set_sign(Sign::Pos);
            }
            out.push(Self::new(v, state.cc.clone(), state.precision));
        }
        out
    }

    fn run_task(state: &Self::State, task: &str, values: &[Self]) -> (Self, std::time::Duration) {
        let start = Instant::now();
        let cc = state.cc.clone();
        let result = match task {
            "add" => task_for_two_args(values, |a, b| {
                Self::new(a.inner.add(b.inner(), a.p(), a.rm()), cc.clone(), a.p())
            }),
            "sub" => task_for_two_args(values, |a, b| {
                Self::new(a.inner.sub(b.inner(), a.p(), a.rm()), cc.clone(), a.p())
            }),
            "mul" => task_for_two_args(values, |a, b| {
                Self::new(a.inner.mul(b.inner(), a.p(), a.rm()), cc.clone(), a.p())
            }),
            "div" => task_for_two_args(values, |a, b| {
                Self::new(a.inner.div(b.inner(), a.p(), a.rm()), cc.clone(), a.p())
            }),
            "sqrt" => task_for_one_arg(values, |a, b| {
                Self::new(b.inner().sqrt(a.p(), a.rm()), cc.clone(), a.p())
            }),
            "cbrt" => task_for_one_arg(values, |a, b| {
                Self::new(b.inner().cbrt(a.p(), a.rm()), cc.clone(), a.p())
            }),
            "ln" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().ln(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "exp" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().exp(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "pow" => task_for_two_args(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(
                    a.inner().pow(b.inner(), a.p(), a.rm(), &mut c),
                    cc.clone(),
                    a.p(),
                )
            }),
            "sin" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().sin(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "asin" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().asin(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "cos" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().cos(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "acos" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().acos(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "tan" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().tan(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "atan" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().atan(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "sinh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().sinh(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "asinh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().asinh(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "cosh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().cosh(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "acosh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().acosh(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "tanh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().tanh(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            "atanh" => task_for_one_arg(values, |a, b| {
                let mut c = a.cc_mut();
                Self::new(b.inner().atanh(a.p(), a.rm(), &mut c), cc.clone(), a.p())
            }),
            _ => panic!("unknown task: {task}"),
        };
        (result, start.elapsed())
    }
}

impl ZenithFloat {
    fn inner(&self) -> &ExactNum {
        &self.inner
    }
}
