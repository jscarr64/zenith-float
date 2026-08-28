//! Workload sizes and reduction loops (aligned with bigfloat-bench).

pub const TASKS: &[&str] = &[
    "add", "sub", "mul", "div", "sqrt", "cbrt", "ln", "exp", "pow", "sin", "asin", "cos", "acos",
    "tan", "atan", "sinh", "asinh", "cosh", "acosh", "tanh", "atanh",
];

/// `(batch_size, exp_from_base10, exp_to_base10, force_positive)`.
pub fn workload(task: &str) -> (usize, i32, i32, bool) {
    match task {
        "add" | "sub" | "mul" | "div" => (1_000_000, -10, 10, false),
        "sqrt" => (100_000, -10, 10, true),
        "cbrt" => (100_000, -10, 10, false),
        "ln" => (10_000, -10, 10, true),
        "exp" => (10_000, -10, 3, false),
        "pow" => (10_000, -5, 5, false),
        "sin" | "cos" | "tan" | "sinh" | "cosh" | "tanh" => (10_000, -10, 3, false),
        "asin" | "acos" | "atan" | "atanh" => (10_000, -10, 0, false),
        "asinh" => (10_000, -10, 10, false),
        "acosh" => (10_000, 1, 10, true),
        _ => panic!("unknown task: {task}"),
    }
}

pub fn task_for_two_args<T, F>(values: &[T], mut op: F) -> T
where
    T: Clone,
    F: FnMut(T, &T) -> T,
{
    let mut acc = values[0].clone();
    let (left, right) = values.split_at(values.len() / 2);
    for _ in 0..2 {
        for (_a, b) in left.iter().zip(right) {
            acc = op(acc, b);
        }
    }
    acc
}

pub fn task_for_one_arg<T, F>(values: &[T], mut op: F) -> T
where
    T: Clone,
    F: FnMut(T, &T) -> T,
{
    let mut acc = values[0].clone();
    for v in values {
        acc = op(acc, v);
    }
    acc
}
