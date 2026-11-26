//! Fold PR 1 (RELEASE-2.0-CHECKLIST 4.1; fold plan row L): the facade keeps `zenith_float::Ball` as the
//! zenith-float-num ball (ExactNum midpoint and radius). The bound-layer ball (ExactNum midpoint, `Mag` radius)
//! is reachable only as `zenith_float::bound::ball::Ball`, behind the default-on `bound` feature, never at the
//! facade root.
use core::any::TypeId;

#[test]
fn fold_pr1_facade_ball_still_num_ball() {
    assert_eq!(
        TypeId::of::<zenith_float::Ball>(),
        TypeId::of::<zenith_float_num::Ball>()
    );
}

#[cfg(feature = "bound")]
#[test]
fn fold_pr1_bound_ball_is_bound_crate_ball() {
    assert_eq!(
        TypeId::of::<zenith_float::bound::ball::Ball>(),
        TypeId::of::<zenith_float_bound::ball::Ball>()
    );
    assert_ne!(
        TypeId::of::<zenith_float::bound::ball::Ball>(),
        TypeId::of::<zenith_float::Ball>()
    );
}
