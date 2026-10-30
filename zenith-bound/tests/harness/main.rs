//! The ONE harness integration test target of zenith-bound: shared support (tests/common), its self-tests and, later, the suites built on it.
//! tests/common is included here and nowhere else (see the rule in tests/common/mod.rs).
#[path = "../common/mod.rs"]
mod common;

mod ball_boundary;
mod ball_goldens_selftest;
mod edges_selftest;
mod enclosure_selftest;
mod fuzz_mag;
mod fuzz_mag_l0a5;
mod fuzz_mag_props;
mod fuzz_mag_props_selftest;
mod generators_selftest;
mod goldens_selftest;
mod grid_digests;
mod mag_carry;
mod mag_rows;
mod mag_sort;
mod modular_selftest;
mod points_selftest;
mod replay_selftest;
mod rng_selftest;
mod seeds_selftest;
mod shrink_selftest;
mod tightness_selftest;
mod w1_translation;

// Per-slice harness files (S0, 2026-10-24): each belongs to its slice alone; this list is final.
mod slice_arith;
mod slice_bridge;
mod slice_carlson;
mod slice_cball;
mod slice_cball_s1_cover;
mod slice_cholesky;
mod slice_eigen;
mod slice_elementary;
mod slice_gls;
mod slice_hypergeom;
mod slice_lu;
mod slice_matrix;
mod slice_qr;
mod slice_quadrature;
mod slice_zf_special;
