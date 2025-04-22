//! Integration tests for constant depth

use std::f64::consts::PI;

use crate::bathymetry::ConstantDepth;
use crate::current::ConstantCurrent;
use crate::datatype::{Point, RayState, WaveNumber};
use crate::ray::ManyRays;

// import the helper functions and types for integration testing
use crate::tests::helper::*;

#[test]
/// rays in deep water with constant depth
///
/// ## Initial conditions:
///
/// 12 rays where i = 0, 1, 2, ..., 11
///
/// k = 0.05 m^-1
///
/// kx0 = k * cos(pi*i/6)
///
/// ky0 = k * sin(pi*i/6)
///
/// h = 2000 m
///
/// (u, v) = (0, 0)
///
/// ## Description:
///
/// Rays propagate from left to right over constant bathymetry in deep water.
/// Initial direction increments of 30 degrees.
///
/// ## Expected behavior:
/// The ray path goes straight from the center and the values of (kx, ky) are
/// equal to the initial condition at all points.
fn constant_depth_deep() {
    // load the data
    let bathymetry_data = ConstantDepth::new(2000.0);
    let current_data = ConstantCurrent::new(0.0, 0.0); // default (u, v) = (0, 0)

    // create 12 rays starting at the same point and in a circle with angle pi/6 between them
    let init_rays: Vec<RayState<f64>> = (0..12)
        .map(|i| {
            RayState::new(
                Point::new(50_000.0, 25_000.0),
                WaveNumber::new(
                    0.05 * (PI * i as f64 / 6.0).cos(),
                    0.05 * (PI * i as f64 / 6.0).sin(),
                ),
            )
        })
        .collect();

    let rays = ManyRays::new(&bathymetry_data, &current_data, &init_rays);

    let results = rays.trace_many(0.0, 5000.0, 1.0);

    for (i, ray) in results.iter().flatten().enumerate() {
        let (_, data) = &ray.get();
        let target_ray = init_rays.get(i).unwrap().to_owned();
        let kx = target_ray.wave_number().kx();
        let ky = target_ray.wave_number().ky();

        // x
        if (kx - 0.0).abs() < f64::EPSILON {
            assert!(same(data, XINDEX));
        } else if kx.is_sign_positive() {
            assert!(increase(data, XINDEX));
        } else {
            assert!(decrease(data, XINDEX));
        }

        // y
        if (ky - 0.0).abs() < f64::EPSILON {
            assert!(same(data, YINDEX));
        } else if ky.is_sign_positive() {
            assert!(increase(data, YINDEX));
        } else {
            assert!(decrease(data, YINDEX));
        }

        // kx and ky will be the same
        assert!(same(data, KX_INDEX));
        assert!(same(data, KY_INDEX));
    }
}

#[test]
/// rays in shallow water with constant depth
///
/// ## Initial conditions:
///
/// 12 rays where i = 0, 1, 2, ..., 11
///
/// k = 0.05 m^-1
///
/// kx0 = k * cos(pi*i/6)
///
/// ky0 = k * sin(pi*i/6)
///
/// h = 10 m
///
/// (u, v) = (0, 0)
///
/// ## Description:
/// Rays propagate from the center over constant bathymetry in shallow water.
/// Initial direction increments of 30 degrees.
///
/// ## Expected behavior:
/// The ray path goes straight from the center and the values of (kx, ky) are
/// equal to the initial condition at all points.
fn constant_depth_shallow() {
    let bathymetry_data = ConstantDepth::new(10.0);
    let current_data = ConstantCurrent::new(0.0, 0.0); // default (u, v) = (0, 0)

    // create 12 rays starting at the same point and in a circle with angle pi/6 between them
    let init_rays: Vec<RayState<f64>> = (0..12)
        .map(|i| {
            RayState::new(
                Point::new(50_000.0, 25_000.0),
                WaveNumber::new(
                    0.05 * (PI * i as f64 / 6.0).cos(),
                    0.05 * (PI * i as f64 / 6.0).sin(),
                ),
            )
        })
        .collect();

    let rays = ManyRays::new(&bathymetry_data, &current_data, &init_rays);

    let results = rays.trace_many(0.0, 5000.0, 1.0);

    for (i, ray) in results.iter().flatten().enumerate() {
        let (_, data) = &ray.get();
        let target_ray = init_rays.get(i).unwrap().to_owned();
        let kx = target_ray.wave_number().kx();
        let ky = target_ray.wave_number().ky();
        // x
        if (kx - 0.0).abs() < f64::EPSILON {
            assert!(same(data, XINDEX));
        } else if kx.is_sign_positive() {
            assert!(increase(data, XINDEX));
        } else {
            assert!(decrease(data, XINDEX));
        }

        // y
        if (ky - 0.0).abs() < f64::EPSILON {
            assert!(same(data, YINDEX));
        } else if ky.is_sign_positive() {
            assert!(increase(data, YINDEX));
        } else {
            assert!(decrease(data, YINDEX));
        }

        // kx and ky will be the same
        assert!(same(data, KX_INDEX));
        assert!(same(data, KY_INDEX));
    }
}
