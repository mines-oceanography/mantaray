//! integration tests for linear beaches

use std::f64::consts::PI;

use crate::{
    bathymetry::ConstantSlope,
    current::ConstantCurrent,
    datatype::{Point, RayState, WaveNumber},
    ray::ManyRays,
};

use crate::tests::helper::*;

#[test]
/// test a linear beach on the right side of domain starting in shallow water
///
/// ## Bathymetry
/// `ConstantSlope` object with initial conditions:
/// - `x0 = 0 m`
/// - `y0 = 0 m`
/// - `h0 = 100 m`
/// - `dhdx =  -0.05`
/// - `dhdy = 0`
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 0 m`
/// - `y = 0 m`
/// - `kx = k * cos(PI/6)`
/// - `ky = k * sin(PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 0 m`
/// - `y = 0 m`
/// - `kx = k * cos(-PI/6)`
/// - `ky = k * sin(-PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 0 m`
/// - `y = 0 m`
/// - `kx = k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from left to right starting in shallow water. The
/// rays start on a beach and immediately start to curve towards the beach,
/// so the kx values of each ray will increase.
///
/// ## Expected behavior
/// The rays will curve towards the beach
fn test_linear_beach_right() {
    let bathymetry_data = ConstantSlope::builder()
        .x0(0.0)
        .y0(0.0)
        .h0(100.0)
        .dhdx(-0.05)
        .dhdy(0.0)
        .build()
        .unwrap();

    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let up_ray = RayState::new(
        Point::new(0.0, 0.0),
        WaveNumber::new(k * (PI / 6.0).cos(), k * (PI / 6.0).sin()),
    );

    let down_ray = RayState::new(
        Point::new(0.0, 0.0),
        WaveNumber::new(k * (-PI / 6.0).cos(), k * (-PI / 6.0).sin()),
    );

    let straight_ray = RayState::new(Point::new(0.0, 0.0), WaveNumber::new(k, 0.0));

    let initial_rays = vec![up_ray, down_ray, straight_ray];

    let waves = ManyRays::new(&bathymetry_data, &current_data, &initial_rays);

    let results = waves.trace_many(0.0, 1_000.0, 1.0);

    let mut results_iter = results.iter().flatten();

    // order is up, down, straight
    let up_result = results_iter.next().unwrap();
    let down_result = results_iter.next().unwrap();
    let straight_result = results_iter.next().unwrap();
    assert!(results_iter.next().is_none());

    // verify up ray
    let (_, data) = up_result.get();
    assert!(increase(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KY_INDEX));
    assert!(increase(data, KX_INDEX));

    // verify the down ray
    let (_, data) = down_result.get();
    assert!(increase(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KY_INDEX));
    assert!(increase(data, KX_INDEX));

    // verify the straight ray
    let (_, data) = straight_result.get();
    assert!(increase(data, XINDEX));
    assert!(same(data, YINDEX));
    assert!(same(data, KY_INDEX));
    assert!(increase(data, KX_INDEX));
}

#[test]
/// test a linear beach on the left side of domain starting in shallow water
///
/// ## Bathymetry
/// `ConstantSlope` object with initial conditions:
/// - `x0 = 0 m`
/// - `y0 = 0 m`
/// - `h0 = 100 m`
/// - `dhdx =  0.05`
/// - `dhdy = 0`
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 2000 m`
/// - `y = 0 m`
/// - `kx = -k * cos(PI/6)`
/// - `ky = k * sin(PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 2000 m`
/// - `y = 0 m`
/// - `kx = -k * cos(-PI/6)`
/// - `ky = k * sin(-PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 2000 m`
/// - `y = 100 m`
/// - `kx = k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from right to left starting in shallow water. The
/// rays start on a beach and immediately start to curve towards the beach,
/// so the kx values of each ray will decrease.
///
/// ## Expected behavior
/// The rays will curve towards the beach
fn test_linear_beach_left() {
    let bathymetry_data = ConstantSlope::builder()
        .x0(0.0)
        .y0(0.0)
        .h0(0.0)
        .dhdx(0.05)
        .dhdy(0.0)
        .build()
        .unwrap();

    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let up_ray = RayState::new(
        Point::new(2_000.0, 0.0),
        WaveNumber::new(-k * (PI / 6.0).cos(), k * (PI / 6.0).sin()),
    );

    let down_ray = RayState::new(
        Point::new(2_000.0, 0.0),
        WaveNumber::new(-k * (-PI / 6.0).cos(), k * (-PI / 6.0).sin()),
    );

    let straight_ray = RayState::new(Point::new(2_000.0, 100.0), WaveNumber::new(-k, 0.0)); // The y value is 100 to avoid a floating point error

    let initial_rays = vec![up_ray, down_ray, straight_ray];

    let waves = ManyRays::new(&bathymetry_data, &current_data, &initial_rays);

    let results = waves.trace_many(0.0, 1_000.0, 1.0);

    let mut results_iter = results.iter().flatten();

    // order is up, down, straight
    let up_result = results_iter.next().unwrap();
    let down_result = results_iter.next().unwrap();
    let straight_result = results_iter.next().unwrap();
    assert!(results_iter.next().is_none());

    // verify up ray
    let (_, data) = up_result.get();
    assert!(decrease(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KY_INDEX));
    assert!(decrease(data, KX_INDEX));

    // verify the down ray
    let (_, data) = down_result.get();
    assert!(decrease(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KY_INDEX));
    assert!(decrease(data, KX_INDEX));

    // verify the straight ray
    let (_, data) = straight_result.get();
    assert!(decrease(data, XINDEX));
    assert!(same(data, YINDEX));
    assert!(same(data, KY_INDEX));
    assert!(decrease(data, KX_INDEX));
}

#[test]
/// test a linear beach on the top side of domain starting in shallow water
///
/// ## Bathymetry
/// `ConstantSlope` object with initial conditions:
/// - `x0 = 0 m`
/// - `y0 = 0 m`
/// - `h0 = 100 m`
/// - `dhdx =  0`
/// - `dhdy = -0.05`
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted left)
/// - `x = 0 m`
/// - `y = 0 m`
/// - `kx = k * cos(4PI/6)`
/// - `ky = k * sin(4PI/6)`
///
/// ### ray 2 (slanted right)
/// - `x = 0 m`
/// - `y = 0 m`
/// - `kx = k * cos(2PI/6)`
/// - `ky = k * sin(2PI/6)`
///
/// ### ray 3 (vertical)
/// - `x = 100 m`
/// - `y = 0 m`
/// - `kx = 0`
/// - `ky = k`
///
/// ## Description
/// The 3 rays propagate from bottom to top starting in shallow water. The
/// rays start on a beach and immediately start to curve towards the beach,
/// so the ky values of each ray will increase.
///
/// ## Expected behavior
/// The rays will curve towards the beach
fn test_linear_beach_top() {
    let bathymetry_data = ConstantSlope::builder()
        .x0(0.0)
        .y0(0.0)
        .h0(100.0)
        .dhdx(0.0)
        .dhdy(-0.05)
        .build()
        .unwrap();

    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let left_ray = RayState::new(
        Point::new(0.0, 0.0),
        WaveNumber::new(k * (4.0 * PI / 6.0).cos(), k * (4.0 * PI / 6.0).sin()),
    );

    let right_ray = RayState::new(
        Point::new(0.0, 0.0),
        WaveNumber::new(k * (2.0 * PI / 6.0).cos(), k * (2.0 * PI / 6.0).sin()),
    );

    let vertical_ray = RayState::new(Point::new(100.0, 0.0), WaveNumber::new(0.0, k)); // the x value is 100 to avoid a floating point error

    let initial_rays = vec![left_ray, right_ray, vertical_ray];

    let waves = ManyRays::new(&bathymetry_data, &current_data, &initial_rays);

    let results = waves.trace_many(0.0, 1_000.0, 1.0);

    let mut results_iter = results.iter().flatten();

    // order is left, right, vertical
    let left_result = results_iter.next().unwrap();
    let right_result = results_iter.next().unwrap();
    let vertical_result = results_iter.next().unwrap();
    assert!(results_iter.next().is_none());

    // verify left ray
    let (_, data) = left_result.get();
    assert!(decrease(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(increase(data, KY_INDEX));

    // verify the down ray
    let (_, data) = right_result.get();
    assert!(increase(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(increase(data, KY_INDEX));

    // verify the straight ray
    let (_, data) = vertical_result.get();
    assert!(same(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(increase(data, KY_INDEX));
}

#[test]
/// test a linear beach on the top side of domain starting in shallow water
///
/// ## Bathymetry
/// `ConstantSlope` object with initial conditions:
/// - `x0 = 0 m`
/// - `y0 = 0 m`
/// - `h0 = 0 m`
/// - `dhdx =  0`
/// - `dhdy = 0.05`
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted left)
/// - `x = 0 m`
/// - `y = 2000 m`
/// - `kx = k * cos(4PI/6)`
/// - `ky = -k * sin(4PI/6)`
///
/// ### ray 2 (slanted right)
/// - `x = 0 m`
/// - `y = 2000 m`
/// - `kx = k * cos(2PI/6)`
/// - `ky = -k * sin(2PI/6)`
///
/// ### ray 3 (vertical)
/// - `x = 100 m`
/// - `y = 2000 m`
/// - `kx = 0`
/// - `ky = -k`
///
/// ## Description
/// The 3 rays propagate from top to bottom starting in shallow water. The
/// rays start on a beach and immediately start to curve towards the beach,
/// so the ky values of each ray will decrease.
///
/// ## Expected behavior
/// The rays will curve towards the beach
fn test_linear_beach_bottom() {
    let bathymetry_data = ConstantSlope::builder()
        .x0(0.0)
        .y0(0.0)
        .h0(0.0)
        .dhdx(0.0)
        .dhdy(0.05)
        .build()
        .unwrap();

    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let left_ray = RayState::new(
        Point::new(0.0, 2_000.0),
        WaveNumber::new(k * (4.0 * PI / 6.0).cos(), -k * (4.0 * PI / 6.0).sin()),
    );

    let right_ray = RayState::new(
        Point::new(0.0, 2_000.0),
        WaveNumber::new(k * (2.0 * PI / 6.0).cos(), -k * (2.0 * PI / 6.0).sin()),
    );

    let vertical_ray = RayState::new(Point::new(100.0, 2_000.0), WaveNumber::new(0.0, -k)); // the x value is 100 to prevent floating point error

    let initial_rays = vec![left_ray, right_ray, vertical_ray];

    let waves = ManyRays::new(&bathymetry_data, &current_data, &initial_rays);

    let results = waves.trace_many(0.0, 1_000.0, 1.0);

    let mut results_iter = results.iter().flatten();

    // order is left, right, vertical
    let left_result = results_iter.next().unwrap();
    let right_result = results_iter.next().unwrap();
    let vertical_result = results_iter.next().unwrap();
    assert!(results_iter.next().is_none());

    // verify left ray
    let (_, data) = left_result.get();
    assert!(decrease(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(decrease(data, KY_INDEX));

    // verify the down ray
    let (_, data) = right_result.get();
    assert!(increase(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(decrease(data, KY_INDEX));

    // verify the straight ray
    let (_, data) = vertical_result.get();
    assert!(same(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(decrease(data, KY_INDEX));
}
