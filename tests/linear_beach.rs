use std::f64::consts::PI;

use mantaray::{
    self, bathymetry::CartesianNetcdf3, current::ConstantCurrent,
    io::utility::create_netcdf3_bathymetry, ray::ManyRays,
};
use tempfile::NamedTempFile;

mod helper;
use helper::*;

#[test]
/// test a linear beach on the right side of domain starting in deep water
///
/// ## Bathymetry file
/// `Lx = 100 km`
///
/// `Ly = 50 km`
///  
/// `dx = dy = 500 m`
///
/// `h0 =  2000 m`
///
/// `h = h0 if x < 50km else h0 - 0.05 * (x - 50km)`
///
/// *minimum depth in the file is 0 m*
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 10km`
/// - `y = 1km`
/// - `kx = k * cos(PI/6)`
/// - `ky = k * sin(PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 10km`
/// - `y = 49 km`
/// - `kx = k * cos(-PI/6)`
/// - `ky = k * sin(-PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 10km`
/// - `y = 25km`
/// - `kx = k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from left to right first in deep water. Then the rays
/// will reach a beach and start to curve towards the beach, so the kx values of
/// each ray will increase.
///
/// ## Expected behavior
/// The rays will go straight in the deep water, but they will begin the curve
/// towards the beach when the interaction with the bathymetry is larger than
/// the difference between fp values.
fn test_deep_linear_beach_right() {
    // create file for the bathymetry data
    let temp_path = NamedTempFile::new().unwrap().into_temp_path();

    // beach on right
    create_netcdf3_bathymetry(&temp_path, 200, 100, 500.0, 500.0, |x, _y| {
        if x < 50_000.0 {
            2000.0
        } else {
            if 2000.0 - 0.05 * (x - 50_000.0) as f64 >= 0.0 {
                2000.0 - 0.05 * (x - 50_000.0) as f64
            } else {
                0.0
            }
        }
    });

    let bathymetry_data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth").unwrap();
    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let up_ray = (
        10_000.0,
        1_000.0,
        k * (PI / 6.0).cos(),
        k * (PI / 6.0).sin(),
    );

    let down_ray = (
        10_000.0,
        49_000.0,
        k * (-PI / 6.0).cos(),
        k * (-PI / 6.0).sin(),
    );

    let straight_ray = (10_000.0, 25_000.0, k, 0.0);

    let initial_rays = vec![up_ray, down_ray, straight_ray];

    let waves = ManyRays::new(&bathymetry_data, Some(&current_data), &initial_rays);

    let results = waves.trace_many(0.0, 100_000.0, 1.0);

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

    // kx same until beach, where it will be >= to start
    assert!(can_increase_after(data, KX_INDEX, |state| state[0] >= 50_000.0));
    // verify last kx should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // verify the down ray
    let (_, data) = down_result.get();
    assert!(increase(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KY_INDEX));

    // kx same until beach, where it will be >= to start
    assert!(can_increase_after(data, KX_INDEX, |state| state[0] >= 50_000.0));
    // verify last kx should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // verify the straight ray
    let (_, data) = straight_result.get();
    assert!(increase(data, XINDEX));
    assert!(same(data, YINDEX));
    assert!(same(data, KY_INDEX));

    // kx same until beach, where it will be >= to start
    assert!(can_increase_after(data, KX_INDEX, |state| state[0] >= 50_000.0));
    // verify last kx should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );
}

#[test]
/// test a linear beach on the left side of domain starting in deep water
///
/// ## Bathymetry file
/// `Lx = 100 km`
///
/// `Ly = 50 km`
///  
/// `dx = dy = 500 m`
///
/// `h0 =  2000 m`
///
/// `h = h0 if x > 50km else 0.05 * (x - 10km)`
///
/// *minimum depth in the file is 0 m*
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 90km`
/// - `y = 1km`
/// - `kx = k * cos(5PI/6)`
/// - `ky = k * sin(5PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 90km`
/// - `y = 49 km`
/// - `kx = k * cos(-5PI/6)`
/// - `ky = k * sin(-5PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 90km`
/// - `y = 25km`
/// - `kx = -k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from right to left first in deep water. Then the rays
/// will reach a beach and start to curve towards the beach, so the kx values of
/// each ray will increase.
///
/// ## Expected behavior
/// The rays will go straight in the deep water, but they will begin the curve
/// towards the beach when the interaction with the bathymetry is larger than
/// the difference between fp values.
fn test_deep_linear_beach_left() {
    // create file for the bathymetry data
    let temp_path = NamedTempFile::new().unwrap().into_temp_path();

    // beach on right
    create_netcdf3_bathymetry(&temp_path, 200, 100, 500.0, 500.0, |x, _y| {
        if x > 50_000.0 {
            2000.0
        } else if x < 10_000.0 {
            0.0
        } else {
            0.05 * (x - 10_000.0) as f64
        }
    });

    let bathymetry_data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth").unwrap();
    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let up_ray = (
        90_000.0,
        1_000.0,
        k * (5.0 * PI / 6.0).cos(),
        k * (5.0 * PI / 6.0).sin(),
    );

    let down_ray = (
        90_000.0,
        49_000.0,
        k * (-5.0 * PI / 6.0).cos(),
        k * (-5.0 * PI / 6.0).sin(),
    );

    let straight_ray = (90_000.0, 25_000.0, -k, 0.0);

    let initial_rays = vec![up_ray, down_ray, straight_ray];

    let waves = ManyRays::new(&bathymetry_data, Some(&current_data), &initial_rays);

    let results = waves.trace_many(0.0, 100_000.0, 1.0);

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

    // kx same until beach, where it will be <= to start
    assert!(can_decrease_after(data, KX_INDEX, |state| state[0] <= 50_000.0));
    // verify last kx should be less than first (started negative)
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // verify the down ray
    let (_, data) = down_result.get();
    assert!(decrease(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KY_INDEX));

    // kx same until beach, where it will be <= to start
    assert!(can_decrease_after(data, KX_INDEX, |state| state[0] <= 50_000.0));
    // verify last kx should be less than first (started negative)
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // verify the straight ray
    let (_, data) = straight_result.get();
    assert!(decrease(data, XINDEX));
    assert!(same(data, YINDEX));
    assert!(same(data, KY_INDEX));

    // kx same until beach, where it will be <= to start
    assert!(can_decrease_after(data, KX_INDEX, |state| state[0] <= 50_000.0));
    // verify last kx should be less than first (started negative)
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );
}
