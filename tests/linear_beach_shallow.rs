// FIXME

use std::{f64::consts::PI, path::Path};

use mantaray::{
    self,
    bathymetry::CartesianNetcdf3,
    current::ConstantCurrent,
    io::utility::create_netcdf3_bathymetry,
    ray::{output_or_append_to_tsv_file, ManyRays},
};
use tempfile::NamedTempFile;

mod helper;
use helper::*;

#[test]
/// test a linear beach on the right side of domain starting in deep water
///
/// ## Bathymetry file
/// `Lx = 4 km`
///
/// `Ly = 2 km`
///  
/// `dx = dy = 40 m`
///
/// `h0 =  100 m`
///
/// `h = h0 if x < 1km else h0 - 0.05 * (x - 1km)`
///
/// *minimum depth in the file is 0 m*
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 100 m`
/// - `y = 100 m`
/// - `kx = k * cos(PI/6)`
/// - `ky = k * sin(PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 100 m`
/// - `y = 1900 m`
/// - `kx = k * cos(-PI/6)`
/// - `ky = k * sin(-PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 100 m`
/// - `y = 1000 m`
/// - `kx = k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from left to right starting in shallow water. Then the rays
/// will reach a beach and immediately start to curve towards the beach, so the kx values of
/// each ray will increase.
///
/// ## Expected behavior
/// The rays will go straight in the deep water, but they will begin the curve
/// towards the beach when the interaction with the bathymetry is larger than
/// the difference between fp values. Since the rays start in shallow water, the
/// change in wave number due to the beach will instantly be noticeable.
fn test_shallow_linear_beach_right() {
    // create file for the bathymetry data
    let temp_path = NamedTempFile::new().unwrap().into_temp_path();

    // beach on right
    create_netcdf3_bathymetry(&temp_path, 80, 40, 50.0, 50.0, |x, _y| {
        if x < 1_000.0 {
            100.0
        } else if x < 3_000.0 {
            100.0 - 0.05 * (x - 1_000.0) as f64
        } else {
            0.0
        }
    });

    let bathymetry_data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth").unwrap();
    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let up_ray = (100.0, 100.0, k * (PI / 6.0).cos(), k * (PI / 6.0).sin());

    let down_ray = (100.0, 1_900.0, k * (-PI / 6.0).cos(), k * (-PI / 6.0).sin());

    let straight_ray = (100.0, 1_000.0, k, 0.0);

    let initial_rays = vec![up_ray, down_ray, straight_ray];

    let waves = ManyRays::new(&bathymetry_data, Some(&current_data), &initial_rays);

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

    // kx same until beach, where it will be >= to start
    assert!(can_increase_after(data, KX_INDEX, |state| state[0] > 975.0));
    // verify last kx should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // determine the x value where the kx is increasing
    let mut last_kx = data[0][2];
    for state in data.iter().skip(1) {
        let x = state[0];
        let kx = state[2];
        if kx != last_kx {
            println!("The kx value changes at x = {}", x);
            break;
        }
        last_kx = kx;
    }

    // verify the down ray
    let (_, data) = down_result.get();
    assert!(increase(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KY_INDEX));

    // kx same until beach, where it will be >= to start
    assert!(can_increase_after(data, KX_INDEX, |state| state[0] > 975.0));
    // verify last kx should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // determine the x value where the kx is increasing
    let mut last_kx = data[0][2];
    for state in data.iter().skip(1) {
        let x = state[0];
        let kx = state[2];
        if kx != last_kx {
            println!("The kx value changes at x = {}", x);
            break;
        }
        last_kx = kx;
    }

    // verify the straight ray
    let (_, data) = straight_result.get();
    assert!(increase(data, XINDEX));
    assert!(same(data, YINDEX));
    assert!(same(data, KY_INDEX));

    // kx same until beach, where it will be >= to start
    assert!(can_increase_after(data, KX_INDEX, |state| state[0] > 975.0));
    // verify last kx should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // determine the x value where the kx is increasing
    let mut last_kx = data[0][2];
    for state in data.iter().skip(1) {
        let x = state[0];
        let kx = state[2];
        if kx != last_kx {
            println!("The kx value changes at x = {}", x);
            break;
        }
        last_kx = kx;
    }
}

#[test]
/// test a linear beach on the right side of domain starting in deep water
///
/// ## Bathymetry file
/// `Lx = 4 km`
///
/// `Ly = 2 km`
///  
/// `dx = dy = 40 m`
///
/// `h0 =  100 m`
///
/// `h = h0 if x > 3km else 0.05 * (x - 1km)`
///
/// *minimum depth in the file is 0 m*
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 3500 m`
/// - `y = 100 m`
/// - `kx = k * cos(5PI/6)`
/// - `ky = k * sin(5PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 3500 m`
/// - `y = 1900 m`
/// - `kx = k * cos(-5PI/6)`
/// - `ky = k * sin(-5PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 3500 m`
/// - `y = 1000 m`
/// - `kx = k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from right to left starting in shallow water. Then the rays
/// will reach a beach and immediately start to curve towards the beach, so the kx values of
/// each ray will decrease, since they start negative.
///
/// ## Expected behavior
/// The rays will go straight in the deep water, but they will begin the curve
/// towards the beach when the interaction with the bathymetry is larger than
/// the difference between fp values. Since the rays start in shallow water, the
/// change in wave number due to the beach will instantly be noticeable.
fn test_shallow_linear_beach_left() {
    // create file for the bathymetry data
    let temp_path = NamedTempFile::new().unwrap().into_temp_path();

    // beach on right
    create_netcdf3_bathymetry(&temp_path, 80, 40, 50.0, 50.0, |x, _y| {
        if x < 1_000.0 {
            0.0
        } else if x < 3_000.0 {
            0.05 * (x - 1_000.0) as f64
        } else {
            100.0
        }
    });

    let bathymetry_data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth").unwrap();
    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let up_ray = (
        3_500.0,
        100.0,
        k * (5.0 * PI / 6.0).cos(),
        k * (5.0 * PI / 6.0).sin(),
    );

    let down_ray = (
        3_500.0,
        1_900.0,
        k * (-5.0 * PI / 6.0).cos(),
        k * (-5.0 * PI / 6.0).sin(),
    );

    let straight_ray = (3_500.0, 1_000.0, -k, 0.0);

    let initial_rays = vec![up_ray, down_ray, straight_ray];

    let waves = ManyRays::new(&bathymetry_data, Some(&current_data), &initial_rays);

    let results = waves.trace_many(0.0, 1_000.0, 1.0);

    let mut results_iter = results.iter().flatten();

    // for result in results_iter.clone() {
    //     output_or_append_to_tsv_file(Path::new("shallow_left.txt"), result);
    // }

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

    // kx same until beach, where it will be less than the start (because
    // starting in negative direction)
    assert!(can_decrease_after(data, KX_INDEX, |state| state[0] < 3_025.0));
    // verify last kx should be less than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // determine the x value where the kx is increasing
    let mut last_kx = data[0][2];
    for state in data.iter().skip(1) {
        let x = state[0];
        let kx = state[2];
        if kx != last_kx {
            println!("The kx value changes at x = {}", x);
            break;
        }
        last_kx = kx;
    }

    // verify the down ray
    let (_, data) = down_result.get();
    assert!(decrease(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KY_INDEX));

    assert!(can_decrease_after(data, KX_INDEX, |state| state[0] < 3_025.0));
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // determine the x value where the kx is increasing
    let mut last_kx = data[0][2];
    for state in data.iter().skip(1) {
        let x = state[0];
        let kx = state[2];
        if kx != last_kx {
            println!("The kx value changes at x = {}", x);
            break;
        }
        last_kx = kx;
    }

    // verify the straight ray
    let (_, data) = straight_result.get();
    assert!(decrease(data, XINDEX));
    assert!(same(data, YINDEX));
    assert!(same(data, KY_INDEX));

    assert!(can_decrease_after(data, KX_INDEX, |state| state[0] < 3_025.0));
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KX_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KX_INDEX]
    );

    // determine the x value where the kx is increasing
    let mut last_kx = data[0][2];
    for state in data.iter().skip(1) {
        let x = state[0];
        let kx = state[2];
        if kx != last_kx {
            println!("The kx value changes at x = {}", x);
            break;
        }
        last_kx = kx;
    }
}

#[test]
/// test a linear beach on the right side of domain starting in deep water
///
/// ## Bathymetry file
/// `Lx = 2 km`
///
/// `Ly = 4 km`
///  
/// `dx = dy = 40 m`
///
/// `h0 =  100 m`
///
/// `h = h0 if y < 1km else h0 - 0.05 * (y - 1km)`
///
/// *minimum depth in the file is 0 m*
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted up)
/// - `x = 100 m`
/// - `y = 100 m`
/// - `kx = k * cos(2PI/6)`
/// - `ky = k * sin(2PI/6)`
///
/// ### ray 2 (slanted down)
/// - `x = 1900 m`
/// - `y = 100 m`
/// - `kx = k * cos(4PI/6)`
/// - `ky = k * sin(4PI/6)`
///
/// ### ray 3 (horizontal)
/// - `x = 1000 m`
/// - `y = 100 m`
/// - `kx = 0`
/// - `ky = k`
///
/// ## Description
/// The 3 rays propagate from bottom to top starting in shallow water. Then the rays
/// will reach a beach and immediately start to curve towards the beach, so the ky values of
/// each ray will increase.
///
/// ## Expected behavior
/// The rays will go straight in the deep water, but they will begin the curve
/// towards the beach when the interaction with the bathymetry is larger than
/// the difference between fp values. Since the rays start in shallow water, the
/// change in wave number due to the beach will instantly be noticeable.
fn test_shallow_linear_beach_top() {
    // create file for the bathymetry data
    let temp_path = NamedTempFile::new().unwrap().into_temp_path();

    // beach on right
    create_netcdf3_bathymetry(&temp_path, 40, 80, 50.0, 50.0, |_x, y| {
        if y < 1_000.0 {
            100.0
        } else if y < 3_000.0 {
            100.0 - 0.05 * (y - 1_000.0) as f64
        } else {
            0.0
        }
    });

    let bathymetry_data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth").unwrap();
    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let right_ray = (
        100.0,
        100.0,
        k * (2.0 * PI / 6.0).cos(),
        k * (2.0 * PI / 6.0).sin(),
    );

    let left_ray = (
        1_900.0,
        100.0,
        k * (4.0 * PI / 6.0).cos(),
        k * (4.0 * PI / 6.0).sin(),
    );

    let vertical_ray = (1_000.0, 100.0, 0.0, k);

    let initial_rays = vec![left_ray, right_ray, vertical_ray];

    let waves = ManyRays::new(&bathymetry_data, Some(&current_data), &initial_rays);

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

    // ky same until beach, where it will increase
    assert!(can_increase_after(data, KY_INDEX, |state| state[YINDEX] > 975.0));
    // verify last ky should be greater than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KY_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KY_INDEX]
    );

    // determine the y value where the ky is increasing
    let mut last_ky = data[0][3];
    for state in data.iter().skip(1) {
        let y = state[1];
        let ky = state[3];
        if ky != last_ky {
            println!("The ky value changes at y = {}", y);
            break;
        }
        last_ky = ky;
    }

    // verify the right ray
    let (_, data) = right_result.get();
    assert!(increase(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(can_increase_after(data, KY_INDEX, |state| state[YINDEX] > 975.0));
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KY_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KY_INDEX]
    );

    // determine the y value where the ky is increasing
    let mut last_ky = data[0][3];
    for state in data.iter().skip(1) {
        let y = state[1];
        let ky = state[3];
        if ky != last_ky {
            println!("The ky value changes at y = {}", y);
            break;
        }
        last_ky = ky;
    }

    // verify the vertical ray
    let (_, data) = vertical_result.get();
    assert!(same(data, XINDEX));
    assert!(increase(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(can_increase_after(data, KY_INDEX, |state| state[YINDEX] > 975.0));
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KY_INDEX]
            > data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KY_INDEX]
    );

    // determine the y value where the ky is increasing
    let mut last_ky = data[0][3];
    for state in data.iter().skip(1) {
        let y = state[1];
        let ky = state[3];
        if ky != last_ky {
            println!("The ky value changes at y = {}", y);
            break;
        }
        last_ky = ky;
    }
}

#[test]
/// test a linear beach on the right side of domain starting in deep water
///
/// ## Bathymetry file
/// `Lx = 2 km`
///
/// `Ly = 4 km`
///  
/// `dx = dy = 40 m`
///
/// `h0 =  100 m`
///
/// `h = h0 if y < 1km else h0 - 0.05 * (y - 1km)`
///
/// *minimum depth in the file is 0 m*
///
/// ## Initial conditions
/// 3 rays with different angles and starting locations `k = 0.05;`
///
/// ### ray 1 (slanted right)
/// - `x = 100 m`
/// - `y = 3900 m`
/// - `kx = k * cos(-2PI/6)`
/// - `ky = k * sin(-2PI/6)`
///
/// ### ray 2 (slanted left)
/// - `x = 1900 m`
/// - `y = 3900 m`
/// - `kx = k * cos(-4PI/6)`
/// - `ky = k * sin(-4PI/6)`
///
/// ### ray 3 (vertical)
/// - `x = 1000 m`
/// - `y = 3900 m`
/// - `kx = -k`
/// - `ky = 0`
///
/// ## Description
/// The 3 rays propagate from top to bottom starting in shallow water. Then the rays
/// will reach a beach and immediately start to curve towards the beach, so the ky values of
/// each ray will decrease because ky starts negative.
///
/// ## Expected behavior
/// The rays will go straight in the deep water, but they will begin the curve
/// towards the beach when the interaction with the bathymetry is larger than
/// the difference between fp values. Since the rays start in shallow water, the
/// change in wave number due to the beach will instantly be noticeable.
fn test_shallow_linear_beach_bottom() {
    // create file for the bathymetry data
    let temp_path = NamedTempFile::new().unwrap().into_temp_path();

    let temp_path = Path::new("shallow_bottom.nc");

    // beach on right
    create_netcdf3_bathymetry(&temp_path, 40, 80, 50.0, 50.0, |_x, y| {
        if y < 1_000.0 {
            0.0
        } else if y < 3_000.0 {
            0.05 * (y - 1_000.0) as f64
        } else {
            100.0
        }
    });

    let bathymetry_data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth").unwrap();
    let current_data = ConstantCurrent::new(0.0, 0.0);

    let k = 0.05;

    let right_ray = (
        100.0,
        3_900.0,
        k * (-2.0 * PI / 6.0).cos(),
        k * (-2.0 * PI / 6.0).sin(),
    );

    let left_ray = (
        1_900.0,
        3_900.0,
        k * (-4.0 * PI / 6.0).cos(),
        k * (-4.0 * PI / 6.0).sin(),
    );

    let vertical_ray = (1_000.0, 3_900.0, 0.0, -k);

    let initial_rays = vec![left_ray, right_ray, vertical_ray];

    let waves = ManyRays::new(&bathymetry_data, Some(&current_data), &initial_rays);

    let results = waves.trace_many(0.0, 1_000.0, 1.0);

    let mut results_iter = results.iter().flatten();

    for result in results_iter.clone() {
        output_or_append_to_tsv_file(Path::new("shallow_bottom.txt"), result);
    }

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

    // ky same until beach, where it will decrease
    assert!(can_decrease_after(data, KY_INDEX, |state| state[YINDEX] > 975.0));
    // verify last ky should be less than first
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KY_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KY_INDEX]
    );

    // determine the y value where the ky is decreasing
    let mut last_ky = data[0][3];
    for state in data.iter().skip(1) {
        let y = state[1];
        let ky = state[3];
        if ky != last_ky {
            println!("The ky value changes at y = {}", y);
            break;
        }
        last_ky = ky;
    }

    // verify the right ray
    let (_, data) = right_result.get();
    assert!(increase(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(can_decrease_after(data, KY_INDEX, |state| state[YINDEX] > 975.0));
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KY_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KY_INDEX]
    );

    // determine the y value where the ky is decreasing
    let mut last_ky = data[0][3];
    for state in data.iter().skip(1) {
        let y = state[1];
        let ky = state[3];
        if ky != last_ky {
            println!("The ky value changes at y = {}", y);
            break;
        }
        last_ky = ky;
    }

    // verify the vertical ray
    let (_, data) = vertical_result.get();
    assert!(same(data, XINDEX));
    assert!(decrease(data, YINDEX));
    assert!(same(data, KX_INDEX));
    assert!(can_decrease_after(data, KY_INDEX, |state| state[YINDEX] > 975.0));
    assert!(
        data.iter().filter(|v| !v[0].is_nan()).last().unwrap()[KY_INDEX]
            < data.iter().filter(|v| !v[0].is_nan()).next().unwrap()[KY_INDEX]
    );

    // determine the y value where the ky is increasing
    let mut last_ky = data[0][3];
    for state in data.iter().skip(1) {
        let y = state[1];
        let ky = state[3];
        if ky != last_ky {
            println!("The ky value changes at y = {}", y);
            break;
        }
        last_ky = ky;
    }
}