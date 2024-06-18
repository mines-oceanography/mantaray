use tempfile::NamedTempFile;

use mantaray::bathymetry::CartesianNetcdf3;
use mantaray::io::utility::create_netcdf3_bathymetry;
use mantaray::ray::ManyRays;

#[cfg(test)]
#[test]
/// Horizontal ray in deep water with constant depth
/// ## Initial conditions:
///
/// kx0 = 0.05 m^-1
///
/// ky0 = 0
///
/// h = 2000 m
///
/// (u, v) = (0, 0)
///
/// ## Description:
///
/// Rays propagate from left to right over constant bathymetry in deep water.
/// Initial direction theta0 = 0.
///
/// ## Expected behavior:
///
/// The ray path goes from left to right, parallel to the x-axis and the values
/// of (kx, ky) are equal to the initial condition at all points.
fn horizontal_constant_depth_deep() {
    /// h = 2000 m
    fn depth_fn(_x: f32, _y: f32) -> f64 {
        2000.0
    }

    // temporary file
    let tmp_file = NamedTempFile::new().unwrap();
    let tmp_path = tmp_file.into_temp_path();

    // x length = 100km, y length = 50 km, and dx = dy = 0.5 km
    create_netcdf3_bathymetry(&tmp_path, 200, 100, 500.0, 500.0, depth_fn);

    // load the data
    let bathymetry_data = CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();
    let current_data = None; // default (u, v) = (0, 0)

    // rays propagate from left side, kx = 0.05 m^-1, and ky = 0.0
    let init_rays: Vec<(f64, f64, f64, f64)> = (1..=9)
        .map(|i| (5000.0, (i * 5000) as f64, 0.05, 0.0))
        .collect();

    let ystart: Vec<f64> = init_rays.clone().iter().map(|(_, y, _, _)| *y).collect();

    let rays = ManyRays::new(&bathymetry_data, current_data, &init_rays);

    let results = rays.trace_many(0.0, 10000.0, 1.0);

    for (i, ray) in results.iter().flatten().enumerate() {
        let (_, data) = &ray.get();

        // x value increases
        let mut last_x = data[0][0];
        for r in data.iter().skip(1) {
            assert!(r[0] > last_x);
            last_x = r[0];
        }

        // y, kx, ky values stay the same
        data.iter().for_each(|r| assert_eq!(r[1], ystart[i])); // y
        data.iter().for_each(|r| assert_eq!(r[2], 0.05)); // kx
        data.iter().for_each(|r| assert_eq!(r[3], 0.0)); // ky
    }
}

#[cfg(test)]
#[test]
/// Horizontal ray in shallow water with constant depth
///
/// ## Initial conditions:
///
/// kx0 = 0.05 m^-1
///
///  ky0 = 0
///
/// h = 10 m
///
/// (u, v) = (0, 0)
///
/// ## Description:
/// Rays propagate from left to right over constant bathymetry in shallow water.
/// Initial direction theta0 = 0.
///
/// ## Expected behavior:
/// The ray path goes from left to right, parallel to the x-axis and the values
/// of (kx, ky) are equal to the initial condition at all points.
fn horizontal_constant_depth_shallow() {
    /// h = 10 m
    fn depth_fn(_x: f32, _y: f32) -> f64 {
        10.0
    }

    // temporary file
    let tmp_file = NamedTempFile::new().unwrap();
    let tmp_path = tmp_file.into_temp_path();

    // x length = 100km, y length = 50 km, and dx = dy = 0.5 km
    create_netcdf3_bathymetry(&tmp_path, 200, 100, 500.0, 500.0, depth_fn);

    // load the data
    let bathymetry_data = CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();
    let current_data = None;

    // rays propagate from left side, kx = 0.05 m^-1, and ky = 0.0
    let init_rays: Vec<(f64, f64, f64, f64)> = (1..=9)
        .map(|i| (5000.0, (i * 5000) as f64, 0.05, 0.0))
        .collect();

    let ystart: Vec<f64> = init_rays.clone().iter().map(|(_, y, _, _)| *y).collect();

    let rays = ManyRays::new(&bathymetry_data, current_data, &init_rays);

    let results = rays.trace_many(0.0, 1000.0, 1.0);

    for (i, ray) in results.iter().flatten().enumerate() {
        let (_, data) = &ray.get();

        // x value increases
        let mut last_x = data[0][0];
        for r in data.iter().skip(1) {
            assert!(r[0] > last_x);
            last_x = r[0];
        }

        // y, kx, ky values stay the same
        data.iter().for_each(|r| assert_eq!(r[1], ystart[i])); // y
        data.iter().for_each(|r| assert_eq!(r[2], 0.05)); // kx
        data.iter().for_each(|r| assert_eq!(r[3], 0.0)); // ky
    }
}

#[cfg(test)]
#[test]
/// Slanted ray in deep water with constant depth
///
/// ## Initial conditions:
///
/// kx0 = 0.05 m^-1
///
/// ky0 = 0.04 m^-1
///
/// h = 2000 m
///
/// (u, v) = (0, 0)
///
/// ## Description:
/// Rays propagate from left to right over constant bathymetry in deep water.
/// Initial direction theta0 approximately 38.7 degrees
///
/// ## Expected behavior:
/// The ray path goes from left to right at an angle theta approximately 38.7
/// degrees. The values of (kx, ky) are equal to the initial condition at all
/// points
fn slanted_constant_depth_deep() {
    // h = 2000 m (using more concise closure notation)
    let depth_fn = |_, _| 2000.0;

    // temporary file
    let tmp_file = NamedTempFile::new().unwrap();
    let tmp_path = tmp_file.into_temp_path();

    // x length = 100km, y length = 50 km, and dx = dy = 0.5 km
    create_netcdf3_bathymetry(&tmp_path, 200, 100, 500.0, 500.0, depth_fn);

    // load the data
    let bathymetry_data = CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();
    let current_data = None;

    // rays propagate from left side, kx = 0.05 m^-1, and ky = 0.0
    let init_rays: Vec<(f64, f64, f64, f64)> = (1..=9)
        .map(|i| (5000.0, (i * 5000) as f64, 0.05, 0.04))
        .collect();

    let rays = ManyRays::new(&bathymetry_data, current_data, &init_rays);

    let results = rays.trace_many(0.0, 10000.0, 1.0);

    for ray in results.iter().flatten() {
        let (_, data) = &ray.get();

        // x and y value increases
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter().skip(1) {
            if r[0].is_nan() || r[1].is_nan() {
                break; // reached end of domain
            }
            assert!(r[0] > last_x);
            assert!(r[1] > last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // kx and ky values stay the same (filtering out NaNs)
        data.iter()
            .filter(|v| !v[0].is_nan())
            .for_each(|r| assert_eq!(r[2], 0.05)); // kx
        data.iter()
            .filter(|v| !v[0].is_nan())
            .for_each(|r| assert_eq!(r[3], 0.04)); // ky
    }
}

#[cfg(test)]
#[test]
/// Slanted ray in shallow water with constant depth
///
/// ## Initial conditions:
///
/// kx0 = 0.05 m^-1
///
/// ky0 = 0.04 m^-1
///
/// h = 10 m
///
/// (u, v) = (0, 0)
///
/// ## Description:
/// Rays propagate from left to right over constant bathymetry in shallow water.
/// Initial direction theta0 approximately 38.7 degrees.
///
/// ## Expected behavior:
/// The ray path goes from left to right at an angle theta approximately
/// 38.7 degrees. The values of (kx, ky) are equal to the initial condition
/// at all points
fn slanted_constant_depth_shallow() {
    // h = 10 m (using more concise closure notation)
    let depth_fn = |_, _| 10.0;

    // temporary file
    let tmp_file = NamedTempFile::new().unwrap();
    let tmp_path = tmp_file.into_temp_path();

    // x length = 100km, y length = 50 km, and dx = dy = 0.5 km
    create_netcdf3_bathymetry(&tmp_path, 200, 100, 500.0, 500.0, depth_fn);

    // load the data
    let bathymetry_data = CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();
    let current_data = None;

    // rays propagate from left side, kx = 0.05 m^-1, and ky = 0.0
    let init_rays: Vec<(f64, f64, f64, f64)> = (1..=9)
        .map(|i| (5000.0, (i * 5000) as f64, 0.05, 0.04))
        .collect();

    let rays = ManyRays::new(&bathymetry_data, current_data, &init_rays);

    let results = rays.trace_many(0.0, 1000.0, 1.0);

    for ray in results.iter().flatten() {
        let (_, data) = &ray.get();

        // x and y value increases
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter().skip(1) {
            if r[0].is_nan() || r[1].is_nan() {
                break; // reached end of domain
            }
            assert!(r[0] > last_x);
            assert!(r[1] > last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // kx and ky values stay the same (filtering out NaNs)
        data.iter()
            .filter(|v| !v[0].is_nan())
            .for_each(|r| assert_eq!(r[2], 0.05)); // kx
        data.iter()
            .filter(|v| !v[0].is_nan())
            .for_each(|r| assert_eq!(r[3], 0.04)); // ky
    }
}
