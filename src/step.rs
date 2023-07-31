//! Calculate next latitude and longitude position
//! 
//! # Note
//! This requires the use of the crate proj, which was for me tricky to set up.
//! What worked for me was setting up a conda environment with rust, proj,
//! pkgconfig. Then add the location of the environment to the path. export the
//! location of proj.pc to PKG_CONFIG_PATHS. Then I created the build.rs file,
//! which solved the last two linker errors once the code was able to compile.
//! 

use proj::Proj;
use std::f32::consts::PI;

use crate::error::Error;


/// Calculate next latitude and longitude given azimuth and distance according
/// to the projection
/// 
/// The input `lat`, `lon` are converted to meters according to the
/// `projection`. Then a new x, y point is calculated using the azimuth and
/// distance. Then, the new latitude and longitude are calculated by taking the
/// inverse of `projection` on the new x and y.
/// 
/// # Arguments:
/// 
/// `lat`: `&f32`
/// - Latitude of the starting point in degrees.
/// 
/// `lon`: `&f32`
/// - Longitude of the starting point in degrees.
/// 
/// `azimuth`: `&f32`
/// - Direction of travel in degrees clockwise to north.
/// 
/// `distance`: `&f32`
/// - distance in meters to trace forward
/// 
/// `projection`: `&Proj`
/// - struct representing the map projection converting between EPSG:4326 and
///   EPSG:3857: latitude and longitude to x and y in meters.
/// 
/// # Returns
/// `Result<(f32, f32), Error>`
/// - `Ok((f32, f32))`: the new latitude and longitude coordinates
/// - `Err(Error)`: there was an error in `proj::Proj::project`
/// 
/// # Errors
/// `Error::ProjectionError`: this error is returned when `proj::Proj::project`
/// returns an error.
/// 
/// # Note
/// This function is only tested converting between EPSG:4326 and EPSG:3857
/// right now, but I think it could also work for different projections to the same units in the
/// future.
fn trace_forward(lat: &f32, lon: &f32, azimuth: &f32, distance: &f32, projection: &Proj) -> Result<(f32, f32), Error> {

    // convert lat, lon to x, y in meters
    let (x, y) = latlon_to_m(lat, lon, projection, false)?;

    // use distance and azimuth to find new point
    let (x_new, y_new) = (x + distance * (azimuth * PI / 180.0).sin(), y + distance * (azimuth * PI / 180.0).cos());

    // convert new point to lat, lon
    let (lon_new, lat_new) = m_to_latlon(&x_new, &y_new, projection, true)?;

    Ok((lat_new, lon_new))

}

/// Convert latitude and longitude coordinates to meters.
/// 
/// # Arguments:
/// 
/// `lat`: `&f32`
/// - latitude in degrees
/// 
/// `lon`: `&f32`
/// - longitude in degrees
/// 
/// `projection`: `&Proj`
/// - struct representing the map projection converting between latitude and longitude to meters.
/// 
/// `inverse`: `bool`
/// - if `true` using the inverse of the given `projection` instead.
/// 
/// # Returns
/// `Result<(f32, f32), Error>`
/// - `Ok((f32, f32))`: the new x and y coordinates in meters
/// - `Err(Error)`: there was an error in `proj::Proj::project`
/// 
/// # Errors
/// `Error::ProjectionError`: this error is returned when `proj::Proj::project`
/// returns an error.
fn latlon_to_m(lat: &f32, lon: &f32, projection: &Proj, inverse: bool) -> Result<(f32, f32), Error> {
    match projection.project((*lon, *lat), inverse) {
        Ok(p) => Ok(p),
        Err(_) => Err(Error::ProjectionError), 
    }
}

/// Convert meters to latitude and longitude coordinates.
/// 
/// # Arguments:
/// 
/// `x`: `&f32`
/// - x in meters
/// 
/// `y`: `&f32`
/// - y in meters
/// 
/// `projection`: `&Proj`
/// - struct representing the map projection converting between latitude and longitude to meters.
/// 
/// `inverse`: `bool`
/// - if `true` using the inverse of the given `projection` instead.
/// 
/// # Returns
/// `Result<(f32, f32), Error>`
/// - `Ok((f32, f32))`: the new latitude and longitude coordinates in degrees
/// - `Err(Error)`: there was an error in `proj::Proj::project`
/// 
/// # Errors
/// `Error::ProjectionError`: this error is returned when `proj::Proj::project`
/// returns an error.
fn m_to_latlon(x: &f32, y: &f32, projection: &Proj, inverse: bool) -> Result<(f32, f32), Error> {
    match projection.project((*x, *y), inverse) {
        Ok(p) => Ok(p),
        Err(_) => Err(Error::ProjectionError),
    }
}

#[test]
/// test output is the same as calculated from https://epsg.io/transform#s_srs=4326&t_srs=3857&x=45.0000000&y=45.0000000
fn test_latlon_to_m() {
    let from = "EPSG:4326"; // lat, lon (WGS 84)
    let to = "EPSG:3857"; // x, y (web mercator)
    let projection = Proj::new_known_crs(&from, &to, None).unwrap();

    let (x, y) = latlon_to_m(&45.0, &45.0, &projection, false).unwrap();

    assert!((x - 5009377.085697311).abs() < f32::EPSILON);
    assert!((y - 5621521.486192066).abs() < f32::EPSILON);

}

#[test]
/// inverse of `test_latlon_to_m`
fn test_m_to_latlon() {
    let from = "EPSG:4326"; // lat, lon (WGS 84)
    let to = "EPSG:3857"; // x, y (web mercator)
    let projection = Proj::new_known_crs(&from, &to, None).unwrap();

    let (x, y) = m_to_latlon(&5009377.085697311, &5621521.486192066, &projection, true).unwrap();

    assert!((x - 45.0).abs() < f32::EPSILON);
    assert!((y - 45.0).abs() < f32::EPSILON);

}

#[test]
/// tests trace_forward function using calculated case from https://epsg.io/transform#s_srs=3857&t_srs=4326&x=5431995.3474380&y=6527829.2732287
fn test_trace_forward () {

    let from = "EPSG:4326"; // lat, lon (WGS 84)
    let to = "EPSG:3857"; // x, y (web mercator)
    let latlon_to_m = Proj::new_known_crs(&from, &to, None).unwrap();

    let result = trace_forward(&45.0, &45.0, &25.0, &1_000_000.0, &latlon_to_m).unwrap();

    // it seems that proj returns 6 decimal places for lat, lon
    assert!((result.0 - 50.468606).abs() < f32::EPSILON, "Expected 50.468606, recieved {}", result.0);
    assert!((result.1 - 48.796444).abs() < f32::EPSILON, "Expected 48.796444, recieved {}", result.0);

}