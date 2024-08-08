#![allow(unused)]

// extern crate std;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::Path;
use std::str;

use ode_solvers::dop_shared::SolverResult;
use pyo3::prelude::*;

use crate::bathymetry::CartesianNetcdf3;
use crate::current::CartesianCurrent;
use crate::datatype::{Point, Ray, RayState, WaveNumber};
use crate::ray::{ManyRays, SingleRay};

/// A Python module implemented in Rust.
#[pymodule]
fn _mantaray(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(single_ray, m)?)?;
    m.add_function(wrap_pyfunction!(ray_tracing, m)?)?;
    Ok(())
}

#[pyfunction]
fn single_ray(
    x0: f64,
    y0: f64,
    kx0: f64,
    ky0: f64,
    duration: f64,
    step_size: f64,
    bathymetry_filename: String,
    current_filename: String,
) -> PyResult<(Vec<(f64, f64, f64, f64, f64)>)> {
    let bathymetry = CartesianNetcdf3::open(Path::new(&bathymetry_filename), "x", "y", "depth")
        .expect("could not open bathymetry file");
    let current = CartesianCurrent::open(Path::new(&current_filename), "x", "y", "u", "v");
    let initial_state = RayState::new(Point::new(x0, y0), WaveNumber::new(kx0, ky0));
    let wave = SingleRay::new(&bathymetry, &current, &initial_state);
    let res = wave.trace_individual(0.0, duration, step_size).unwrap();
    let (t, s) = res.get();
    let ans: Vec<_> = t
        .iter()
        .zip(s.iter())
        .map(|(t, s)| (*t, s[0], s[1], s[2], s[3]))
        .collect();
    Ok(ans)
}

#[pyfunction]
fn ray_tracing(
    x0: Vec<f64>,
    y0: Vec<f64>,
    kx0: Vec<f64>,
    ky0: Vec<f64>,
    duration: f64,
    step_size: f64,
    bathymetry_filename: String,
    current_filename: String,
) -> PyResult<(Vec<Vec<(f64, f64, f64, f64, f64)>>)> {
    let bathymetry = CartesianNetcdf3::open(Path::new(&bathymetry_filename), "x", "y", "depth")
        .expect("could not open bathymetry file");
    let current = CartesianCurrent::open(Path::new(&current_filename), "x", "y", "u", "v");
    let init_cond = x0
        .iter()
        .zip(y0.iter())
        .zip(kx0.iter().zip(ky0.iter()))
        .map(|((x, y), (kx, ky))| RayState::new(Point::new(*x, *y), WaveNumber::new(*kx, *ky)))
        .collect::<Vec<RayState<f64>>>();
    let waves = ManyRays::new(&bathymetry, &current, &init_cond);
    let res = waves.trace_many(0.0, duration, step_size);
    let rays: Vec<Vec<(f64, f64, f64, f64, f64)>> = res
        .iter()
        .filter_map(|r| r.as_ref())
        .map(|r| {
            let (t, s) = r.get();
            t.iter()
                .zip(s.iter())
                .map(|(t, s)| (*t, s[0], s[1], s[2], s[3]))
                .collect::<Vec<_>>()
        })
        .collect();
    Ok(rays)
}

/*
#[no_mangle]
pub unsafe extern "C" fn single_ray(
    bathymetry_path: *const c_char,
    x0: f64,
    y0: f64,
    kx0: f64,
    ky0: f64,
    end_time: f64,
    step_size: f64,
) -> i32 {
    let bytes = unsafe { CStr::from_ptr(bathymetry_path).to_bytes() };
    let str_slice = str::from_utf8(bytes).unwrap();
    let path = Path::new(str_slice);
    let bathymetry = CartesianFile::new(path);
    let wave = SingleRay::new(&bathymetry, x0, y0, kx0, ky0);
    let res = wave.trace_individual(0.0, end_time, step_size).unwrap();
    dbg!(&res);
    42
}
*/
