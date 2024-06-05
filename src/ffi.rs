#![allow(unused)]

// extern crate std;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::Path;
use std::str;

use ode_solvers::dop_shared::SolverResult;
use pyo3::prelude::*;

use crate::bathymetry::CartesianFile;
use crate::ray::SingleRay;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn single_ray(
    x0: f64,
    y0: f64,
    kx0: f64,
    ky0: f64,
    duration: f64,
    step_size: f64,
    filename: String,
) -> PyResult<(Vec<(f64, f64, f64, f64, f64)>)> {
    let bathymetry = CartesianFile::new(Path::new(&filename));
    let wave = SingleRay::new(&bathymetry, x0, y0, kx0, ky0);
    let res = wave.trace_individual(0.0, duration, step_size).unwrap();
    let (t, s) = res.get();
    let ans: Vec<_> = t
        .iter()
        .zip(s.iter())
        .map(|(t, s)| (*t, s[0], s[1], s[2], s[3]))
        .collect();
    Ok(ans)
}

/// A Python module implemented in Rust.
#[pymodule]
fn ray_tracing(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(single_ray, m)?)?;
    Ok(())
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
