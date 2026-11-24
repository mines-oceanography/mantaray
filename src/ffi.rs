#![allow(unused)]

// extern crate std;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::Path;
use std::str;

use ode_solvers::dop_shared::SolverResult;
use numpy::{PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::bathymetry::CartesianNetcdf3;
use crate::current::{CartesianCurrent, ConstantCurrent, CurrentData};
use crate::datatype::{Point, Ray, RayState, WaveNumber};
use crate::ray::{ManyRays, SingleRay};

/// A Python module implemented in Rust.
#[pymodule]
fn _mantaray(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(single_ray, m)?)?;
    m.add_function(wrap_pyfunction!(ray_tracing, m)?)?;
    m.add_function(wrap_pyfunction!(ray_tracing_arrays, m)?)?;
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

/// Ray tracing where bathymetry (and optionally current) are given as arrays,
/// with no files involved. `depth`, `u` and `v` have shape `(len(y), len(x))`.
/// Rays that fail during integration are returned as a single NaN row so that
/// the ray index is preserved.
#[pyfunction]
#[pyo3(signature = (x0, y0, kx0, ky0, duration, step_size, x, y, depth, u=None, v=None))]
#[allow(clippy::too_many_arguments)]
fn ray_tracing_arrays(
    py: Python<'_>,
    x0: Vec<f64>,
    y0: Vec<f64>,
    kx0: Vec<f64>,
    ky0: Vec<f64>,
    duration: f64,
    step_size: f64,
    x: PyReadonlyArray1<'_, f64>,
    y: PyReadonlyArray1<'_, f64>,
    depth: PyReadonlyArray2<'_, f64>,
    u: Option<PyReadonlyArray2<'_, f64>>,
    v: Option<PyReadonlyArray2<'_, f64>>,
) -> PyResult<Vec<Vec<(f64, f64, f64, f64, f64)>>> {
    let err = |e: crate::error::Error| PyValueError::new_err(e.to_string());
    if [y0.len(), kx0.len(), ky0.len()].iter().any(|n| *n != x0.len()) {
        return Err(PyValueError::new_err("x0, y0, kx0 and ky0 must have the same length"));
    }
    let (nx, ny) = (x.len(), y.len());
    if depth.shape() != [ny, nx] {
        return Err(PyValueError::new_err("depth must have shape (len(y), len(x))"));
    }
    let xs = x.as_array().to_vec();
    let ys = y.as_array().to_vec();
    let bathymetry =
        CartesianNetcdf3::from_arrays(xs.clone(), ys.clone(), depth.as_array().iter().copied().collect())
            .map_err(err)?;

    let current: Box<dyn CurrentData> = match (u, v) {
        (Some(u), Some(v)) => {
            if u.shape() != [ny, nx] || v.shape() != [ny, nx] {
                return Err(PyValueError::new_err("u and v must have shape (len(y), len(x))"));
            }
            Box::new(
                CartesianCurrent::from_arrays(
                    xs,
                    ys,
                    u.as_array().iter().copied().collect(),
                    v.as_array().iter().copied().collect(),
                )
                .map_err(err)?,
            )
        }
        (None, None) => Box::new(ConstantCurrent::new(0.0, 0.0)),
        _ => return Err(PyValueError::new_err("u and v must be given together")),
    };

    let init_cond = x0
        .iter()
        .zip(y0.iter())
        .zip(kx0.iter().zip(ky0.iter()))
        .map(|((x, y), (kx, ky))| RayState::new(Point::new(*x, *y), WaveNumber::new(*kx, *ky)))
        .collect::<Vec<RayState<f64>>>();

    let res = py.detach(|| {
        ManyRays::new(&bathymetry, current.as_ref(), &init_cond).trace_many(0.0, duration, step_size)
    });
    Ok(res
        .iter()
        .map(|r| match r {
            Some(r) => {
                let (t, s) = r.get();
                t.iter()
                    .zip(s.iter())
                    .map(|(t, s)| (*t, s[0], s[1], s[2], s[3]))
                    .collect()
            }
            None => vec![(f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN)],
        })
        .collect())
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
