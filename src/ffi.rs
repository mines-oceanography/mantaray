#![allow(unused)]

extern crate std;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::Path;
use std::str;

use crate::bathymetry::BathymetryData;
use crate::bathymetry::CartesianFile;
use crate::ray::SingleRay;

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
    for r in res.1.iter() {
        println!("{}", r);
    }
    42
}
