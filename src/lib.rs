//! Ray tracing ocean waves
//!
//! This library uses ode_solvers, netcdf3, nalgebra, and thiserror.
//!
//! As of 2023-08-10, the library contains a module `ray.rs` that has a struct
//! `SingleRay`. `SingleRay` has a method `trace_individual` to create a `WaveRayPath`, perform `ode_solvers`
//! Runge-Kutta4 algorithm, and return the result. This `lib.rs` module contains a `WaveRayPath`
//! struct that contains either a ConstantDepth, ArrayDepth, or CartesianFile.
//! The struct implements the ode_solvers `system` method, which is defined with
//! the helper function `odes`. The `odes` function uses the `group_velocity`,
//! `gradient`, and `dk_vector_dt` to calculate the derivatives at the current
//! state. The Rk4 is used similar to the
//! [examples](https://srenevey.github.io/ode-solvers/examples/kepler_orbit.html).
//!
//! There is also a file output_to_file, which runs the Rk4, then saves the
//! output to a file. There is a folder named support that contains the python
//! file plot_ode_solvers, which plots a single ray. This will also contain a
//! plotting tool for many rays in the future.
//!
//! This only does one ray at the moment and verified for constant depth waves,
//! but in the future, it will include variable depth, ray bundles, and current.

// enforce documentation
#![deny(missing_docs)]

mod bathymetry; 
mod current;
mod datatype;
mod error;
mod ffi;
mod interpolator;
mod io;
mod ray;
mod ray_result;
#[cfg(test)]
mod tests;
/// cbindgen:ignore
mod wave_ray_path;

#[allow(unused_imports)]
use datatype::{Coordinate, Current, Point};
#[allow(unused_imports)]
pub(crate) use wave_ray_path::State;
