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
mod error;
mod interpolator;
mod ray;
mod wave_ray_path;

/// A point in 2D cartesian space
///
/// A `Point` is composed by `x` and `y`, expected to be in meters.
struct Point<T> {
    x: T,
    y: T,
}

#[allow(dead_code)]
impl<T> Point<T> {
    /// Create a new `Point` with the given `x` and `y` coordinates.
    ///
    /// # Examples
    /// ```
    /// let p = Point::new(1.0, 2.0);
    /// assert_eq!(*p.x(), 1.0);
    /// assert_eq!(*p.y(), 2.0);
    /// ```
    fn new(x: T, y: T) -> Self {
        Point { x, y }
    }

    /// Get the `x` coordinate of the `Point`.
    ///
    /// # Examples
    /// ```
    /// let p = Point::new(1.0, 2.0);
    /// assert_eq!(*p.x(), 1.0);
    /// ```
    fn x(&self) -> &T {
        &self.x
    }

    /// Get the `y` coordinate of the `Point`.
    ///
    /// # Examples
    /// ```
    /// let p = Point::new(1.0, 2.0);
    /// assert_eq!(*p.y(), 2.0);
    /// ```
    fn y(&self) -> &T {
        &self.y
    }
}

/// A 2D geolocation in a 2D space
///
/// A `Coordinate` is composed by `lat` and `lon`, expected to be in decimal
/// degrees. For instance, the latitude of the North Pole is 90, and a
/// latitude of -10.5 is equivalent to 10 degrees and 30 minutes South.
struct Coordinate<T> {
    lat: T,
    lon: T,
}

#[allow(dead_code)]
impl<T> Coordinate<T> {
    /// Create a new `Coordinate` with the given `lat` and `lon` coordinates.
    ///
    /// # Examples
    /// ```
    /// let c = Coordinate::new(1.0, 2.0);
    /// assert_eq!(*c.lat(), 1.0);
    /// assert_eq!(*c.lon(), 2.0);
    /// ```
    fn new(lon: T, lat: T) -> Self {
        Coordinate { lat, lon }
    }

    fn lat(&self) -> &T {
        &self.lat
    }

    fn lon(&self) -> &T {
        &self.lon
    }
}

/// The current in a 2D cartesian point
///
/// A `Current` is composed by `u` and `v`, expected to be in meters per
/// second.
struct Current<T> {
    u: T,
    v: T,
}

#[allow(dead_code)]
impl<T> Current<T> {
    fn new(u: T, v: T) -> Self {
        Current { u, v }
    }

    fn u(&self) -> &T {
        &self.u
    }

    fn v(&self) -> &T {
        &self.v
    }
}
