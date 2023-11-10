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

use bathymetry::BathymetryData;
use current::CurrentData;
use derive_builder::Builder;
use ode_solvers::*;
use std::fs::File;
use std::io::Write;

mod error;

mod interpolator;

mod bathymetry;
mod current;

mod ray;

use error::Error;

// Define constant gravity
const G: f64 = 9.8;

impl<'a> WaveRayPath<'a> {
    /// Construct a new `WaveRayPath`
    ///
    /// # Arguments:
    ///
    /// `depth_data`: `&'a dyn BathymetryData`
    /// - a variable that implements the `BathymetryData` trait's `get_depth`
    ///   methods. Note that the lifetime requires that the `WaveRayPath` struct
    ///   will only live as long as the `BathymetryData` is available.
    ///
    /// `step_size`: `f64`
    /// - the change in time step used during Rk4 integration.
    ///
    /// Returns:
    /// `Self` : the newly created `WaveRayPath`
    pub fn new(depth_data: &'a dyn BathymetryData) -> Self {
        WaveRayPath {
            bathy_data: depth_data,
            current_data: None,
        }
    }

    pub fn new_with_current(
        depth_data: &'a dyn BathymetryData,
        current_data: Option<&'a dyn CurrentData>,
    ) -> Self {
        WaveRayPath {
            bathy_data: depth_data,
            current_data,
        }
    }

    /// build design method
    ///
    /// Used to create builder object then set each argument individually.
    ///
    /// # Returns
    /// `WaveRayPathBuilder<'a>` : the default WaveRayPathBuilder
    pub fn build() -> WaveRayPathBuilder<'a> {
        WaveRayPathBuilder::default()
    }

    /// Get the depth at the point x, y
    ///
    /// # Arguments
    /// `x` : `&f32`
    /// - the x point in meters
    ///
    /// `y` : `&f32`
    /// - the y point in meters
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : the depth at the x, y point
    /// - `Err(Error)` : an error while getting the depth
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : this error is returned when the
    /// `x` or `y` input give an out of bounds output.
    /// - `Error::InvalidArgument` : this error is returned from
    ///   `interpolator::bilinear` due to incorrect argument passed.
    pub fn depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
        let depth = self.bathy_data.get_depth(x, y)?;
        //println!("{:?}", depth);
        Ok(depth)
    }

    /// get the depth and gradient at point x, y
    pub fn depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error> {
        let h_dh = self.bathy_data.get_depth_and_gradient(x, y)?;
        //println!("The depth is: {}\nThe dh/dx is {}\nThe dh/dy is {}", h_dh.0, h_dh.1.0, h_dh.1.1);
        Ok(h_dh)
    }

    pub fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64), Error> {
        let (u, v) = self.current_data.unwrap().current(x, y)?;
        Ok((u, v))
    }

    /// Calculates system of odes from the given state
    ///
    /// The state is defined by x, y, kx, and ky. Then, the group velocity, depth,
    /// and depth gradient are calculated. The derivatives of the inputs are
    /// calculated using the equations in notes.md.
    ///
    /// # Arguments
    /// `x` : `&f64`
    /// - the x coordinate in meters
    ///
    /// `y` : `&f64`
    /// - the y coordinate in meters
    ///
    /// `kx` : `&f64`
    /// - x component of wavenumber vector
    ///
    /// `ky` : `&f64`
    /// - y component of wavenumber vector
    ///
    /// # Returns
    /// `Result<(f64, f64, f64, f64), Error>`
    /// - `Ok((f64, f64, f64, f64))` : a tuple of floats corresponding to (dxdt, dydt, dkxdt, dkydt).
    /// - `Err(Error)` : an error occured either getting the depth, or calculating the group velocity.
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : this error is returned when the
    /// `x` or `y` input give an out of bounds output.
    /// - `Error::InvalidArgument` : this error is returned from
    ///   `interpolator::bilinear` due to incorrect argument passed.
    /// `Error::ArgumentOutOfBounds`
    /// - If k is negative, group velocity will return this error.
    pub fn odes(
        &self,
        x: &f64,
        y: &f64,
        kx: &f64,
        ky: &f64,
    ) -> Result<(f64, f64, f64, f64), Error> {
        let (h, (dhdx, dhdy)) = self.depth_and_gradient(&(*x as f32), &(*y as f32))?;

        let (u, v, dudx, dvdx, dudy, dvdy) = if self.current_data.is_none() {
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        } else {
            let (a, b) = self.current(x, y)?;
            (a, b, 0.0, 0.0, 0.0, 0.0)
        };

        let h = h as f64;

        let k_mag = (kx * kx + ky * ky).sqrt();
        let k_dir = ky.atan2(*kx);

        let cg = group_velocity(&k_mag, &h)?;
        let cgx = cg * k_dir.cos() + u;
        let cgy = cg * k_dir.sin() + v;

        let dxdt = cgx;
        let dydt = cgy;

        let (dkxdt_bathy, dkydt_bathy) = dk_vector_dt(&k_mag, &h, &(dhdx as f64), &(dhdy as f64));

        let dkxdt = dkxdt_bathy - kx * dudx - ky * dvdx;
        let dkydt = dkydt_bathy - kx * dudy - ky * dvdy;

        Ok((dxdt, dydt, dkxdt, dkydt))
    }
}

/// Calculates the group velocity
///
/// # Arguments
///
/// `k` : `&f64`
/// - The wavenumber \[m^-1\] should always be positive.
///
/// `h` : `&f64`
/// - The depth \[m\] in this case should be positive.
///
/// # Returns
///
/// `Result<f64, Error>`
///
/// - `Ok(f64)` : returns the calculated group velocity as a float. Note: if `d`
///   is less then 0, it will return `f64::NAN`. In the future, this will return
///   either an error or warning.
///
/// - `Err(Error::ArgumentOutOfBounds)` : returns this error if k <= 0.
///
/// # Errors
///
/// `Error::ArgumentOutOfBounds`
/// - If k is negative, group velocity will return this error.
///
pub(crate) fn group_velocity(k: &f64, h: &f64) -> Result<f64, Error> {
    if *h <= 0.0 {
        return Ok(f64::NAN); // FIXME: should this also return an error?
    }
    if *k <= 0.0 {
        return Err(Error::ArgumentOutOfBounds);
    }
    let cg = (G / 2.0)
        * (((k * h).tanh() + (k * h) / (k * h).cosh().powi(2)) / (k * G * (k * h).tanh()).sqrt());
    // println!("The group velocity is: {}", cg);
    Ok(cg)
}

/// calculate the derivative of the wavenumber vector with respect to time
///
/// # Arguments
/// `k_mag` : `&f64`
/// - the magnitude of the wavenumber
///
/// `h` : &f64`
/// - the depth of the water
///
/// `dhdx` : `&f64`
/// - the partial of depth with respect to x
///
/// `dhdy` : `&f64`
/// - the partial of depth with respect to y
///
/// # Returns
/// `(f64, f64)` : values cooresponding to (dkx/dt, dky/dt)
fn dk_vector_dt(k_mag: &f64, h: &f64, dhdx: &f64, dhdy: &f64) -> (f64, f64) {
    // TODO: label this as dkdt bathy
    let dkxdt = (0.5) * k_mag * 1.0 / (k_mag * h).sinh() * 1.0 / (k_mag * h).cosh()
        * (G * k_mag * (k_mag * h).tanh()).sqrt()
        * dhdx;
    let dkydt = (0.5) * k_mag * 1.0 / (k_mag * h).sinh() * 1.0 / (k_mag * h).cosh()
        * (G * k_mag * (k_mag * h).tanh()).sqrt()
        * dhdy;

    //println!("The value for dkx/dt is {}", dkxdt);

    (dkxdt, dkydt)
}

type State = Vector4<f64>;
type Time = f64;

#[derive(Builder)]
/// A struct that stores the bathymetry/depth data related to an individual ray.
struct WaveRayPath<'a> {
    /// A reference to a pointer to a struct that implements the depth trait. The
    /// lifetime of WaveRayPath is restricted to the lifetime of this variable.
    bathy_data: &'a dyn BathymetryData,
    #[builder(default = "None")]
    current_data: Option<&'a dyn CurrentData>, // TODO: change this to implement a default builder
}

impl<'a> ode_solvers::System<State> for WaveRayPath<'a> {
    fn system(&self, t: Time, s: &State, ds: &mut State) {
        let (dxdt, dydt, dkxdt, dkydt) = match self.odes(&s[0], &s[1], &s[2], &s[3]) {
            Err(e) => {
                println!(
                    "ERROR at time {}: \"{}\". Setting all further output to NaN.",
                    t, e
                );
                (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
            }
            Ok(v) => v,
        };

        ds[0] = dxdt;
        ds[1] = dydt;
        ds[2] = dkxdt;
        ds[3] = dkydt;
    }
}

#[cfg(test)]
/// tests for constant depth
mod test_constant_bathymetry {
    use crate::{
        bathymetry::ArrayDepth, bathymetry::ConstantDepth, group_velocity, BathymetryData,
    };
    use crate::{dk_vector_dt, State, WaveRayPath};
    use ode_solvers::*;

    /// Runs ode solvers on the given check cases
    ///
    /// # Arguments
    /// `depth_data` : `Box<dyn depth>`
    /// - either ConstantDepth or ArrayDepth
    ///
    /// `check_axis` : `[(f64, f64, f64, f64); 4]`
    /// - these are an array of kx, ky, x_final, y_final
    ///
    /// # Panics
    /// If there is an error during integration of ode_solvers, this function will panic
    fn run_check_ode_solvers(
        depth_data: &dyn BathymetryData,
        check_axis: [(f64, f64, f64, f64); 4],
    ) {
        for (kx, ky, xf, yf) in check_axis {
            let system = WaveRayPath::new(depth_data);
            let y0 = State::new(0.0, 0.0, kx, ky);
            let mut stepper = Rk4::new(system, 0.0, y0, 1.0, 1.0);
            if stepper.integrate().is_ok() {
                let last_state = stepper.y_out().last().unwrap();
                assert!(
                    (last_state.x - xf).abs() < f64::EPSILON // super close, so I will take the values it gives as accurate
               && (last_state.y - yf).abs() < f64::EPSILON,
                    "expected xf: {}, actual: {} \nexpected yf: {}, actual: {}",
                    xf,
                    last_state.x,
                    yf,
                    last_state.y
                );
            } else {
                panic!("Error during ode_solvers integration")
            }
        }
    }

    #[test]
    /// testing group velocity function against values generated by wolfram alpha
    fn test_group_velocity() {
        let results = [
            (1.0, 1.565247584249853),
            (3.0, 0.9036961141150639),
            (5.0, 0.7),
            (10.0, 0.4949747468305833),
        ];
        for (k, ans) in results {
            assert!(
                (group_velocity(&k, &1000.0).unwrap() - ans).abs() < 1.0e-4,
                "k: {}, ans: {}",
                k,
                ans
            );
        }
    }

    #[test]
    /// verifying a negative k passed to group_velocity will return an error.
    fn test_negative_k() {
        assert!(group_velocity(&-1.0, &1000.0).is_err());
        assert!(group_velocity(&-12.0, &1000.0).is_err())
    }

    #[test]
    /// testing ode on simple cases worked out by hand
    fn test_odes() {
        let results = [
            // (kx, ky, dxdt, dydt)
            (1.0, 0.0, 1.565247584249853, 0.0),
            (0.0, 1.0, 0.0, 1.565247584249853),
            (-1.0, 0.0, -1.565247584249853, 0.0),
            (0.0, -1.0, 0.0, -1.565247584249853),
            // (0.0, 0.0, 0.0, 0.0) // this would cause panic
        ];

        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let system = WaveRayPath::new(data);

        for (kx, ky, ans_dxdt, ans_dydt) in results {
            let (dxdt, dydt, _, _) = system.odes(&0.0, &0.0, &kx, &ky).unwrap();
            assert!(
                (ans_dxdt - dxdt).abs() < 1.0e-4 && (ans_dydt - dydt).abs() < 1.0e-4,
                "ans_dxdt: {}, ans_dydt: {}, dxdt: {}, dydt: {}, kx: {}, ky: {}",
                ans_dxdt,
                ans_dydt,
                dxdt,
                dydt,
                kx,
                ky
            );
        }
    }

    #[test]
    /// all outputs should be NaN if k starts out of bounds
    fn test_zero_k() {
        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let system = WaveRayPath::new(data);
        let y0 = State::new(0.0, 0.0, 0.0, 0.0);

        let t0 = 0.0;
        let tf = 10.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        assert!(stepper.y_out().last().unwrap().x.is_nan());
        assert!(stepper.y_out().last().unwrap().y.is_nan());
        assert!(stepper.y_out().last().unwrap().z.is_nan());
        assert!(stepper.y_out().last().unwrap().w.is_nan());
    }

    #[test]
    fn test_zero_h() {
        let data: &dyn BathymetryData = &ConstantDepth::new(0.0);
        let system = WaveRayPath::new(data);
        let y0 = State::new(0.0, 0.0, 1.0, 1.0);

        let t0 = 0.0;
        let tf = 10.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        assert!(stepper.y_out().last().unwrap().x.is_nan());
        assert!(stepper.y_out().last().unwrap().y.is_nan());
        assert!(stepper.y_out().last().unwrap().z.is_nan());
        assert!(stepper.y_out().last().unwrap().w.is_nan());
    }

    #[test]
    /// Testing the ode_solvers Rk4 function only in the kx or ky direction
    fn test_axis() {
        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        // answers should be the square root of gravity
        let check_axis = [
            (0.0, 1.0, 0.0, (9.8_f64).sqrt() / 2.0),
            (1.0, 0.0, (9.8_f64).sqrt() / 2.0, 0.0),
            (0.0, -1.0, 0.0, -(9.8_f64).sqrt() / 2.0),
            (-1.0, 0.0, -(9.8_f64).sqrt() / 2.0, 0.0),
        ];

        run_check_ode_solvers(data, check_axis)
    }

    #[test]
    /// check that function accepts array
    fn test_array_as_parameter() {
        let data: &dyn BathymetryData = &ArrayDepth::new(vec![
            vec![1000.0, 1000.0, 1000.0],
            vec![1000.0, 1000.0, 1000.0],
            vec![1000.0, 1000.0, 1000.0],
        ]);
        // answers should be the square root of gravity divided by 2
        let check_axis = [
            (0.0, 1.0, 0.0, (9.8_f64).sqrt() / 2.0),
            (1.0, 0.0, (9.8_f64).sqrt() / 2.0, 0.0),
            (0.0, -1.0, 0.0, -(9.8_f64).sqrt() / 2.0),
            (-1.0, 0.0, -(9.8_f64).sqrt() / 2.0, 0.0),
        ];

        run_check_ode_solvers(data, check_axis)
    }

    #[test]
    /// if x input is NAN, the output x should be NaN. if k is zero, it will still error.
    fn test_x_nan() {
        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let system = WaveRayPath::new(data);
        let nan = f64::NAN;
        let y0 = State::new(nan, 0.0, 1.0, 0.0);

        let t0 = 0.0;
        let tf = 1.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        assert!(stepper.y_out().last().unwrap().x.is_nan());
    }

    #[test]
    /// if y input is NAN, the output x should be NaN. if k is zero, it will still error.
    fn test_y_nan() {
        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let system = WaveRayPath::new(data);
        let nan = f64::NAN;
        let y0 = State::new(0.0, nan, 1.0, 0.0);

        let t0 = 0.0;
        let tf = 1.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        assert!(stepper.y_out().last().unwrap().y.is_nan());
    }

    #[test]
    /// if either k input is NAN, the output x and y should be NaN.
    fn test_kx_nan() {
        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let system = WaveRayPath::new(data);
        let nan = f64::NAN;
        let y0 = State::new(0.0, 0.0, nan, 0.0);

        let t0 = 0.0;
        let tf = 1.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        assert!(stepper.y_out().last().unwrap().x.is_nan());
        assert!(stepper.y_out().last().unwrap().y.is_nan());
    }

    #[test]
    /// if either k input is NAN, the output x and y should be NaN.
    fn test_ky_nan() {
        let data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let system = WaveRayPath::new(data);
        let nan = f64::NAN;
        let y0 = State::new(0.0, 0.0, 0.0, nan);

        let t0 = 0.0;
        let tf = 1.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        assert!(stepper.y_out().last().unwrap().x.is_nan());
        assert!(stepper.y_out().last().unwrap().y.is_nan());
    }

    #[test]
    /// test when d / wavelenth < 1 / 20
    fn test_shallow() {
        let data: &dyn BathymetryData = &ConstantDepth::new(0.1);
        // the approximation is the square root of gravity * h, but are not, they get closer as d approaches 0.
        let check_axis = [
            // the numbers very close to zero are likely due to switching between f32 and f64
            (
                0.0,
                1.0,
                0.00000000000000006031543168844801,
                0.9850257515953494,
            ), // should be 0.0, 0.9899494936611665
            (1.0, 0.0, 0.9850257515953494, 0.0),
            (
                0.0,
                -1.0,
                0.00000000000000006031543168844801,
                -0.9850257515953494,
            ),
            (
                -1.0,
                0.0,
                -0.9850257515953494,
                0.00000000000000012063086337689602,
            ),
        ];

        run_check_ode_solvers(data, check_axis)
    }

    #[test]
    /// If the bathymetry array index is out of range, it will return nan.
    fn out_of_range_give_nan() {
        let data: &dyn BathymetryData =
            &ArrayDepth::new(vec![vec![1000.0, 1000.0], vec![1000.0, 1000.0]]);
        let system = WaveRayPath::new(data);
        let y0 = State::new(0.0, 0.0, 0.0, 1.0);

        let t0 = 0.0;
        let tf = 10.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        let last_step = stepper.y_out().last().unwrap();

        assert!(last_step.x.is_nan() && last_step.y.is_nan());
    }

    #[test]
    // test the k derivative function
    fn test_dk_deep() {
        let k_mag = 1000.0;
        let h = 1000.0;
        let dhdx = 0.2;
        let dhdy = 0.2;

        let ans = dk_vector_dt(&k_mag, &h, &dhdx, &dhdy);

        assert!(
            (ans.0 - 0.0).abs() < f64::EPSILON,
            "Expected 0, got {}",
            ans.0
        );
        assert!(
            (ans.1 - 0.0).abs() < f64::EPSILON,
            "Expected 0, got {}",
            ans.1
        )
    }
}

/// tests for constant current
#[cfg(test)]
mod test_current {
    use crate::{
        bathymetry::{BathymetryData, ConstantDepth},
        current::{ConstantCurrent, CurrentData},
        WaveRayPath,
    };

    #[test]
    /// this test I added by copying a test from the module
    /// test_constant_current and using the WaveRayPath from the builder. I am
    /// comparing the results using the function `odes` because it uses both the
    /// current and the bathymetry and will make sure both work.
    fn test_waveraypath_builder() {
        let bd = ConstantDepth::new(1000.0);
        let cd = ConstantCurrent::new(0.0, 0.0);

        // build pattern with supplying current data
        let wave = WaveRayPath::build()
            .bathy_data(&bd) // TODO: should bathy_data also be optional?
            .current_data(Some(&cd))
            .build()
            .unwrap();

        // build pattern without supplying current data
        let wave2 = WaveRayPath::build().bathy_data(&bd).build().unwrap();

        let results = [
            // (kx, ky, dxdt, dydt)
            (1.0, 0.0, 1.565247584249853, 0.0),
            (0.0, 1.0, 0.0, 1.565247584249853),
            (-1.0, 0.0, -1.565247584249853, 0.0),
            (0.0, -1.0, 0.0, -1.565247584249853),
            // (0.0, 0.0, 0.0, 0.0) // this would cause panic
        ];

        // check the first wave
        for (kx, ky, ans_dxdt, ans_dydt) in results {
            let (dxdt, dydt, _, _) = wave.odes(&0.0, &0.0, &kx, &ky).unwrap();
            assert!(
                (ans_dxdt - dxdt).abs() < 1.0e-4 && (ans_dydt - dydt).abs() < 1.0e-4,
                "ans_dxdt: {}, ans_dydt: {}, dxdt: {}, dydt: {}, kx: {}, ky: {}",
                ans_dxdt,
                ans_dydt,
                dxdt,
                dydt,
                kx,
                ky
            );
        }

        // check the second wave
        for (kx, ky, ans_dxdt, ans_dydt) in results {
            let (dxdt, dydt, _, _) = wave2.odes(&0.0, &0.0, &kx, &ky).unwrap();
            assert!(
                (ans_dxdt - dxdt).abs() < 1.0e-4 && (ans_dydt - dydt).abs() < 1.0e-4,
                "ans_dxdt: {}, ans_dydt: {}, dxdt: {}, dydt: {}, kx: {}, ky: {}",
                ans_dxdt,
                ans_dydt,
                dxdt,
                dydt,
                kx,
                ky
            );
        }
    }

    #[test]
    fn test_constant_depth_current() {
        // test case: initial group velocity in x axis only test 1, -1 for both
        // u and v individually the results should be equal to the original
        // without current plus or minus one. dk/dt will be zero in both x and y

        // these results are copied from test_odes in mod test_constant_bathymetry, but with 1 added in the correct place
        let results = [
            // (kx, ky, dxdt, dydt)
            (1.0, 0.0, 1.565247584249853 + 1.0, 0.0), // u = 1, v = 0
            (1.0, 0.0, 1.565247584249853 - 1.0, 0.0), // u = -1, v = 0
            (1.0, 0.0, 1.565247584249853, 0.0 + 1.0), // u = 0, v = 1
            (1.0, 0.0, 1.565247584249853, 0.0 - 1.0), // u = 0, v = -1
        ];

        let bathy_data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let current_data_1: &dyn CurrentData = &ConstantCurrent::new(1.0, 0.0);
        let current_data_2: &dyn CurrentData = &ConstantCurrent::new(-1.0, 0.0);
        let current_data_3: &dyn CurrentData = &ConstantCurrent::new(0.0, 1.0);
        let current_data_4: &dyn CurrentData = &ConstantCurrent::new(0.0, -1.0);

        // TODO: need a default builder

        for (i, (kx, ky, ans_dxdt, ans_dydt)) in results.iter().enumerate() {
            let system = match i {
                0 => WaveRayPath::new_with_current(bathy_data, Some(current_data_1)),
                1 => WaveRayPath::new_with_current(bathy_data, Some(current_data_2)),
                2 => WaveRayPath::new_with_current(bathy_data, Some(current_data_3)),
                3 => WaveRayPath::new_with_current(bathy_data, Some(current_data_4)),
                _ => WaveRayPath::new_with_current(bathy_data, None),
            };
            let (dxdt, dydt, _, _) = system.odes(&0.0, &0.0, &kx, &ky).unwrap();
            assert!(
                (ans_dxdt - dxdt).abs() < f64::EPSILON && (ans_dydt - dydt).abs() < f64::EPSILON,
                "ans_dxdt: {}, ans_dydt: {}, dxdt: {}, dydt: {}, kx: {}, ky: {}",
                ans_dxdt,
                ans_dydt,
                dxdt,
                dydt,
                kx,
                ky
            );
        }
    }
}
