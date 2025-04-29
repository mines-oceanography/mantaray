//! Struct and methods for calculating the ray path of a wave.
//!
//! The `WaveRayPath` struct contains references to the depth and current (or
//! none if using defaults), and contains methods to calculate the change in
//! state using the equations from `notes.md`.
//!
//! The `WaveRayPath` struct implements `ode_solvers`'s `system` and `solout`
//! functions to use RK4 to numerically integrate these equations given the
//! initial conditions.

use derive_builder::Builder;
use ode_solvers::*;

use crate::bathymetry::BathymetryData;
use crate::bathymetry::DEFAULT_BATHYMETRY;
use crate::current::CurrentData;
use crate::current::DEFAULT_CURRENT;
use crate::datatype::Point;
use crate::error::Error;
use crate::error::Result;

/// constant for gravity
const G: f64 = 9.8;

/// state of the ray system for `ode_solvers`
/// the values in the state are x, y, kx, ky
/// for example: `State::new(x, y, kx, ky)`
pub(crate) type State = Vector4<f64>;

/// time in seconds for `ode_solvers` to use
pub(crate) type Time = f64;

#[derive(Builder)]
/// Stores the bathymetry and current data and calculates the system of odes
/// that define the ray tracing.
///
/// This struct can be created in two ways:
/// - using the `new` method, for example, `WaveRayPath::new(Some(&depth_data),
///   Some(&current_data))`
/// - using the `builder` method, for example,
///   `WaveRayPath::builder().bathymetry_data(&depth_data).current_data(&current_data).build().unwrap()`
pub(crate) struct WaveRayPath<'a> {
    #[builder(default = "&DEFAULT_BATHYMETRY")]
    /// A reference to a BathymetryData trait object. If this is None, the depth
    /// will be set to 2000 m.
    bathymetry_data: &'a dyn BathymetryData,
    #[builder(default = "&DEFAULT_CURRENT")]
    /// Optional reference to a CurrentData trait object. If this is None, the
    /// current will be set to 0 m/s.
    current_data: &'a dyn CurrentData,
}

#[allow(dead_code)]
impl<'a> WaveRayPath<'a> {
    /// Construct a new `WaveRayPath`
    ///
    /// # Arguments:
    ///
    /// `depth_data`: `Option<&'a dyn BathymetryData>`
    /// - an optional variable that implements the `BathymetryData` trait's
    ///   `depth` methods. If this is `None`, the depth will be set to 2000
    ///   m.
    ///
    /// `current_data`: `Option<&'a dyn CurrentData>`
    /// - an optional variable that implements the `CurrentData` trait's
    ///  `current_and_gradient` method. If this is `None`, the current will be
    ///  set to 0 m/s.
    ///
    /// Returns: `Self` : the newly created `WaveRayPath`
    pub(crate) fn new(bathymetry_data: &'a dyn BathymetryData, current_data: &'a dyn CurrentData) -> Self {
        WaveRayPath {
            bathymetry_data,
            current_data,
        }
    }

    /// create a new `WaveRayPathBuilder` using the builder method
    ///
    /// Used to create builder object then set each argument individually. For
    /// example,
    /// `WaveRayPath::builder().bathymetry_data(&depth_data).current_data(&current_data).build().unwrap()`
    ///
    /// # Returns
    /// `WaveRayPathBuilder<'a>` : the default WaveRayPathBuilder
    pub(crate) fn builder() -> WaveRayPathBuilder<'a> {
        WaveRayPathBuilder::default()
    }

    /// Calculates system of odes from the given state
    ///
    /// See `notes.md` for more information
    ///
    /// # Arguments
    /// `x` : `&f64`
    /// - the x coordinate in meters
    ///
    /// `y` : `&f64`
    /// - the y coordinate in meters
    ///
    /// `kx` : `&f64`
    /// - x component of wavenumber vector \[m^-1\]
    ///
    /// `ky` : `&f64`
    /// - y component of wavenumber vector \[m^-1\]
    ///
    /// # Returns
    /// `Result<(f64, f64, f64, f64)>`
    /// - `Ok((f64, f64, f64, f64))` : a tuple of floats corresponding to (dxdt, dydt, dkxdt, dkydt).
    /// - `Err(Error)` : an error occurred either getting the depth, or calculating the group velocity.
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : this error is returned when the `x` or `y`
    /// input give an out of bounds output.
    /// - `Error::InvalidArgument` : this error is returned from
    ///   `interpolator::bilinear` due to incorrect argument passed.
    /// `Error::ArgumentOutOfBounds`
    /// - If k is negative, group velocity will return this error.
    pub(crate) fn odes(&self, x: &f64, y: &f64, kx: &f64, ky: &f64) -> Result<(f64, f64, f64, f64)> {
        let point = crate::Point::new(*x, *y);
        let (h, dh) = self
            .bathymetry_data
            .depth_and_gradient(&Point::new(*x as f32, *y as f32))?;

        let h = h as f64;
        let dhdx = *dh.dx() as f64;
        let dhdy = *dh.dy() as f64;

        // get the current and gradient from the current data or use default.
        let (current, (du, dv)) = self.current_data.current_and_gradient(&point)?;

        // magnitude and direction of the wavenumber
        let k = (kx * kx + ky * ky).sqrt();
        let theta = ky.atan2(*kx);

        // calculate the group velocity
        let cg = self.group_velocity(&k, &h)?;
        let cgx = cg * theta.cos() + current.u();
        let cgy = cg * theta.sin() + current.v();

        let dxdt = cgx;
        let dydt = cgy;

        // calculate dk/dt
        let (dkxdt_bathy, dkydt_bathy) = self.dkdt_bathy(&k, &h, &dhdx, &dhdy);

        let dkxdt = dkxdt_bathy - kx * du.dx() - ky * dv.dx();
        let dkydt = dkydt_bathy - kx * du.dy() - ky * dv.dy();

        Ok((dxdt, dydt, dkxdt, dkydt))
    }

    /// Calculates the group velocity
    ///
    /// # Arguments
    ///
    /// `k` : `&f64`
    /// - The magnitude of the wavenumber \[m^-1\] should always be positive.
    ///
    /// `h` : `&f64`
    /// - The depth \[m\] in this case should be positive.
    ///
    /// # Returns
    ///
    /// `Result<f64>`
    ///
    /// - `Ok(f64)` : returns the calculated group velocity as a float. Note: if
    ///   `h` is less then 0, it will return `f64::NAN`. This represents the
    ///   wave have no water to propagate in.
    ///
    /// - `Err(Error::ArgumentOutOfBounds)` : returns this error if k <= 0.
    ///
    /// # Errors
    ///
    /// `Error::ArgumentOutOfBounds`
    /// - If k is negative, group velocity will return this error.
    ///
    pub(crate) fn group_velocity(&self, k: &f64, h: &f64) -> Result<f64> {
        if *h <= 0.0 {
            return Ok(f64::NAN);
        }
        if *k <= 0.0 {
            return Err(Error::ArgumentOutOfBounds);
        }
        let cg = (G / 2.0)
            * (((k * h).tanh() + (k * h) / (k * h).cosh().powi(2))
                / (k * G * (k * h).tanh()).sqrt());
        Ok(cg)
    }

    /// calculate the derivative of the wavenumber vector with respect to time
    ///
    /// # Arguments
    /// `k` : `&f64`
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
    /// `(f64, f64)` : values corresponding to (dkx/dt, dky/dt)
    fn dkdt_bathy(&self, k: &f64, h: &f64, dhdx: &f64, dhdy: &f64) -> (f64, f64) {
        let dkxdt_bathy = (-0.5) * k * 1.0 / (k * h).sinh() * 1.0 / (k * h).cosh()
            * (G * k * (k * h).tanh()).sqrt()
            * dhdx;
        let dkydt_bathy = (-0.5) * k * 1.0 / (k * h).sinh() * 1.0 / (k * h).cosh()
            * (G * k * (k * h).tanh()).sqrt()
            * dhdy;

        (dkxdt_bathy, dkydt_bathy)
    }
}

impl<'a> ode_solvers::System<Time, State> for WaveRayPath<'a> {
    fn system(&self, _t: Time, s: &State, ds: &mut State) {
        // calculate the derivatives using the system of odes
        let (dxdt, dydt, dkxdt, dkydt) = match self.odes(&s[0], &s[1], &s[2], &s[3]) {
            Err(_) => {
                // Error at time t. Setting all further output to NaN.
                (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
            }
            Ok(v) => v,
        };

        ds[0] = dxdt;
        ds[1] = dydt;
        ds[2] = dkxdt;
        ds[3] = dkydt;
    }

    fn solout(&mut self, _x: Time, y: &State, dy: &State) -> bool {
        if (dy[0].is_nan() && dy[1].is_nan() && dy[2].is_nan() && dy[3].is_nan())
            || (y[0].is_nan() && y[1].is_nan() && y[2].is_nan() && y[3].is_nan())
        {
            // NaN in derivatives or output. Likely reached end of current or
            // bathy domain. Stopping integration.
            true
        } else {
            false
        }
    }
}

#[cfg(test)]
/// tests for constant depth
mod test_constant_bathymetry {
    use crate::current::ConstantCurrent;
    use crate::wave_ray_path::{State, WaveRayPath};
    use crate::{bathymetry::ArrayDepth, bathymetry::BathymetryData, bathymetry::ConstantDepth};
    use ode_solvers::*;

    /// Runs ode solvers on the given check cases
    ///
    /// # Arguments
    /// `depth_data` : `Box<dyn depth>`
    /// - either ConstantDepth or ArrayDepth
    ///
    /// `check_axis` : `[(f64, f64, f64, f64); 4]`
    /// - these are an array of kx, ky, final x, and final y
    ///
    /// # Panics
    /// If there is an error during integration of ode_solvers, this function will panic
    fn run_check_ode_solvers(
        depth_data: &dyn BathymetryData,
        check_axis: [(f64, f64, f64, f64); 4],
    ) {
        for (kx, ky, xf, yf) in check_axis {
            let current = ConstantCurrent::new(0.0, 0.0);
            let system = WaveRayPath::new(depth_data, &current);
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
        let depth = ConstantDepth::new(1000.0);
        let current = ConstantCurrent::new(0.0, 0.0);

        let wave_ray_path = WaveRayPath::new(&depth, &current);
        let results = [
            (1.0, 1.565247584249853),
            (3.0, 0.9036961141150639),
            (5.0, 0.7),
            (10.0, 0.4949747468305833),
        ];
        for (k, ans) in results {
            assert!(
                (wave_ray_path.group_velocity(&k, &1000.0).unwrap() - ans).abs() < 1.0e-4,
                "k: {}, ans: {}",
                k,
                ans
            );
        }
    }

    #[test]
    /// verifying a negative k passed to group_velocity will return an error.
    fn test_negative_k() {
        let depth = ConstantDepth::new(1000.0);
        let current = ConstantCurrent::new(0.0, 0.0);

        let wave_ray_path = WaveRayPath::new(&depth, &current);
        assert!(wave_ray_path.group_velocity(&-1.0, &1000.0).is_err());
        assert!(wave_ray_path.group_velocity(&-12.0, &1000.0).is_err())
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);

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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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
    /// test when d / wavelength < 1 / 20
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
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
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

        let depth = ConstantDepth::new(1000.0);
        let current = ConstantCurrent::new(0.0, 0.0);

        let wave_ray_path = WaveRayPath::new(&depth, &current);

        let ans = wave_ray_path.dkdt_bathy(&k_mag, &h, &dhdx, &dhdy);

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

    #[test]
    // test the solout function stops integration early
    fn test_solout() {
        let data: &dyn BathymetryData =
            &ArrayDepth::new(vec![vec![1000.0, 1000.0], vec![1000.0, 1000.0]]);
        let current = ConstantCurrent::new(0.0, 0.0);

        let system = WaveRayPath::new(data, &current);
        let y0 = State::new(0.0, 0.0, 0.0, 1.0);

        let t0 = 0.0;
        let tf = 10.0;
        let step_size = 1.0;

        let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
        let _ = stepper.integrate();

        // the integration stops much before the final time (10s)
        assert_eq!(*(&stepper.results().get().0.len()), 3);

        // the last time stamp should de NAN
        let last_step = stepper.y_out().last().unwrap();

        assert!(last_step.x.is_nan());
        assert!(last_step.y.is_nan());
        assert!(last_step.z.is_nan());
        assert!(last_step.w.is_nan());
    }
}

/// tests for constant current
#[cfg(test)]
mod test_current {
    use crate::{
        bathymetry::{BathymetryData, ConstantDepth},
        current::{ConstantCurrent, CurrentData},
        wave_ray_path::WaveRayPath,
    };

    #[test]
    /// this test I added by copying a test from the module
    /// test_constant_current and using the WaveRayPath from the builder. I am
    /// comparing the results using the function `odes` because it uses both the
    /// current and the bathymetry and will make sure both work.
    fn test_wave_ray_path_builder() {
        let bd = ConstantDepth::new(1000.0);
        let cd = ConstantCurrent::new(0.0, 0.0);

        // build pattern with supplying current data
        let wave = WaveRayPath::builder()
            .bathymetry_data(&bd)
            .current_data(&cd)
            .build()
            .unwrap();

        // build pattern without supplying current data
        let wave2 = WaveRayPath::builder().bathymetry_data(&bd).build().unwrap();

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
            (1.0, 0.0, 1.565247584249853 + 1.0, 0.0 + 1.0), // u = 1, v = 1
            (1.0, 0.0, 1.565247584249853 - 1.0, 0.0 - 1.0), // u = -1, v = -1
        ];

        let bathy_data: &dyn BathymetryData = &ConstantDepth::new(1000.0);
        let current_data_1: &dyn CurrentData = &ConstantCurrent::new(1.0, 0.0);
        let current_data_2: &dyn CurrentData = &ConstantCurrent::new(-1.0, 0.0);
        let current_data_3: &dyn CurrentData = &ConstantCurrent::new(0.0, 1.0);
        let current_data_4: &dyn CurrentData = &ConstantCurrent::new(0.0, -1.0);
        let current_data_5: &dyn CurrentData = &ConstantCurrent::new(1.0, 1.0);
        let current_data_6: &dyn CurrentData = &ConstantCurrent::new(-1.0, -1.0);

        for (i, (kx, ky, ans_dxdt, ans_dydt)) in results.iter().enumerate() {
            let system = match i {
                0 => WaveRayPath::builder()
                    .bathymetry_data(bathy_data)
                    .current_data(current_data_1)
                    .build()
                    .unwrap(),
                1 => WaveRayPath::builder()
                    .bathymetry_data(bathy_data)
                    .current_data(current_data_2)
                    .build()
                    .unwrap(),
                2 => WaveRayPath::builder()
                    .bathymetry_data(bathy_data)
                    .current_data(current_data_3)
                    .build()
                    .unwrap(),
                3 => WaveRayPath::builder()
                    .bathymetry_data(bathy_data)
                    .current_data(current_data_4)
                    .build()
                    .unwrap(),
                4 => WaveRayPath::builder()
                    .bathymetry_data(bathy_data)
                    .current_data(current_data_5)
                    .build()
                    .unwrap(),
                5 => WaveRayPath::builder()
                    .bathymetry_data(bathy_data)
                    .current_data(current_data_6)
                    .build()
                    .unwrap(),
                _ => panic!("Index out of range"),
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
