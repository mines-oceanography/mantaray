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
use ode_solvers::*;
use std::io::Write;
use std::fs::File;

mod error;

mod interpolator;

mod bathymetry;

mod ray;

use error::Error;

// Define constant gravity
const G: f64 = 9.8;


/// Run the Rk4 stepper.integrate() and save the results to a space separated file
/// 
/// # Arguments
/// `system` : `WaveRayPath`
/// - struct that calculates the derivatives for each state
/// 
/// `t0` : `f64`
/// - initial time
/// 
/// `y0` : `State`
/// - initial state
/// 
/// `tf` : `f64`
/// - final time
/// 
/// `step_size` : `f64`
/// - size of time increment delta t
/// 
/// # Returns
/// `Result<File, Error>`
/// - `Ok(File)` : return the file created
/// - `Err(Error::Undefined)` : error during integration
/// 
/// # Errors
/// `Error::Undefined` : this is a placeholder for an integration error during the Rk4 algorithm.
/// 
/// # Panics
/// The function will panic if it can not read or write to the file.
fn output_to_file(
    system: WaveRayPath,
    t0: f64,
    y0: State,
    tf: f64,
    step_size: f64,
) -> Result<File, Error> {
    let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
    stepper.integrate()?;

    let y = stepper.y_out();

    let mut file = File::create("y_out.txt")?;
    writeln!(&mut file, "t x y kx ky")?;
    for (i, x) in stepper.x_out().iter().enumerate() {
        write!(&mut file, "{x} ")?;
        for elem in y[i].iter() {
            write!(&mut file, "{elem} ")?;
        }
        writeln!(&mut file, " ")?;
    }
    Ok(file)
}

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
   pub fn new(depth_data: &'a dyn BathymetryData, step_size: f64) -> Self {
      WaveRayPath { data: depth_data, step_size: step_size }
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
      let depth = self.data.get_depth(x, y)?;
      println!("{:?}", depth);
      Ok(depth)
   }

   /// calculates the depth gradient at a given x, y using finite differences.
   /// 
   /// # Arguments:
   /// `x`: `&f32`
   /// -  x-coordinate of the point
   /// 
   /// `y`: `&f32`
   /// - y-coordinate of the point
   /// 
   /// `dx`: `&f32`
   /// - change in the x direction
   /// 
   /// `dy`: `&f32`
   /// - change in the y direction
   /// 
   /// # Returns:
   /// `Result<(f32, f32), Error`
   /// - `Ok((f32,f32))` : values representing the depth gradient (dh/dx, dh/dy)
   /// -`Err(Error)` : there was an error when getting the depth
   /// 
   /// # Errors
   /// - `Error::IndexOutOfBounds` : this error is returned when the `x` or `y`
   /// input give an out of bounds output.
   /// - `Error::InvalidArgument` : this error is returned from
   ///   `interpolator::bilinear` due to incorrect argument passed.
   /// 
   /// # Note
   /// It would be more efficient to calculate the gradient while calculating
   /// the depth in the first place. That way the data is only looped over once.
   pub fn gradient(&self, x: &f32, y: &f32, dx: &f32, dy: &f32) -> Result<(f32, f32), Error> {

      let x_grad = if *dx==0.0 {
         0.0
      } else {
         (self.data.get_depth(&(x + dx), y)? - self.data.get_depth(&(x - dx), y)?) / (2.0 * dx)
      };

      let y_grad = if *dy==0.0 {
         0.0
      } else {
         (self.data.get_depth(x, &(y + dy))? - self.data.get_depth(x, &(y - dy))?) / (2.0 * dy)
      };

      Ok((x_grad, y_grad))

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
   pub fn odes(&self, x: &f64, y: &f64, kx: &f64, ky: &f64) -> Result<(f64, f64, f64, f64), Error> {

      let h = self.depth(&(*x as f32), &(*y as f32))? as f64;

      let k_mag = (kx*kx + ky*ky).sqrt();
      let k_dir = ky.atan2(*kx);

      let cg = group_velocity(&k_mag, &h)?;
      let cgx = cg * k_dir.cos();
      let cgy = cg * k_dir.sin();

      let dxdt = cgx;
      let dydt = cgy;

      let dx = (dxdt * self.step_size).abs();
      let dy = (dydt * self.step_size).abs();

      let grad_of_h = self.gradient(&(*x as f32), &(*y as f32), &(dx as f32), &(dy as f32))?;

      let (dkxdt, dkydt) = dk_vector_dt(&k_mag, &h, &(grad_of_h.0 as f64), &(grad_of_h.1 as f64));

      Ok((dxdt, dydt, dkxdt, dkydt))
   }

}

struct ConstantDepth {
   h: f32,
}

impl BathymetryData for ConstantDepth {
   fn get_depth(&self, _x: &f32, _y: &f32) -> Result<f32, Error> {
       Ok(self.h)
   }
}

struct ArrayDepth {
   array: Vec<Vec<f32>>,
}

// FIXME: the program is not crashing, but NAN isn't that useful. maybe use option: some or none
impl BathymetryData for ArrayDepth {
   fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
      if *x as usize >= self.array.len() || *y as usize >= self.array.len() {
         return Ok(f32::NAN);
      }
       Ok(self.array[*x as usize][*y as usize]) // FIXME: since x and y are floats, they are truncated or rounded to a usize. I probably want a better interpolation estimate
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
   if *h < 0.0 {
      return Ok(f64::NAN); // FIXME: should this also return an error?
   }
   if *k <= 0.0 {
      return Err(Error::ArgumentOutOfBounds);
   }
   Ok( (G / 2.0) * ( ((k*h).tanh() + (k*h)/(k*h).cosh().powi(2)) / (k*G*(k*h).tanh()).sqrt() ) )
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
   let dkxdt = (-0.5 * G * k_mag / ( k_mag * h ).tanh().sqrt()) * dhdx;
   let dkydt = (-0.5 * G * k_mag / ( k_mag * h ).tanh().sqrt()) * dhdy;

   (dkxdt, dkydt)
}


type State = Vector4<f64>;
type Time = f64;

/// A struct that stores the bathymetry/depth data related to an individual ray.
struct WaveRayPath<'a> {
   /// A reference to a pointer to a struct that implements the depth trait. The
   /// lifetime of WaveRayPath is restricted to the lifetime of this variable.
   data: &'a dyn BathymetryData,
   /// the size of delta t used during Rk4 integration.
   step_size: f64,
}

impl<'a> ode_solvers::System<State> for WaveRayPath<'a> {
   fn system(&self, t: Time, s: &State, ds: &mut State) { 

       let (dxdt, dydt, dkxdt, dkydt) = 
       match self.odes(&s[0], &s[1], &s[2], &s[3]) {
         Err(e) => {
            println!("ERROR at time {}: \"{}\". Setting all further output to NaN.", t, e);
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
/// tests for constant depth, constant group velocity
mod test_constant_cg {
   use crate::{State, WaveRayPath, dk_vector_dt};
   use ode_solvers::*;
   use crate::{group_velocity, ConstantDepth, ArrayDepth, output_to_file, BathymetryData};

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
   fn run_check_ode_solvers(depth_data: &dyn BathymetryData, check_axis: [(f64, f64, f64, f64); 4]) {
      for (kx, ky, xf, yf) in check_axis {
         let system = WaveRayPath::new(depth_data, 1.0);
         let y0 = State::new(0.0, 0.0, kx, ky);
         let mut stepper = Rk4::new(system, 0.0, y0, 1.0, 1.0);
         if stepper.integrate().is_ok() {
            let last_state = stepper.y_out().last().unwrap();
            assert!(
               (last_state.x - xf).abs() < f64::EPSILON // super close, so I will take the values it gives as accurate
               && (last_state.y - yf).abs() < f64::EPSILON,
               "expected xf: {}, actual: {} \nexpected yf: {}, actual: {}",
               xf, last_state.x, yf, last_state.y
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
         (10.0, 0.4949747468305833)
      ];
      for (k, ans) in results {
         assert!((group_velocity(&k, &1000.0).unwrap() - ans).abs() < 1.0e-4, "k: {}, ans: {}", k, ans);
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

      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);

      for (kx, ky, ans_dxdt, ans_dydt) in results {
         let (dxdt, dydt, _, _) = system.odes(&0.0, &0.0, &kx, &ky).unwrap();
         assert!(
            (ans_dxdt - dxdt).abs() < 1.0e-4
            && (ans_dydt -dydt).abs() < 1.0e-4,
            "ans_dxdt: {}, ans_dydt: {}, dxdt: {}, dydt: {}, kx: {}, ky: {}",
            ans_dxdt, ans_dydt, dxdt, dydt, kx, ky
         );
      }
   }

   #[test]
   /// all outputs should be NaN if k starts out of bounds
   fn test_zero_k() {
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);
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
   /// Testing the ode_solvers Rk4 function only in the kx or ky direction
   fn test_axis() {
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      // answers should be the square root of gravity
      let check_axis = [
         (0.0, 1.0, 0.0, (9.8_f64).sqrt()/2.0),
         (1.0, 0.0, (9.8_f64).sqrt()/2.0, 0.0),
         (0.0, -1.0, 0.0, -(9.8_f64).sqrt()/2.0),
         (-1.0, 0.0, -(9.8_f64).sqrt()/2.0, 0.0)
      ];

      run_check_ode_solvers(data, check_axis)
   }

   #[test]
   /// check that function accepts array
   fn test_array_as_parameter() {
      let data: &dyn BathymetryData = &ArrayDepth { array: vec![
         vec![1000.0, 1000.0, 1000.0],
         vec![1000.0, 1000.0, 1000.0],
         vec![1000.0, 1000.0, 1000.0]
      ] };
      // answers should be the square root of gravity divided by 2
      let check_axis = [
         (0.0, 1.0, 0.0, (9.8_f64).sqrt()/2.0),
         (1.0, 0.0, (9.8_f64).sqrt()/2.0, 0.0),
         (0.0, -1.0, 0.0, -(9.8_f64).sqrt()/2.0),
         (-1.0, 0.0, -(9.8_f64).sqrt()/2.0, 0.0)
      ];

         run_check_ode_solvers(data, check_axis)
   }

   #[test]
   /// if x input is NAN, the output x should be NaN. if k is zero, it will still error.
   fn test_x_nan() {
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);
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
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);
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
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);
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
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);
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
      let data : &dyn BathymetryData = &ConstantDepth { h: 0.1 };
      // the approximation is the square root of gravity * h, but are not, they get closer as d approaches 0.
      let check_axis = [
         // the numbers very close to zero are likely due to switching between f32 and f64
         (0.0, 1.0, 0.00000000000000006031543168844801, 0.9850257515953494), // should be 0.0, 0.9899494936611665
         (1.0, 0.0, 0.9850257515953494, 0.0),
         (0.0, -1.0, 0.00000000000000006031543168844801, -0.9850257515953494),
         (-1.0, 0.0, -0.9850257515953494, 0.00000000000000012063086337689602)
      ];

      run_check_ode_solvers(data, check_axis)
   }

   #[test]
   /// If the bathymetry array index is out of range, it will return nan.
   fn out_of_range_give_nan() {
      let data: &dyn BathymetryData = &ArrayDepth { array: vec![
         vec![1000.0, 1000.0],
         vec![1000.0, 1000.0]
      ] };
      let system = WaveRayPath::new(data, 1.0);
      let y0 = State::new(0.0, 0.0, 0.0, 1.0);
   
      let t0 = 0.0;
      let tf = 10.0;
      let step_size = 1.0;
   
      let mut stepper = Rk4::new(system, t0, y0, tf, step_size);
      let _ = stepper.integrate();
   
      let last_step = stepper.y_out().last().unwrap();
   
      assert!(
         last_step.x.is_nan() && last_step.y.is_nan()
      );
   }

   #[test]
   /// test writing a file
   fn write_file() {
      let data : &dyn BathymetryData = &ConstantDepth { h: 1000.0 };
      let system = WaveRayPath::new(data, 1.0);
      let y0 = State::new(0.0, 0.0, 1.0, -1.0);
   
      let t0 = 0.0;
      let tf = 10.0;
      let step_size = 1.0;
   
      assert!(output_to_file(system, t0, y0, tf, step_size).is_ok());
   }

   #[test]
   // test the k derivative function
   fn test_dk() {
      let k_mag = 150.0;
      let h = 3.5;
      let dhdx = 0.2;
      let dhdy = 0.2;

      let ans = dk_vector_dt(&k_mag, &h, &dhdx, &dhdy);

      assert!((ans.0 - -147.0).abs() < f64::EPSILON, "Expected -147, got {}", ans.0);
      assert!((ans.1 - -147.0).abs() < f64::EPSILON, "Expected -147, got {}", ans.1)

   }

}
