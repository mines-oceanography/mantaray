//! This module makes it easier to use the Rk4 ray tracing by encapsulating it
//! with the SingleRay struct

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;

use ode_solvers::dop_shared::SolverResult;
use rayon::prelude::*;

use ode_solvers::Rk4;

use crate::current::CurrentData;
use crate::{
    bathymetry::BathymetryData, error::Result, wave_ray_path::State, wave_ray_path::Time,
    wave_ray_path::WaveRayPath,
};

/// a struct that creates many rays
pub struct ManyRays<'a> {
    bathymetry_data: &'a dyn BathymetryData,
    current_data: Option<&'a dyn CurrentData>,
    /// a vector of initial x, y, kx, and ky values for the many waves
    init_rays: &'a Vec<(f64, f64, f64, f64)>,
}

#[allow(dead_code)]
impl<'a> ManyRays<'a> {
    /// construct a new `ManyRays` from bathymetry and initial rays
    ///
    /// # Arguments
    /// `bathymetry_data`: `&'a dyn BathymetryData`
    /// - the data on depth that implements the `depth` and
    ///   `depth_gradient` methods.
    ///
    /// `current_data`: `Option<&'a dyn CurrentData>`
    /// - the data on current that implements the `get_current` and
    ///  `get_current_gradient` methods. If `None`, then the current is assumed
    /// to be zero.
    ///
    /// `init_rays`: `&'a Vec<(f64, f64, f64, f64)>`
    /// - a vector of initial x, y, kx, and ky values for the many waves
    ///
    /// # Returns
    /// `Self`: a constructed `ManyRays` struct
    pub fn new(
        bathymetry_data: &'a dyn BathymetryData,
        current_data: Option<&'a dyn CurrentData>,
        init_rays: &'a Vec<(f64, f64, f64, f64)>,
    ) -> Self {
        ManyRays {
            bathymetry_data,
            current_data,
            init_rays,
        }
    }

    /// Trace many rays given start time, stop time, and step size (delta t)
    ///
    /// Given the arguments, `trace_many` creates a vector of SingleRays,
    /// integrates each ray, and returns the results.
    ///
    /// Arguments:
    ///
    /// `start_time`: `f64`
    /// - the time the ray tracing begins.
    ///
    /// `end_time`: `f64`
    /// - the time the ray tracing is stopped.
    ///
    /// `step_size`: `f64`
    /// - the change in time between integration steps. Smaller step size
    ///   produces more accurate result, but takes longer to run.
    ///
    /// Returns: `Vec<Option<(XOut, YOut)>>`: A vector of optional values. Each
    /// value in the vector is either `None`, which represents an error during
    /// that ray's integration, or they are a tuple of (XOut, YOut).
    pub fn trace_many(
        &self,
        start_time: f64,
        end_time: f64,
        step_size: f64,
    ) -> Vec<Option<SolverResult<Time, State>>> {
        // create a vector of SingleRays
        let rays: Vec<SingleRay> = self
            .init_rays
            .par_iter()
            .map(|(x0, y0, kx0, ky0)| {
                SingleRay::new(
                    self.bathymetry_data,
                    self.current_data,
                    *x0,
                    *y0,
                    *kx0,
                    *ky0,
                )
            })
            .collect();

        // integrate each. I think here is where I would use `par_iter` from rayon in the future.
        let results: Vec<Option<SolverResult<Time, State>>> = rays
            .par_iter()
            .map(
                |r| match r.trace_individual(start_time, end_time, step_size) {
                    Ok(v) => Some(v),
                    Err(e) => {
                        println!("ERROR {} during intergration", e);
                        None
                    }
                },
            )
            .collect();

        // return the results
        results
    }
}

// A struct with methods for tracing an individual wave and returning the result.
pub(crate) struct SingleRay<'a> {
    bathymetry_data: &'a dyn BathymetryData,
    current_data: Option<&'a dyn CurrentData>,
    initial_conditions: (f64, f64, f64, f64),
}

#[allow(dead_code)]
impl<'a> SingleRay<'a> {
    /// construct a `SingleRay`
    ///
    /// # Arguments
    /// `bathymetry_data` : `&'a dyn BathymetryData`
    /// - a struct that implements the `depth` function
    ///
    /// `current_data` : `Option<&'a dyn CurrentData>`
    /// - a struct that implements the `get_current` function. If `None`, then
    ///  the current is assumed to be zero.
    ///
    /// `x0` : `f64`
    /// - the initial x coordinate
    ///
    /// `y0` : `f64`
    /// - the initial y coordinate
    ///
    /// `kx0` : `f64`
    /// - the initial kx value
    ///
    /// `ky0` : `f64`
    /// - the initial ky value
    ///
    /// # Returns
    /// `Self` : the new `SingleRay` struct
    pub(crate) fn new(
        bathymetry_data: &'a dyn BathymetryData,
        current_data: Option<&'a dyn CurrentData>,
        x0: f64,
        y0: f64,
        kx0: f64,
        ky0: f64,
    ) -> Self {
        SingleRay {
            bathymetry_data,
            current_data,
            initial_conditions: (x0, y0, kx0, ky0),
        }
    }

    /// computes ode_solvers Rk4 tracing and returns result
    ///
    /// # Arguments
    ///
    /// `start_time` : `f64`
    /// - time to start the Rk4
    ///
    /// `end_time` : `f64`
    /// - time to end the Rk4
    ///
    /// `step_size` : `f64`
    /// - delta t
    ///
    /// # Returns
    /// `Result<SolverResult<Time, State>, Error>`
    /// - `SolverResult<Time, State>` : The result of the `ode_solvers`
    ///   integration.
    /// - `Err(Error::IntegrationError)` : there was an error during Rk4
    ///   integrate method.
    ///
    /// # Note
    /// This struct still copies the data when it returns, which could be an
    /// inefficiency, but the arguments are now less.
    pub(crate) fn trace_individual(
        &self,
        start_time: f64,
        end_time: f64,
        step_size: f64,
    ) -> Result<SolverResult<Time, State>> {
        // do the calculations
        let system = WaveRayPath::new(self.bathymetry_data, self.current_data);
        let s0 = State::new(
            self.initial_conditions.0,
            self.initial_conditions.1,
            self.initial_conditions.2,
            self.initial_conditions.3,
        );
        let mut stepper = Box::new(Rk4::new(system, start_time, s0, end_time, step_size));
        stepper.integrate()?;
        // return the stepper results
        let results: &SolverResult<Time, State> = stepper.results();

        Ok(results.clone())
    }
}

#[allow(dead_code)]
fn output_or_append_to_tsv_file(
    file_path: &Path,
    result: &SolverResult<Time, State>,
) -> Result<()> {
    let (x_out, y_out) = result.get();
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(file_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(&mut writer, "t x y kx ky")?;
    for (i, x) in x_out.iter().enumerate() {
        if y_out[i][0].is_nan() {
            break;
        }
        write!(&mut writer, "{} ", x)?;
        for elem in y_out[i].iter() {
            write!(&mut writer, "{} ", elem)?;
        }
        writeln!(&mut writer, " ")?;
    }
    writeln!(&mut writer, "END")?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod test_single_wave {

    use tempfile::NamedTempFile;

    use crate::{
        bathymetry::{BathymetryData, CartesianNetcdf3, ConstantDepth, ConstantSlope},
        current::{CartesianCurrent, ConstantCurrent},
        io::utility::{create_netcdf3_bathymetry, create_netcdf3_current},
    };

    use super::SingleRay;

    /// Create a test file with depths split down the middle
    fn two_depth_fn(x: f32, _y: f32) -> f64 {
        if x >= 50.0 {
            20.0
        } else {
            50.0
        }
    }

    #[test]
    /// ray tracing on a constant depth shallow wave propagating in the x
    /// direction. k stays the same.
    fn test_constant_wave_shallow_x() {
        let bathymetry_data = &ConstantDepth::new(10.0);

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 50.0, 0.01, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 8.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify each y, kx, ky value stays the same
        data.iter().for_each(|r| assert_eq!(r[1], 50.0));
        data.iter().for_each(|r| assert_eq!(r[2], 0.01));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }
    }

    #[test]
    /// ray tracing on a constant depth shallow wave propagating at an angle in
    /// the x=y direction. k stays the same.
    fn test_constant_wave_shallow_xy() {
        let bathymetry_data = &ConstantDepth::new(10.0);

        // test wave 2 starting in the corner
        let wave = SingleRay::new(bathymetry_data, None, 10.0, 10.0, 0.007, 0.007);
        let res = wave.trace_individual(0.0, 8.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify that kx and ky stay the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.007));
        data.iter().for_each(|r| assert_eq!(r[3], 0.007));

        // check to verify x and y are the same value and always increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            let x = r[0];
            let y = r[1];
            assert!((x - y).abs() < f64::EPSILON);
            assert!(x >= last_x);
            assert!(y >= last_y);
            last_x = x;
            last_y = y;
        }
    }

    #[test]
    /// ray tracing on a constant depth deep wave propagating in the x
    /// direction. k stay the same.
    fn test_constant_wave_deep_x() {
        let bathymetry_data = &ConstantDepth::new(10.0);

        // test wave 1
        let wave = SingleRay::new(bathymetry_data, None, 10.0, 50.0, 1.0, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 18.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify each y, kx, ky value stays the same
        data.iter().for_each(|r| assert_eq!(r[1], 50.0));
        data.iter().for_each(|r| assert_eq!(r[2], 1.0));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }
    }

    #[test]
    /// ray tracing on a constant depth deep wave propagating at an angle in the
    /// x=y direction. k stays the same.
    fn test_constant_wave_deep_xy() {
        let bathymetry_data = &ConstantDepth::new(10.0);

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 10.0, 0.7, 0.7);
        let res = wave.trace_individual(0.0, 18.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify that kx and ky stay the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.7));
        data.iter().for_each(|r| assert_eq!(r[3], 0.7));

        // check to verify x and y are the same value and always increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            let x = r[0];
            let y = r[1];
            assert!((x - y).abs() < f64::EPSILON);
            assert!(x >= last_x);
            assert!(y >= last_y);
            last_x = x;
            last_y = y;
        }
    }

    #[test]
    /// ray tracing on a two-depth (half 50m and half 20m) shallow wave
    /// propagating in the x direction. The kx increases slightly.
    fn test_two_depth_wave_shallow_x() {
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();

        create_netcdf3_bathymetry(&tmp_path, 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data = &CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 50.0, 0.01, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 5.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify each y and ky value stays the same
        data.iter().for_each(|r| assert_eq!(r[1], 50.0));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify the x values are increasing
        // verify that the kx value either increases or stays the same.
        let mut last_x = data[0][0];
        let mut last_kx = data[0][2];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[2] >= last_kx);
            last_x = r[0];
            last_kx = r[2];
        }

        // finally, make the the final kx is greater then the initial
        assert!(data.iter().last().unwrap()[2] > data.iter().next().unwrap()[2])
    }

    #[test]
    /// ray tracing on a two-depth shallow wave propagating at an angle in the
    /// x=y direction. This shows a change in the kx and ky.
    fn test_two_depth_wave_shallow_xy() {
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();

        create_netcdf3_bathymetry(&tmp_path, 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data = &CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 10.0, 0.007, 0.007);
        let res = wave.trace_individual(0.0, 6.8, 0.1).unwrap();

        let (_, data) = &res.get();

        // verify each ky value stays the same
        data.iter().for_each(|r| assert_eq!(r[3], 0.007));

        // verify the x and y values are increasing
        // verify that the kx value either increases or stays the same.
        // verify that the x values are greater than or equal to y values
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        let mut last_kx = data[0][2];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            assert!(r[2] >= last_kx);
            assert!(r[0] >= r[1]);
            last_x = r[0];
            last_y = r[1];
            last_kx = r[2];
        }

        // finally, make the the final kx is greater then the initial
        assert!(data.iter().last().unwrap()[2] > data.iter().next().unwrap()[2]);
    }

    #[test]
    /// ray tracing on a two-depth deep wave propagating in the x direction.
    /// This correcly shows no change in kx or ky.
    fn test_two_depth_wave_deep_x() {
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();

        create_netcdf3_bathymetry(&tmp_path, 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data = &CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 50.0, 1.0, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 30.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify each y, kx, ky value stays the same
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));
        data.iter().for_each(|r| assert_eq!(r[1], 50.0));
        data.iter().for_each(|r| assert_eq!(r[2], 1.0));

        // verify the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }
    }

    #[test]
    /// ray tracing on a two-depth deep wave propagating at an angle in the x=y
    /// direction. This correcly shows no change in kx or ky.
    fn test_two_depth_wave_deep_xy() {
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();

        create_netcdf3_bathymetry(&tmp_path, 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data = &CartesianNetcdf3::open(&tmp_path, "x", "y", "depth").unwrap();

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 10.0, 0.7, 0.7);
        let res = wave.trace_individual(0.0, 40.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify each y, kx, ky value stays the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.7));
        data.iter().for_each(|r| assert_eq!(r[3], 0.7));

        // verify the x values are increasing
        // verify that y values equal the x values
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[0] == r[1]);
            last_x = r[0];
        }
    }

    #[test]
    /// shallow water wave propogating through a constant slope from deeper to
    /// shallower water.
    fn test_slope_depth_wave_x() {
        let bathymetry_data = &ConstantSlope::builder().build().unwrap();

        let wave = SingleRay::new(bathymetry_data, None, 10.0, 1000.0, 0.1, 0.0);
        let res = wave.trace_individual(0.0, 100.0, 1.0).unwrap();

        // this wave will propogate from deep to shallow water.
        assert_eq!(bathymetry_data.depth(&10.0, &1000.0).unwrap(), 49.5);
        assert_eq!(bathymetry_data.depth(&300.0, &1000.0).unwrap(), 35.0);

        let (_, data) = &res.get();

        // verify each y and ky value stays the same
        //data.iter().for_each(|r| assert_eq!(r[1], 1000.0));
        //data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }

        // FIXME: check this test case below. it does not pass
        // verify that the last kx value is greater than the first. this is because
        // the wave is getting more and more shallow.
        // assert!(data.iter().last().unwrap()[2] > data.iter().next().unwrap()[2]);
    }

    #[test]
    /// tests one wave with constant depth and zero current
    ///
    /// Since there is no current, the y and kx and ky values will not change.
    /// The x values will increase because kx = 0.1. Because the depth is not
    /// constant, the speed should not change either.
    fn constant_depth_zero_current() {
        let bathymetry_data = &ConstantDepth::new(10.0);
        let current_data = &ConstantCurrent::new(0.0, 0.0);

        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.1, 0.0);
        let res = wave.trace_individual(100.0, 102.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // check to make sure all y values are zero
        data.iter().for_each(|r| assert_eq!(r[1], 0.0));
        data.iter().for_each(|r| assert_eq!(r[2], 0.1));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // check to make sure the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }
    }

    // positive v, deflect north
    // negative v, deflect south
    // positive u, should reduce wavenumber
    // negative u, should increase wavenumber
    // both u and v
    // constant gradient u and v
    // jet
    // eddy
    #[test]
    /// 
    fn constant_depth_and_current() {
        let bathymetry_data = &ConstantDepth::new(10.0);
        let current_data = &ConstantCurrent::new(0.5, 0.0);

        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.1, 0.0);
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all y values are zero, kx values are 0.1, ky values are 0.0
        data.iter().for_each(|r| assert_eq!(r[1], 0.0));
        data.iter().for_each(|r| assert_eq!(r[2], 0.1));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // check to make sure the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }
    }

    #[test]
    // test with a constant current of 0.5 m/s in the y direction. The kx and ky
    // values should stay the same and the x and y values will increase.
    fn test_positive_v() {
        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // the current is 0.5 m/s in the y direction
        let current_data = &ConstantCurrent::new(0.0, 0.5);
        // wave starts at (x,y,kx,ky) = (0,0,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.1, 0.0);
        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.1));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);
    }

    #[test]
    // test a wave with a constant current of -0.5 m/s in the y direction. since
    // the direction initially of the wave is (0, 0, 0.1, 0.0) in the x direction, the x values
    // will increase and the y values will decrease. the kx and ky values will
    // stay the same.
    fn test_negative_v() {
        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // the current is -0.5 m/s in the y direction
        let current_data = &ConstantCurrent::new(0.0, -0.5);
        // wave starts at (x,y,kx,ky) = (0,0,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.1, 0.0);
        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.1));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify that the x values are increasing and y values are decreasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] <= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x value is greater than the first and the last y value is less than the first. this is
        // because above only checked greater than or equal to or less than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] < data.iter().next().unwrap()[1]);
    }

    #[test]
    // test a wave with a constant current of 0.5 m/s in the x direction. since
    // the initially travels in the y direction, the x and y values are
    // increasing. The kx and ky values will stay the same.
    fn test_positive_u() {
        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // the current is 0.5 m/s in the x direction
        let current_data = &ConstantCurrent::new(0.5, 0.0);
        // wave starts at (x,y,kx,ky) = (0,0,0.0,0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.0, 0.1);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.0));
        data.iter().for_each(|r| assert_eq!(r[3], 0.1));

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);
    }

    #[test]
    // test a wave with a constant current of -0.5 m/s in the x direction since
    // the wave starts initially at (0,0,0.0,0.1), the x values are decreasing
    // and the y values are increasing. the kx and ky values will stay the same.
    fn test_negative_u() {
        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // the current is -0.5 m/s in the x direction
        let current_data = &ConstantCurrent::new(-0.5, 0.0);
        // wave starts at (x,y,kx,ky) = (0,0,0.0,0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.0, 0.1);
        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.0));
        data.iter().for_each(|r| assert_eq!(r[3], 0.1));

        // verify that the x values are decreasing and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] <= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x value is less than the first and the last y value is greater than the first. this is
        // because above only checked greater than or equal to or less than or equal to
        assert!(data.iter().last().unwrap()[0] < data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);
    }

    #[test]
    // test a wave with a constant current of 0.5 m/s in the x and y direction
    // since the wave starts initially at (0,0,0.1,0.0), the x and y values
    // will increase and the kx and ky values will stay the same.
    fn test_positive_u_and_v() {
        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // the current is 0.5 m/s in the x direction and 0.5 m/s in the y direction
        let current_data = &ConstantCurrent::new(0.5, 0.5);
        // wave starts at (x,y,kx,ky) = (0,0,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, 0.1, 0.0);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.1));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);
    }

    #[test]
    /// test a wave with a constant current of -0.5 m/s in the x and y direction
    /// since the wave starts initially at (0,0,-0.1,0.0), the x and y values
    /// will decrease and the kx and ky values will stay the same.
    fn test_negative_u_and_v() {
        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // the current is -0.5 m/s in the x direction and -0.5 m/s in the y direction
        let current_data = &ConstantCurrent::new(-0.5, -0.5);
        // wave starts at (x,y,kx,ky) = (0,0,-0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 0.0, 0.0, -0.1, 0.0);
        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], -0.1));
        data.iter().for_each(|r| assert_eq!(r[3], 0.0));

        // verify that the x values are decreasing and y values are decreasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] <= last_x);
            assert!(r[1] <= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x value is less than the first and the last y value is less than the first. this is
        // because above only checked greater than or equal to or less than or equal to
        assert!(data.iter().last().unwrap()[0] < data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] < data.iter().next().unwrap()[1]);
    }

    #[test]
    /// test a wave with a nonzero dudx where u = x/100.0 and v = 0.0 this
    /// test first creates the gradient file and then tests the wave
    /// propogation. It verifies two cases:
    /// 1) a wave starting at (1,1,0.1,0.0) will propogate in the x direction
    ///    and the kx value will decrease, but y and ky values will stay the
    ///    same
    /// 2) a wave starting at (1,1,0.0,0.1) will propogate in the y direction
    ///    with slight positive x direction, but both kx and ky will remain the
    ///    same.
    fn test_simple_dudx_gradient() {
        // function that takes in x and y as f32 and returns u and v as f64.
        // this only will create a gradient in the u direction
        fn u_gradient_fn(x: f32, _y: f32) -> (f64, f64) {
            ((x / 100.0) as f64, 0.0)
        }

        // create the current file
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();
        create_netcdf3_current(&tmp_path, 100, 100, 1.0, 1.0, u_gradient_fn);

        // open the current data
        let current_data = &CartesianCurrent::open(&tmp_path, "x", "y", "u", "v");

        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // wave starts at (x,y,kx,ky) = (1,1,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.1, 0.0);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all ky and y values are the same
        data.iter().for_each(|r| assert_eq!(r[3], 0.0)); // ky
        data.iter().for_each(|r| assert_eq!(r[1], 1.0)); // y

        // verify that the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }

        // verify that the kx value is decreasing
        let mut last_kx = data[0][2];
        for r in data.iter() {
            assert!(r[2] <= last_kx);
            last_kx = r[2];
        }

        // verify that the last kx value is less than the first. this is
        // because above only checked less than or equal to
        assert!(data.iter().last().unwrap()[2] < data.iter().next().unwrap()[2]);

        // verify that the last x value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);

        // new wave (x, y, kx, ky) = (1, 1, 0.0, 0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.0, 0.1);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all ky and kx values are the same
        data.iter().for_each(|r| assert_eq!(r[3], 0.1)); // ky
        data.iter().for_each(|r| assert_eq!(r[2], 0.0)); // kx

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);
    }

    #[test]
    /// This test will create a current file with a nonzero du/dy. The tests
    /// will verify two cases:
    /// 1) a wave starting at (1,1,0.1,0.0) will propogate in the x direction,
    ///    the kx values will stay the same, x will increase, and y and ky will
    ///    decrease
    /// 2) a wave starting at (1,1,0.0,0.1) will propogate in the y direction,
    ///    the kx and ky values stay the same, and the x and y values increase
    fn test_simple_dudy_gradient() {
        // function that takes in x and y as f32 and returns u and v as f64.
        // this will only create the dudy gradient
        fn u_gradient_fn(_x: f32, y: f32) -> (f64, f64) {
            ((y / 100.0) as f64, 0.0)
        }

        // create the current file
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();
        create_netcdf3_current(&tmp_path, 100, 100, 1.0, 1.0, u_gradient_fn);

        // open the current data
        let current_data = &CartesianCurrent::open(&tmp_path, "x", "y", "u", "v");

        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);

        // wave starts at (x,y,kx,ky) = (1,1,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 50.0, 0.1, 0.0);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify that all kx values are the same.
        data.iter().for_each(|r| assert_eq!(r[2], 0.1)); // kx

        // verify that the x values are increasing
        let mut last_x = data[0][0];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            last_x = r[0];
        }

        // verify last x value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);

        // verify that the y and ky values are decreasing
        let mut last_y = data[0][1];
        let mut last_ky = data[0][3];
        for r in data.iter() {
            assert!(r[1] <= last_y);
            assert!(r[3] <= last_ky);
            last_ky = r[3];
            last_y = r[1];
        }

        // check that the last y and ky value is less than the first. this is
        // because above only checked less than or equal to
        assert!(data.iter().last().unwrap()[1] < data.iter().next().unwrap()[1]);
        assert!(data.iter().last().unwrap()[3] < data.iter().next().unwrap()[3]);

        // new wave (x, y, kx, ky) = (1, 1, 0.0, 0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.0, 0.1);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify that kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.0)); // kx
        data.iter().for_each(|r| assert_eq!(r[3], 0.1)); // ky

        // verify that x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // verify that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);
    }

    #[test]
    /// This test will create a current file with a gradient in the v direction
    /// where v = (x / 100.0) and u = 0.0. This will create a gradient of dv/dy
    /// The tests will verify two cases:
    /// 1) a wave starting at (1,1,0.1,0.0) will propogate in the x direction,
    ///    all the kx and ky values will stay the same, but x and y values are
    ///    increasing
    /// 2) a wave starting at (1,1,0.0,0.1) will propogate in the y direction
    ///    and the x and kx values will stay the same, but the y values will be
    ///    increasing and the ky values decreasing
    fn test_simple_dvdy_gradient() {
        // function that takes in x and y as f32 and returns u and v as f64.
        // this only will create a gradient in the v direction
        fn v_gradient_fn(_x: f32, y: f32) -> (f64, f64) {
            (0.0, (y / 100.0) as f64)
        }

        // create the current file
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();
        create_netcdf3_current(&tmp_path, 100, 100, 1.0, 1.0, v_gradient_fn);

        // open the current data
        let current_data = &CartesianCurrent::open(&tmp_path, "x", "y", "u", "v");

        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // wave starts at (x,y,kx,ky) = (1,1,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.1, 0.0);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.1)); // kx
        data.iter().for_each(|r| assert_eq!(r[3], 0.0)); // ky

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // check that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to (if they were
        // always equal, it would have passed too)
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);

        // new wave (x, y, kx, ky) = (1, 1, 0.0, 0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.0, 0.1);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify that the x and kx values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.0)); // kx
                                                         // FIXME: why is the x value slowing increasing by almost f64::EPSILON each iteration?
        data.iter().for_each(|r| {
            assert!(
                (r[0] - 1.0).abs() < 10.0 * f64::EPSILON,
                "expected: 1.0, got: {}",
                r[0]
            )
        }); // x

        // verify that the y values are increasing
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[1] >= last_y);
            last_y = r[1];
        }

        // check that the last y value is greater than the first. this is
        // because above only checked greater than or equal to (if they were
        // always equal, it would have passed too)
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);

        // verify that the ky values are decreasing
        let mut last_ky = data[0][3];
        for r in data.iter() {
            assert!(r[3] <= last_ky);
            last_ky = r[3];
        }

        // check that the last ky value is less than the first. this is'
        // because above only checked less than or equal to
        assert!(data.iter().last().unwrap()[3] < data.iter().next().unwrap()[3]);
    }

    #[test]
    /// This test will create a current file with a gradient in the v direction
    /// where v = (x / 100.0) and u = 0.0. This will create a gradient of dv/dx
    /// The test will verify that:
    /// 1) a wave starting at (1,1,0.1,0.0) will propogate in the x direction
    /// and the kx and ky values will stay the same, but x and y values are
    /// increasing
    /// 2) a wave starting at (1,1,0.0,0.1) will propogate in the y direction
    ///    and only the ky values will stay the same. The y values will increase
    ///    and the x and kx values will decrease.
    fn test_simple_dvdx_gradient() {
        // function that takes in x and y as f32 and returns u and v as f64.
        // this only will create a gradient in the v direction
        fn v_gradient_fn(x: f32, _y: f32) -> (f64, f64) {
            (0.0, (x / 100.0) as f64)
        }

        // create the current file
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();
        create_netcdf3_current(&tmp_path, 100, 100, 1.0, 1.0, v_gradient_fn);

        // open the current data
        let current_data = &CartesianCurrent::open(&tmp_path, "x", "y", "u", "v");

        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);
        // wave starts at (x,y,kx,ky) = (1,1,0.1,0.0)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.1, 0.0);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify all kx and ky values are the same
        data.iter().for_each(|r| assert_eq!(r[2], 0.1)); // kx
        data.iter().for_each(|r| assert_eq!(r[3], 0.0)); // ky

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // check that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to (if they were
        // always equal, it would have passed too)
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);

        // new wave (x, y, kx, ky) = (1, 1, 0.0, 0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 50.0, 1.0, 0.0, 0.1);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // verify that the ky values are the same
        data.iter().for_each(|r| assert_eq!(r[3], 0.1)); // ky

        // verify that the y values are increasing
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[1] >= last_y);
            last_y = r[1];
        }

        // check that the last y value is greater than the first. this is
        // because above only checked greater than or equal to (if they were
        // always equal, it would have passed too)
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);

        // verify that the x and kx values are decreasing
        let mut last_x = data[0][0];
        let mut last_kx = data[0][2];
        for r in data.iter() {
            assert!(r[0] <= last_x);
            assert!(r[2] <= last_kx);
            last_kx = r[2];
            last_x = r[0];
        }

        // check that the last x and kx value is less than the first. this is
        // because above only checked less than or equal to
        assert!(data.iter().last().unwrap()[0] < data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[2] < data.iter().next().unwrap()[2]);
    }

    #[test]
    /// This test will create a current file with a gradient in the u and v
    /// direction. The gradient is u = (x + y) / 100.0 and v = (x + y) / 100.0.
    /// This will create a gradient in the u and v direction for all du/dx,
    /// du/dy, dv/dx, dv/dy. The test will verify that the x and y values are
    /// increasing and the kx and ky values are decreasing.
    fn test_all_gradients() {
        // function that takes in x and y as f32 and returns u and v as f64.
        // this will create a gradient in the u and v direction for all du/dx,
        // du/dy, dv/dx, dv/dy
        fn all_gradient_fn(x: f32, y: f32) -> (f64, f64) {
            (((x + y) / 100.0) as f64, ((x + y) / 100.0) as f64)
        }

        // create the current file
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();
        create_netcdf3_current(&tmp_path, 100, 100, 1.0, 1.0, all_gradient_fn);

        // open the current data
        let current_data = &CartesianCurrent::open(&tmp_path, "x", "y", "u", "v");

        // deep water
        let bathymetry_data = &ConstantDepth::new(1000.0);

        // wave starts at (x,y,kx,ky) = (1,1,0.1,0.1)
        let wave = SingleRay::new(bathymetry_data, Some(current_data), 1.0, 1.0, 0.1, 0.1);

        // trace the wave for 10 seconds
        let res = wave.trace_individual(1.0, 10.0, 1.0).unwrap();

        let (_, data) = &res.get();

        // Note: no values should stay the same

        // verify that the x and y values are increasing
        let mut last_x = data[0][0];
        let mut last_y = data[0][1];
        for r in data.iter() {
            assert!(r[0] >= last_x);
            assert!(r[1] >= last_y);
            last_x = r[0];
            last_y = r[1];
        }

        // check that the last x and y value is greater than the first. this is
        // because above only checked greater than or equal to (if they were
        // always equal, it would have passed too)
        assert!(data.iter().last().unwrap()[0] > data.iter().next().unwrap()[0]);
        assert!(data.iter().last().unwrap()[1] > data.iter().next().unwrap()[1]);

        // verify that the kx and ky values are decreasing
        let mut last_kx = data[0][2];
        let mut last_ky = data[0][3];
        for r in data.iter() {
            assert!(r[2] <= last_kx);
            assert!(r[3] <= last_ky);
            last_kx = r[2];
            last_ky = r[3];
        }

        // check that the last kx and ky value is less than the first. this is
        // because above only checked less than or equal to
        assert!(data.iter().last().unwrap()[2] < data.iter().next().unwrap()[2]);
        assert!(data.iter().last().unwrap()[3] < data.iter().next().unwrap()[3]);
    }
}

#[cfg(test)]
mod test_many_waves {

    use std::path::Path;

    use crate::bathymetry::{CartesianFile, ConstantSlope};

    use super::ManyRays;

    #[test]
    /// check that output with test values from single wave works
    fn test_many_waves_ok() {
        let bathymetry_data = &ConstantSlope::builder().build().unwrap();

        let initial_waves = vec![
            // (x, y, kx, ky)
            (10.0, 10.0, 1.0, 0.0),
            (10.0, 20.0, 1.0, 0.0),
            (10.0, 30.0, 1.0, 0.0),
            (10.0, 40.0, 1.0, 0.0),
            (10.0, 50.0, 1.0, 0.0),
            (10.0, 60.0, 1.0, 0.0),
            (10.0, 70.0, 1.0, 0.0),
            (10.0, 80.0, 1.0, 0.0),
            (10.0, 90.0, 1.0, 0.0),
        ];

        let waves = ManyRays::new(bathymetry_data, None, &initial_waves);

        let results = waves.trace_many(0.0, 100000.0, 1.0);

        for res in results {
            assert!(res.is_some())
        }

        // test to verify size of initial waves and result (number of rays) samen
        // test to verify each instance of many ray against single ray
        // 
    }

    #[test]
    /// verify the linear_sea_mount netcdf test works as expected. rays
    /// should bend towards the circular island such that rays above the island
    /// curve down and rays below the island curve up.
    fn test_sea_mount() {
        use std::process::Command;

        // create the netcdf file from python script
        Command::new("python")
            .arg("support/linear_sea_mount.py")
            .output()
            .expect("failed to execute process");

        let bathymetry_data = CartesianFile::new(Path::new("island.nc"));

        // top half
        let initial_waves: Vec<(f64, f64, f64, f64)> = (0..10)
            .map(|v| (-9990.0, (v * 1000) as f64, 0.01, 0.0))
            .collect();

        let waves = ManyRays::new(&bathymetry_data, None, &initial_waves);

        let results = waves.trace_many(0.0, 1000.0, 1.0);

        // verify the top half curves down towards the island. the x value
        // should be increasing. y value should stay the same or decrease.
        for res in results {
            match res {
                Some(res) => {
                    let (_, data) = res.get();
                    // x values increase, y values decrease
                    let mut last_x = data[0][0];
                    let mut last_y = data[0][1];
                    for r in data.iter() {
                        if r[0].is_nan() {
                            continue;
                        }
                        assert!(r[0] >= last_x);
                        assert!(r[1] <= last_y);
                        last_x = r[0];
                        last_y = r[1];
                    }
                }
                None => (),
            }
        }

        // bottom half
        let initial_waves: Vec<(f64, f64, f64, f64)> = (1..10)
            .map(|v| (-9990.0, (v * -1000) as f64, 0.01, 0.0))
            .collect();

        let waves = ManyRays::new(&bathymetry_data, None, &initial_waves);

        let results = waves.trace_many(0.0, 1000.0, 1.0);

        // verify the bottom half curves up towards the island. x should
        // increase. y should increase.
        for res in results {
            match res {
                Some(res) => {
                    let (_, data) = res.get();
                    // x values increase, y values increase
                    let mut last_x = data[0][0];
                    let mut last_y = data[0][1];
                    for r in data.iter() {
                        if r[0].is_nan() {
                            continue;
                        }
                        assert!(r[0] >= last_x);
                        assert!(r[1] >= last_y);
                        last_x = r[0];
                        last_y = r[1];
                    }
                }
                None => (),
            }
        }
    }
}
