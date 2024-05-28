//! This module makes it easier to use the Rk4 ray tracing by encapsulating it
//! with the SingleRay struct

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;

use ode_solvers::dop_shared::SolverResult;
use rayon::prelude::*;

use ode_solvers::Rk4;

use crate::{
    bathymetry::BathymetryData, error::Error, wave_ray_path::State, wave_ray_path::Time,
    wave_ray_path::WaveRayPath,
};

/// a struct that creates many rays
pub struct ManyRays<'a> {
    bathymetry_data: &'a dyn BathymetryData,
    /// a vector of initial x, y, kx, and ky values for the many waves
    init_rays: &'a Vec<(f64, f64, f64, f64)>,
}

#[allow(dead_code)]
impl<'a> ManyRays<'a> {
    /// construct a new `ManyRays` from bathymetry and initial rays
    ///
    /// # Arguments
    /// `bathymetry_data`: `&'a dyn BathymetryData`
    /// - the data on depth that implements the `get_depth` and
    ///   `get_depth_gradient` methods.
    ///
    /// `init_rays`: `&'a Vec<(f64, f64, f64, f64)>`
    /// - a vector of initial x, y, kx, and ky values for the many waves
    ///
    /// # Returns
    /// `Self`: a constructed `ManyRays` struct
    pub fn new(
        bathymetry_data: &'a dyn BathymetryData,
        init_rays: &'a Vec<(f64, f64, f64, f64)>,
    ) -> Self {
        ManyRays {
            bathymetry_data,
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
            .map(|(x0, y0, kx0, ky0)| SingleRay::new(self.bathymetry_data, *x0, *y0, *kx0, *ky0))
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
struct SingleRay<'a> {
    bathymetry_data: &'a dyn BathymetryData,
    initial_conditions: (f64, f64, f64, f64),
}

#[allow(dead_code)]
impl<'a> SingleRay<'a> {
    /// construct a `SingleRay`
    ///
    /// # Arguments
    /// `bathymetry_data` : `&'a dyn BathymetryData`
    /// - a struct that implements the `get_depth` function
    ///
    /// # Arguments
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
    pub fn new(
        bathymetry_data: &'a dyn BathymetryData,
        x0: f64,
        y0: f64,
        kx0: f64,
        ky0: f64,
    ) -> Self {
        SingleRay {
            bathymetry_data,
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
    pub fn trace_individual(
        &self,
        start_time: f64,
        end_time: f64,
        step_size: f64,
    ) -> Result<SolverResult<Time, State>, Error> {
        // do the calculations
        let system = WaveRayPath::new(self.bathymetry_data, None);
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
) -> Result<(), Error> {
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

    use lockfile::Lockfile;
    use std::path::Path;
    use tempfile::tempdir;

    use crate::{
        {bathymetry::{BathymetryData, CartesianFile, ConstantDepth, ConstantSlope},
        ray_result::RayResults,
    }, io::utility::create_netcdf3_bathymetry};

    use super::SingleRay;

    /// Create a test file with depths split down the middle
    fn two_depth_fn(x: f32, _y: f32) -> f64 {
        if x >= 50.0 {
            20.0
        } else {
            50.0
        }
    }

    fn temp_filename(filename: &str) -> String {
        let tmp_dir = tempdir().unwrap();
        tmp_dir
            .path()
            .join(filename)
            .into_os_string()
            .into_string()
            .unwrap()
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a constant depth
    /// shallow wave propagating in the x direction.
    fn test_constant_wave_shallow_x() {
        let bathymetry_data: &dyn BathymetryData = &ConstantDepth::new(10.0);
        let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 0.01, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 8.0, 1.0).unwrap();

        let filename = temp_filename("constant_depth_shallow_x_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a constant depth
    /// shallow wave propagating at an angle in the x=y direction.
    fn test_constant_wave_shallow_xy() {
        let bathymetry_data: &dyn BathymetryData = &ConstantDepth::new(10.0);

        // test wave 2 starting in the corner
        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.007, 0.007);
        let res = wave.trace_individual(0.0, 8.0, 1.0).unwrap();
        let filename = temp_filename("constant_depth_shallow_x_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a constant depth
    /// deep wave propagating in the x direction.
    fn test_constant_wave_deep_x() {
        let bathymetry_data: &dyn BathymetryData = &ConstantDepth::new(10.0);

        // test wave 1
        let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 1.0, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 18.0, 1.0).unwrap();

        let filename = temp_filename("constant_depth_deep_x_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a constant depth
    /// deep wave propagating at an angle in the x=y direction.
    fn test_constant_wave_deep_xy() {
        let bathymetry_data: &dyn BathymetryData = &ConstantDepth::new(10.0);

        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.7, 0.7);
        let res = wave.trace_individual(0.0, 18.0, 1.0).unwrap();

        let filename = temp_filename("constant_depth_deep_xy_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth shallow
    /// wave propagating in the x direction. The kx increases slightly.
    fn test_two_depth_wave_shallow_x() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_shallow_x.nc")).unwrap();
        create_netcdf3_bathymetry(&lockfile.path(), 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 0.01, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 6.0, 1.0).unwrap();

        let filename = temp_filename("two_depth_shallow_x_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth shallow
    /// wave propagating at an angle in the x=y direction. NOTE: as this
    /// function is written, the output kx and ky do not change. However,
    /// decreasing the step size will make the kx and ky change.
    fn test_two_depth_wave_shallow_xy() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_shallow_xy.nc")).unwrap();
        create_netcdf3_bathymetry(&lockfile.path(), 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.007, 0.007);
        let res = wave.trace_individual(0.0, 7.0, 1.0).unwrap();

        let filename = temp_filename("two_depth_shallow_xy_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth deep wave
    /// propagating in the x direction. This correcly shows no change in kx or ky.
    fn test_two_depth_wave_deep_x() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_deep_x.nc")).unwrap();
        create_netcdf3_bathymetry(&lockfile.path(), 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 1.0, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 30.0, 1.0).unwrap();

        let filename = temp_filename("two_depth_deep_x_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth deep wave
    /// propagating at an angle in the x=y direction. This correcly shows no change in kx or ky.
    fn test_two_depth_wave_deep_xy() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_deep_xy.nc")).unwrap();
        create_netcdf3_bathymetry(&lockfile.path(), 100, 100, 1.0, 1.0, two_depth_fn);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.7, 0.7);
        let res = wave.trace_individual(0.0, 40.0, 1.0).unwrap();

        let filename = temp_filename("two_depth_deep_xy_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }

    #[test]
    /// shallow water
    fn test_slope_depth_wave_x() {
        let bathymetry_data: &dyn BathymetryData = &ConstantSlope::builder().build().unwrap();

        let wave = SingleRay::new(bathymetry_data, 10.0, 1000.0, 0.01, 0.0);
        let res = wave.trace_individual(0.0, 100.0, 1.0).unwrap();

        let filename = temp_filename("slope_depth_x_out.txt");
        let _ = RayResults::from(res).save_file(Path::new(&filename));
    }
}

#[cfg(test)]
mod test_many_waves {

    use crate::bathymetry::{BathymetryData, ConstantSlope};

    use super::ManyRays;

    #[test]
    /// check that output with test values from single wave works
    fn test_many_waves_ok() {
        let bathymetry_data: &dyn BathymetryData = &ConstantSlope::builder().build().unwrap();

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

        let waves = ManyRays::new(bathymetry_data, &initial_waves);

        let results = waves.trace_many(0.0, 55.0, 1.0);

        for res in results {
            assert!(res.is_some())
        }
    }
}
