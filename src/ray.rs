//! This module makes it easier to use the Rk4 ray tracing by encapsulating it
//! with the SingleRay struct

use rayon::prelude::*;

use ode_solvers::{OVector, Rk4};

use crate::{BathymetryData, error::Error, State, WaveRayPath};
use std::fs::{File, OpenOptions};
use std::io::{Write, BufWriter};
use std::path::Path;

// define x_out and y_out types
type XOut = Vec<f64>;
type YOut = Vec<OVector<f64, nalgebra::base::dimension::Const<4>>>;

/// a struct that creates many rays
struct ManyRays<'a> {
    bathymetry_data: &'a dyn BathymetryData,
    /// a vector of initial x, y, kx, and ky values for the many waves
    init_rays: &'a Vec<(f64, f64, f64, f64)>
}

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
    fn new(bathymetry_data: &'a dyn BathymetryData, init_rays: &'a Vec<(f64, f64, f64, f64)>) -> Self {
        ManyRays { bathymetry_data, init_rays }
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
        fn trace_many(&self, start_time: f64, end_time: f64, step_size: f64) -> Vec<Option<(XOut, YOut)>> {

        // create a vector of SingleRays
        let rays: Vec<SingleRay> = self.init_rays
            .par_iter()
            .map(|(x0, y0, kx0, ky0)| {
                SingleRay::new(self.bathymetry_data, *x0, *y0, *kx0, *ky0)
            })
            .collect();

        // integrate each. I think here is where I would use `par_iter` from rayon in the future.
        let results: Vec<Option<(XOut, YOut)>> = rays
            .par_iter()
            .map(|r| {
                match r.trace_individual(start_time, end_time, step_size) {
                    Ok(v) => Some(v),
                    Err(e) => {
                        println!("ERROR {} during intergration", e);
                        None      
                    }
                }
            })
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
    /// `Result<(XOut, YOut), Error>`
    /// - `Ok((XOut, YOut))` : a tuple of the x_out and y_out from the Rk4
    ///   stepper.
    /// - `Err(Error::IntegrationError)` : there was an error during Rk4
    ///   integrate method.
    ///
    /// # Panics
    /// This function will panic if the wave ever propagates out of bounds of
    /// the bathymetry data array during Rk4 integration.
    ///
    /// # Note
    /// This struct still copies the data when it returns, which could be an
    /// inefficiency, but the arguments are now less.
    pub fn trace_individual(
        &self,
        start_time: f64,
        end_time: f64,
        step_size: f64,
    ) -> Result<(XOut, YOut), Error> {
        // do the calculations
        let system = WaveRayPath::new(self.bathymetry_data);
        let s0 = State::new(
            self.initial_conditions.0,
            self.initial_conditions.1,
            self.initial_conditions.2,
            self.initial_conditions.3,
        );
        let mut stepper = Box::new(Rk4::new(system, start_time, s0, end_time, step_size));
        stepper.integrate()?;
        // return the stepper results
        let x_out: &XOut = stepper.x_out();
        let y_out: &YOut = stepper.y_out();

        // FIXME: how to prevent copying this data? This will take longer,
        // since it has to copy the data over, but I don't see another way
        // since stepper is a local variable, and when it goes out of scope,
        // it's references will too. The solution is I need to figure out
        // the type of stepper and have it be stored inside Single as a
        // member variable.
        Ok((x_out.clone(), y_out.clone()))
    }

    /// set the initial conditions
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
    pub fn set_initial_conditions(&mut self, x0: f64, y0: f64, kx0: f64, ky0: f64) {
        self.initial_conditions = (x0, y0, kx0, ky0)
    }
}

// output to space separated file
fn output_to_tsv_file(file_name: &str, x_out: &XOut, y_out: &YOut) -> Result<(), Error> {
    let file = File::create(file_name)?;
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
    writer.flush()?;
    Ok(())
}

fn output_or_append_to_tsv_file(
    file_name: &str,
    x_out: &XOut,
    y_out: &YOut,
) -> Result<(), Error> {
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(file_name)?;
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

    use crate::{BathymetryData, bathymetry::{CartesianFile, ConstantSlope, ConstantDepth}};

    use super::{output_to_tsv_file, SingleRay};

    /// Create a test file with depths split down the middle
    fn create_two_depth_file(path: &Path, x_len: usize, y_len: usize, x_step: f32, y_step: f32) {
        let x_data: Vec<f32> = (0..x_len).map(|x| x as f32 * x_step).collect();
        let y_data: Vec<f32> = (0..y_len).map(|y| y as f32 * y_step).collect();
        let depth_data: Vec<f64> = (0..(x_len * y_len))
            .map(|p| if p % x_len >= 50 { 20.0f64 } else { 50.0f64 })
            .collect();

        // most below copied from the docs
        use netcdf3::{DataSet, FileWriter, Version};
        let y_dim_name: &str = "y";
        let y_var_name: &str = y_dim_name;
        let y_var_len: usize = y_len;

        let x_dim_name: &str = "x";
        let x_var_name: &str = x_dim_name;
        let x_var_len: usize = x_len;

        let depth_var_name: &str = "depth";
        let depth_var_len: usize = depth_data.len();

        // Create the NetCDF-3 definition
        // ------------------------------
        let data_set: DataSet = {
            let mut data_set: DataSet = DataSet::new();
            // Define the dimensions
            data_set.add_fixed_dim(y_dim_name, y_var_len).unwrap();
            data_set.add_fixed_dim(x_dim_name, x_var_len).unwrap();
            // Define the variable
            data_set.add_var_f32(y_var_name, &[y_dim_name]).unwrap();
            data_set.add_var_f32(x_var_name, &[x_var_name]).unwrap();
            data_set
                .add_var_f64(depth_var_name, &[y_dim_name, x_var_name])
                .unwrap();

            data_set
        };

        // ...

        // Create and write the NetCDF-3 file
        // ----------------------------------
        let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
        // Set the NetCDF-3 definition
        file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
        assert_eq!(depth_var_len, x_var_len * y_var_len);
        file_writer.write_var_f32(y_var_name, &y_data[..]).unwrap();
        file_writer.write_var_f32(x_var_name, &x_data[..]).unwrap();
        file_writer
            .write_var_f64(depth_var_name, &depth_data[..])
            .unwrap();
        file_writer.close().unwrap();
        // end of copied from docs
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

        let _ = output_to_tsv_file("constant_depth_shallow_x_out.txt", &res.0, &res.1);
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
        let _ = output_to_tsv_file("constant_depth_shallow_xy_out.txt", &res.0, &res.1);
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

        let _ = output_to_tsv_file("constant_depth_deep_x_out.txt", &res.0, &res.1);
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a constant depth
    /// deep wave propagating at an angle in the x=y direction.
    fn test_constant_wave_deep_xy() {
        let bathymetry_data: &dyn BathymetryData = &ConstantDepth::new(10.0);

        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.7, 0.7);
        let res = wave.trace_individual(0.0, 18.0, 1.0).unwrap();
        let _ = output_to_tsv_file("constant_depth_deep_xy_out.txt", &res.0, &res.1);
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth shallow
    /// wave propagating in the x direction. The kx increases slightly.
    fn test_two_depth_wave_shallow_x() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_shallow_x.nc")).unwrap();
        create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 0.01, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 6.0, 1.0).unwrap();

        let _ = output_to_tsv_file("two_depth_shallow_out_x.txt", &res.0, &res.1);
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth shallow
    /// wave propagating at an angle in the x=y direction. NOTE: as this
    /// function is written, the output kx and ky do not change. However,
    /// decreasing the step size will make the kx and ky change.
    fn test_two_depth_wave_shallow_xy() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_shallow_xy.nc")).unwrap();
        create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.007, 0.007);
        let res = wave.trace_individual(0.0, 7.0, 1.0).unwrap();
        let _ = output_to_tsv_file("two_depth_shallow_xy_out.txt", &res.0, &res.1);
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth deep wave
    /// propagating in the x direction. This correcly shows no change in kx or ky.
    fn test_two_depth_wave_deep_x() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_deep_x.nc")).unwrap();
        create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 1.0, 0.0);

        // make sure the starting point is at least 2 steps away from the edge.
        let res = wave.trace_individual(0.0, 30.0, 1.0).unwrap();

        let _ = output_to_tsv_file("two_depth_deep_x_out.txt", &res.0, &res.1);
    }

    #[test]
    // this test does not check anything yet, but outputs the result to a space separated file
    /// generate an output file showing ray tracing on a two-depth deep wave
    /// propagating at an angle in the x=y direction. This correcly shows no change in kx or ky.
    fn test_two_depth_wave_deep_xy() {
        let lockfile = Lockfile::create(Path::new("tmp_two_depth_deep_xy.nc")).unwrap();
        create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);

        let bathymetry_data: &dyn BathymetryData = &CartesianFile::new(&lockfile.path());

        let wave = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.7, 0.7);
        let res = wave.trace_individual(0.0, 40.0, 1.0).unwrap();
        let _ = output_to_tsv_file("two_depth_deep_xy_out.txt", &res.0, &res.1);
    }

    #[test]
    /// shallow water
    fn test_slope_depth_wave_x() {
        let bathymetry_data: &dyn BathymetryData = &ConstantSlope::builder().build().unwrap();

        let wave = SingleRay::new(bathymetry_data, 10.0, 1000.0, 0.01, 0.0);
        let res = wave.trace_individual(0.0, 100.0, 1.0).unwrap();
        let _ = output_to_tsv_file("slope_depth_x_out.txt", &res.0, &res.1);
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

        let initial_waves = vec![ // (x, y, kx, ky)
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
