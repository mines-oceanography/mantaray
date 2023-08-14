//! This module makes it easier to use the Rk4 ray tracing by encapsulating it
//! with the SingleRay struct

use ode_solvers::{Rk4, OVector};

use crate::{bathymetry::BathymetryData, error::Error, WaveRayPath, State};
use crate::bathymetry::cartesian::{self, CartesianFile};
use std::path::Path;
use std::fs::File;
use std::io::Write;

// define x_out and y_out types
type XOut = Vec<f64>;
type YOut = Vec<OVector<f64, nalgebra::base::dimension::Const<4>>>;

// A struct with methods for tracing an individual wave and returning the result.
struct SingleRay<'a> {
    bathymetry_data: &'a dyn BathymetryData,
    initial_conditions: (f64, f64, f64, f64)
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
    pub fn new(bathymetry_data: &'a dyn BathymetryData, x0: f64, y0: f64, kx0: f64, ky0: f64) -> Self {
        SingleRay {
            bathymetry_data,
            initial_conditions: (x0, y0, kx0, ky0)
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
    pub fn trace_individual( &self, 
            start_time: f64, end_time: f64,
            step_size: f64 ) 
            -> Result<(XOut, YOut), Error>
        {

        // do the calculations
        let system = WaveRayPath::new(self.bathymetry_data);
        let s0 = State::new(
            self.initial_conditions.0, self.initial_conditions.1,
            self.initial_conditions.2, self.initial_conditions.3
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
fn output_to_tsv_file(file_name: &str, x_out: &XOut, y_out: &YOut) -> Result<File, Error> {
    
    let mut file = File::create(file_name).expect("could not open file");
    writeln!(&mut file, "t x y kx ky").expect("could not write to file");
    for (i, x) in x_out.iter().enumerate() {
       write!(&mut file, "{x} ").expect("could not write to file");
       for elem in y_out[i].iter() {
          write!(&mut file, "{elem} ").expect("could not write to file");
       }
       writeln!(&mut file, " ").expect("could not write to file");
    }
    Ok(file)
 
 }

#[cfg(test)]
 mod test_single_wave {

    use std::path::Path;
    use lockfile::Lockfile;

    use crate::bathymetry::{cartesian, BathymetryData};

    use super::{SingleRay, output_to_tsv_file};

        /// Create a constant depth file
        fn create_constant_depth_file(path: &Path, x_len: usize, y_len: usize, x_step: f32, y_step: f32) {
            
            let x_data: Vec<f32> = (0..x_len).map(|x| x as f32 * x_step).collect();
            let y_data: Vec<f32> = (0..y_len).map(|y| y as f32 * y_step).collect();
            let depth_data: Vec<f64> = (0..(x_len*y_len)).map(|_| 10.0f64).collect();

            // most below copied from the docs
            use netcdf3::{FileWriter, DataSet, Version};
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
                data_set.add_var_f64(depth_var_name, &[y_dim_name, x_var_name]).unwrap();
            
                data_set
            };
            
            // ...
            
            // Create and write the NetCDF-3 file
            // ----------------------------------
            let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
            // Set the NetCDF-3 definition
            file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
            assert_eq!(depth_var_len,                     x_var_len * y_var_len);
            file_writer.write_var_f32(y_var_name, &y_data[..]).unwrap();
            file_writer.write_var_f32(x_var_name, &x_data[..]).unwrap();
            file_writer.write_var_f64(depth_var_name, &depth_data[..]).unwrap();
             file_writer.close().unwrap();
            // end of copied from docs

        }

        /// Create a test file with depths split down the middle
        fn create_two_depth_file(path: &Path, x_len: usize, y_len: usize, x_step: f32, y_step: f32) {
            
            let x_data: Vec<f32> = (0..x_len).map(|x| x as f32 * x_step).collect();
            let y_data: Vec<f32> = (0..y_len).map(|y| y as f32 * y_step).collect();
            let depth_data: Vec<f64> = (0..(x_len*y_len))
                .map(|p| {
                    if p % x_len >= 50 {
                        20.0f64
                    } else {
                        50.0f64
                    }
                }
            )
                .collect();

            // most below copied from the docs
            use netcdf3::{FileWriter, DataSet, Version};
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
                data_set.add_var_f64(depth_var_name, &[y_dim_name, x_var_name]).unwrap();
            
                data_set
            };
            
            // ...
            
            // Create and write the NetCDF-3 file
            // ----------------------------------
            let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
            // Set the NetCDF-3 definition
            file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
            assert_eq!(depth_var_len,                     x_var_len * y_var_len);
            file_writer.write_var_f32(y_var_name, &y_data[..]).unwrap();
            file_writer.write_var_f32(x_var_name, &x_data[..]).unwrap();
            file_writer.write_var_f64(depth_var_name, &depth_data[..]).unwrap();
             file_writer.close().unwrap();
            // end of copied from docs

        }


        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_constant_wave_shallow_x() {
            let lockfile = Lockfile::create(Path::new("tmp_constant_depth_shallow_x.nc")).unwrap();
            create_constant_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
        
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        
            let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 0.01, 0.0);

            // make sure the starting point is at least 2 steps away from the edge.
            let res = wave.trace_individual(0.0, 8.0, 1.0).unwrap();
        
            let _ = output_to_tsv_file("constant_depth_shallow_x_out.txt",&res.0, &res.1);
        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_constant_wave_shallow_xy() {
            let lockfile = Lockfile::create(Path::new("tmp_constant_depth_shallow_xy.nc")).unwrap();
            create_constant_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
        
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());

            // test wave 2 starting in the corner
            let wave2 = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.007, 0.007);
            let res2 = wave2.trace_individual(0.0, 8.0, 1.0).unwrap();
            let _ = output_to_tsv_file("constant_depth_shallow_xy_out.txt",&res2.0, &res2.1);
        
        }
        
        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_constant_wave_deep_x() {
            let lockfile = Lockfile::create(Path::new("tmp_constant_depth_deep_x.nc")).unwrap();
            create_constant_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
        
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        
            // test wave 1
            let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 1.0, 0.0);

            // make sure the starting point is at least 2 steps away from the edge.
            let res = wave.trace_individual(0.0, 18.0, 1.0).unwrap();
        
            let _ = output_to_tsv_file("constant_depth_deep_x_out.txt",&res.0, &res.1);

        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_constant_wave_deep_xy() {
            let lockfile = Lockfile::create(Path::new("tmp_constant_depth_deep_xy.nc")).unwrap();
            create_constant_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
        
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        
            let wave2 = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.7, 0.7);
            let res2 = wave2.trace_individual(0.0, 18.0, 1.0).unwrap();
            let _ = output_to_tsv_file("constant_depth_deep_xy_out.txt",&res2.0, &res2.1);

        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_two_depth_wave_shallow_x() {
            let lockfile = Lockfile::create(Path::new("tmp_two_depth_shallow_x.nc")).unwrap();
            create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
                    
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        
            let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 0.01, 0.0);
            
            // make sure the starting point is at least 2 steps away from the edge.
            let res = wave.trace_individual(0.0, 6.0, 1.0).unwrap();
        
            let _ = output_to_tsv_file("two_depth_shallow_out_x.txt",&res.0, &res.1);

        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_two_depth_wave_shallow_xy() {
            let lockfile = Lockfile::create(Path::new("tmp_two_depth_shallow_xy.nc")).unwrap();
            create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
                    
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        

            let wave2 = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.007, 0.007);
            let res2 = wave2.trace_individual(0.0, 7.0, 1.0).unwrap();
            let _ = output_to_tsv_file("two_depth_shallow_xy_out.txt",&res2.0, &res2.1);

        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_two_depth_wave_deep_x() {
            let lockfile = Lockfile::create(Path::new("tmp_two_depth_deep_x.nc")).unwrap();
            create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
                    
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        
            let wave = SingleRay::new(bathymetry_data, 10.0, 50.0, 1.0, 0.0);
            
            // make sure the starting point is at least 2 steps away from the edge.
            let res = wave.trace_individual(0.0, 30.0, 1.0).unwrap();
        
            let _ = output_to_tsv_file("two_depth_deep_x_out.txt",&res.0, &res.1);
        
        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_two_depth_wave_deep_xy() {
            let lockfile = Lockfile::create(Path::new("tmp_two_depth_deep_xy.nc")).unwrap();
            create_two_depth_file(&lockfile.path(), 100, 100, 1.0, 1.0);
                    
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(&lockfile.path());
        
            let wave2 = SingleRay::new(bathymetry_data, 10.0, 10.0, 0.7, 0.7);
            let res2 = wave2.trace_individual(0.0, 40.0, 1.0).unwrap();
            let _ = output_to_tsv_file("two_depth_deep_xy_out.txt",&res2.0, &res2.1);

        }

 }