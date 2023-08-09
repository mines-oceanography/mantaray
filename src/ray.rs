//! This module makes it easier to use the Rk4 ray tracing by encapsulating it with the SingleRay struct

use ode_solvers::{Rk4, OVector};

use crate::{bathymetry::BathymetryData, error::Error, WaveRayPath, State};
use crate::bathymetry::cartesian::{self, CartesianFile};
use std::path::Path;
use std::fs::File;
use std::io::Write;

// A struct with methods for tracing an individual wave and returning the result.
struct SingleRay;

impl SingleRay {
    /// construct a `SingleRay`.
    pub fn new() -> Self {
        SingleRay {}
    }

    /// computes ode_solvers Rk4 tracing and returns result
    /// 
    /// # Arguments
    /// 
    /// `bathymetry_data`: `&'b dyn BathymetryData`
    /// - the data type that implements `get_depth` representing the depth of
    ///   the area.
    /// 
    /// `x0` : `&f64`
    /// - the initial x position
    /// 
    /// `y0` : `&f64`
    /// - the initial y position
    /// 
    /// `kx0` : `&f64`
    /// - the initial kx value
    /// 
    /// `ky0` : `&f64`
    /// - the initial ky value
    /// 
    /// `start_time` : `&f64`
    /// - time to start the Rk4
    /// 
    /// `end_time` :
    /// - time to end the Rk4
    /// 
    /// `step_size` : `&f64`
    /// - delta t
    /// 
    /// # Returns
    /// `Result<(Vec<f64>, Vec<OVector<f64,
    /// nalgebra::base::dimension::Const<4>>>), Error>`
    /// - `Ok((Vec<f64>, Vec<OVector<f64,
    ///   nalgebra::base::dimension::Const<4>>>))` : a tuple of the x_out and
    ///   y_out from the Rk4 stepper.
    /// - `Err(Error::IntegrationError)` : there was an error during Rk4
    ///   integrate method.
    /// 
    /// # Panics
    /// This function will panic if the wave ever propagates out of bounds of
    /// the bathymetry data array during Rk4 integration.
    /// 
    /// # Note
    /// This takes too many arguments, is confusing, and also unoptimized. I
    /// need to move bathymetry_data to being a member variable of the struct.
    pub fn trace_individual<'b>( &self, 
            bathymetry_data: &'b dyn BathymetryData, 
            x0: &f64, y0: &f64, kx0: &f64, ky0: &f64, 
            start_time: &f64, end_time: &f64,
            step_size: &f64 ) 
            -> Result<(Vec<f64>, Vec<OVector<f64, nalgebra::base::dimension::Const<4>>>), Error>
        {

        // do the calculations
        let system = WaveRayPath::new(bathymetry_data, *step_size);
        let y0 = State::new(*x0, *y0, *kx0, *ky0);
        let mut stepper = Box::new(Rk4::new(system, *start_time, y0, *end_time, *step_size));
        if stepper.integrate().is_ok() {

            // return the stepper results
            let x_out: &Vec<f64> = stepper.x_out();
            let y_out: &Vec<OVector<f64, nalgebra::base::dimension::Const<4>>> = stepper.y_out();

            // FIXME: how to prevent copying this data? This will take longer,
            // since it has to copy the data over, but I don't see another way
            // since stepper is a local variable, and when it goes out of scope,
            // it's references will too. The solution is I need to figure out
            // the type of stepper and have it be stored inside Single as a
            // member variable.
            return Ok((x_out.clone(), y_out.clone()));

        } else {
            return Err(Error::IntegrationError);
        }

    }

}

// output to space separated file
fn output_to_tsv_file(file_name: &str, x_out: &Vec<f64>, y_out: &Vec<OVector<f64, nalgebra::base::dimension::Const<4>>>) -> Result<File, Error> {
    
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


 mod test_single_wave {

    use std::path::Path;

    use crate::bathymetry::{cartesian, BathymetryData};

    use super::{SingleRay, output_to_tsv_file};

        /// Create a constant depth file
        fn create_constant_depth_file(path: &Path, x_len: usize, y_len: usize, x_step: f32, y_step: f32) {
            
            let x_data: Vec<f32> = (0..x_len).map(|x| x as f32 * x_step).collect();
            let y_data: Vec<f32> = (0..y_len).map(|y| y as f32 * y_step).collect();
            let depth_data: Vec<f64> = (0..(x_len*y_len)).map(|_| 50.0f64).collect();

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
                        48.0f64
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
        fn test_constant_wave() {
            let wave = SingleRay::new();
        
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(Path::new("constant_depth.nc"));
        
            // make sure the starting point is at least 2 steps away from the edge.
            let res = wave.trace_individual(bathymetry_data, &10.0, &50.0, &0.12, &0.0, &0.0, &18.0, &1.0).unwrap();
        
            let _ = output_to_tsv_file("y_out.txt",&res.0, &res.1);
        }

        #[test]
        // this test does not check anything yet, but outputs the result to a space separated file
        fn test_two_depth_wave() {
            let wave = SingleRay::new();
        
            let bathymetry_data: &dyn BathymetryData = &cartesian::CartesianFile::new(Path::new("two_depth.nc"));
        
            // make sure the starting point is at least 2 steps away from the edge.
            let res = wave.trace_individual(bathymetry_data, &10.0, &50.0, &0.12, &0.0, &0.0, &18.0, &1.0).unwrap();
        
            let _ = output_to_tsv_file("y_out.txt",&res.0, &res.1);
        }

 }