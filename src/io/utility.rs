//! Functions for creating Netcdf3 bathymetry and current files.
//!
//! Note that for the 2d array for depth or velocity, we use `y` variable to
//! represent the row and `x` variable to represent the column.

use std::path::Path;

#[allow(dead_code)]
/// Create a NetCDF3 Bathymetry File
///
/// # Arguments
/// `path` : `&Path` a reference to the path where the file will be created
///
/// `x_num` : `usize` the number of points in the x direction
///
/// `y_num` : `usize` the number of points in the y direction
///
/// `x_step` : `f32` the step size distance between points in the x direction
///
/// `y_step` : `f32` the step size distance between points in the y direction
///
/// `depth_fn` : `fn(f32,f32) -> f64` a function that maps each (x,y) input to
/// the depth, h, at that point.
/// 
/// # Example
/// Create a bathymetry file with a constant depth of 100 m and save to `path`.
///
/// create_netcdf3_bathymetry(&path, 10, 10, 100.0, 100.0, |_, _| 100.0)
pub(crate) fn create_netcdf3_bathymetry(
    path: &Path,
    x_num: usize,
    y_num: usize,
    x_step: f32,
    y_step: f32,
    depth_fn: fn(f32, f32) -> f64,
) {
    let x_data: Vec<f32> = (0..x_num).map(|x| x as f32 * x_step).collect();
    let y_data: Vec<f32> = (0..y_num).map(|y| y as f32 * y_step).collect();

    let mut depth_data: Vec<f64> = Vec::new();
    for y in &y_data {
        for x in &x_data {
            depth_data.push(depth_fn(*x, *y));
        }
    }

    // most below copied from the docs
    use netcdf3::{DataSet, FileWriter, Version};
    let y_dim_name: &str = "y";
    let y_var_name: &str = y_dim_name;
    let y_var_len: usize = y_num;

    let x_dim_name: &str = "x";
    let x_var_name: &str = x_dim_name;
    let x_var_len: usize = x_num;

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

#[allow(dead_code)]
/// Create a NetCDF3 current snapshot (no time)
///
/// # Arguments
/// `path` : `&Path` a reference to the path where the file will be created
///
/// `x_num` : `usize` the number of points in the x direction
///
/// `y_num` : `usize` the number of points in the y direction
///
/// `x_step` : `f32` the step size distance between points in the x direction
///
/// `y_step` : `f32` the step size distance between points in the y direction
///
/// `current_fn` : `fn(f32,f32)->(f64,f64)` a function thatr maps each (x,y)
/// input to the current (u,v) at that point.
pub(crate) fn create_netcdf3_current(
    path: &Path,
    x_num: usize,
    y_num: usize,
    x_step: f32,
    y_step: f32,
    current_fn: fn(f32, f32) -> (f64, f64),
) {
    let x_data: Vec<f32> = (0..x_num).map(|x| x as f32 * x_step).collect();
    let y_data: Vec<f32> = (0..y_num).map(|y| y as f32 * y_step).collect();

    let mut u_data: Vec<f64> = Vec::new();
    let mut v_data: Vec<f64> = Vec::new();

    for y in &y_data {
        for x in &x_data {
            let (u, v) = current_fn(*x, *y);
            u_data.push(u);
            v_data.push(v);
        }
    }

    // most below copied from the docs
    use netcdf3::{DataSet, FileWriter, Version};
    let y_dim_name: &str = "y";
    let y_var_name: &str = y_dim_name;
    let y_var_len: usize = y_num;

    let x_dim_name: &str = "x";
    let x_var_name: &str = x_dim_name;
    let x_var_len: usize = x_num;

    let u_dim_name: &str = "u";
    let u_var_name: &str = u_dim_name;
    let u_var_len: usize = u_data.len();

    let v_dim_name: &str = "v";
    let v_var_name: &str = v_dim_name;
    let v_var_len: usize = v_data.len();

    // Create the NetCDF-3 definition
    // ------------------------------
    assert_eq!(u_var_len, y_var_len * x_var_len);
    assert_eq!(v_var_len, y_var_len * x_var_len);
    let data_set: DataSet = {
        let mut data_set: DataSet = DataSet::new();
        // Define the dimensions
        data_set.add_fixed_dim(y_dim_name, y_var_len).unwrap();
        data_set.add_fixed_dim(x_dim_name, x_var_len).unwrap();
        // Define the variable
        data_set.add_var_f32(y_var_name, &[y_dim_name]).unwrap();
        data_set.add_var_f32(x_var_name, &[x_var_name]).unwrap();
        data_set
            .add_var_f64(u_var_name, &[y_dim_name, x_var_name])
            .unwrap();
        data_set
            .add_var_f64(v_var_name, &[y_dim_name, x_var_name])
            .unwrap();

        data_set
    };

    // Create and write the NetCDF-3 file
    // ----------------------------------
    let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
    // Set the NetCDF-3 definition
    file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
    file_writer.write_var_f32(y_var_name, &y_data[..]).unwrap();
    file_writer.write_var_f32(x_var_name, &x_data[..]).unwrap();
    file_writer.write_var_f64(u_var_name, &u_data[..]).unwrap();
    file_writer.write_var_f64(v_var_name, &v_data[..]).unwrap();
    file_writer.close().unwrap();
    // end of copied from docs
}
