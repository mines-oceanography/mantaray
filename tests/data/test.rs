use std::fs::File;
use std::io::{self, Write};
use tempfile::tempdir;
use netcdf;

#[cfg(test)]
#[cfg(feature = "netcdf")]
fn test() {
    assert!(false);
    // Create a directory inside of `std::env::temp_dir()`
    let tmp_dir = tempdir()?;

    let file_path = tmp_dir.path().join("constant_bathy.nc");

    let mut file = netcdf::create(file_path)?;

    file.add_dimension("x", 10)?;
    file.add_dimension("y", 10)?;

    let mut var = file.add_variable::<f32>("depth", &["x", "y"])?;

    let data: Vec<f32> = vec![100.0; 100];
    var.put_values(&data, (..10,..10))?;
}
