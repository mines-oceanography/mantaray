//! Structures and functions to assist with reading and writing input and output
//!
//! Data types supported:
//! - netcdf4: reading bathymetry file
//! - netcdf3: creating files

mod netcdf;
pub mod utility;

use crate::error::Result;

trait Dataset {
    fn dimension_len(&self, name: &str) -> Result<usize>;
    fn values(&self, name: &str) -> Result<ndarray::ArrayD<f64>>;
    fn get_variable(&self, name: &str, i: usize, j: usize) -> Result<f32>;
}
