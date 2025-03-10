//! Structures and functions to assist with reading and writing input and output
//!
//! Data types supported:
//! - netcdf4: reading bathymetry file
//! - netcdf3: creating files

mod netcdf;
pub mod utility;

use std::collections::HashMap;

use crate::error::Result;

pub(crate) trait Dataset {
    fn dimension_len(&self, name: &str) -> Result<usize>;
    // Better move this to an iterator instead of a vector
    fn varnames(&self) -> Vec<String>;
    fn values(&self, name: &str) -> Result<ndarray::ArrayD<f64>>;
    fn get_variable(&self, name: &str, i: usize, j: usize) -> Result<f32>;
    fn dimensions_order(&self, varname_x: &str, varname_y: &str) -> HashMap<String, String>;
}
