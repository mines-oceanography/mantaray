//! Structures and functions to assist with reading and writing input and output
//!
//! Data types supported:
//! - netcdf4: reading bathymetry file
//! - netcdf3: creating files

mod netcdf;
pub(crate) mod utility;

use std::collections::HashMap;

use crate::error::Result;

pub(crate) trait Dataset {
    /// Get the length of a dimension
    fn dimension_len(&self, name: &str) -> Result<usize>;

    // Better move this to an iterator instead of a vector
    #[allow(dead_code)]
    /// Get the names of the available variables
    fn varnames(&self) -> Vec<String>;

    /// Get the values of a variable
    fn values(&self, name: &str) -> Result<ndarray::ArrayD<f64>>;

    /// Get the value of a variable at a specific index
    fn get_variable(&self, name: &str, i: usize, j: usize) -> Result<f32>;

    /// Get the order of the dimensions for a variable
    fn dimensions_order(&self, varname_x: &str, varname_y: &str) -> HashMap<String, String>;
}
