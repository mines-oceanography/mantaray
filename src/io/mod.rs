//! Structures and functions to assist with reading and writing input and output
//!
//! Data types supported:
//! - netcdf4: reading bathymetry file
//! - netcdf3: creating files

mod netcdf;
pub mod utility;

use std::collections::HashMap;

use tracing::trace;

use crate::error::{Error, Result};

trait Dataset {
    fn dimension_len(&self, name: &str) -> Result<usize>;
    // Better move this to an iterator instead of a vector
    fn varnames(&self) -> Vec<String>;
    fn values(&self, name: &str) -> Result<ndarray::ArrayD<f64>>;
    fn get_variable(&self, name: &str, i: usize, j: usize) -> Result<f32>;
}

#[allow(dead_code)]
#[derive(Debug)]
/// Holds and apply a linear relationship
///
/// This was originally created to handle the dimensions of regular girds,
/// such latitude and longitude in a regularly spaced dataset, providing a
/// cheap conversion between the dimension, such as latitute, to the
/// correspondent index position.
struct LinearFit<T> {
    slope: T,
    intercept: T,
}

impl<T> LinearFit<T>
where
    T: Copy,
    T: std::ops::Sub<Output = T>,
    T: std::ops::Mul<Output = T>,
{
    #[allow(dead_code)]
    /// Predict the index position of a given value
    fn predict(&self, x: T) -> T {
        (x - self.intercept) * self.slope
    }
}

impl LinearFit<f64> {
    /// Create a new LinearFit from a vector of values
    ///
    /// For a given vector of values, calculates the best linear relationship
    /// between the values and its index position with the purpose to estimate
    /// the index position of the closest value to a given target.
    ///
    /// This procedure also validates if a linear relationship is a good
    /// approximation with a threshold of 0.5% of tolerance.
    fn from_fit(x: ndarray::ArrayD<f64>) -> Result<LinearFit<f64>> {
        let dx = &x.slice(ndarray::s![1..]) - &x.slice(ndarray::s![..-1]);
        let slope = dx.mean().unwrap();
        let criteria = ((dx - slope) / slope)
            .abs()
            .into_iter()
            .map(|v| v > 0.005)
            .any(|v| v);
        if criteria {
            return Err(Error::Undefined);
        }
        Ok(LinearFit {
            slope,
            intercept: x[0],
        })
    }

    /*
    fn from_fit(x: Vec<f32>) -> Result<Self> {
        let delta: Vec<_> = x.windows(2).map(|v| v[1] - v[0]).collect();
        let slope = (delta.len() as f32 - 1.0) / delta.iter().sum::<f32>();
        let threshold = delta
            .iter()
            .map(|v| (v - slope).abs() / slope)
            .any(|v| v > 0.005);
        if threshold {
            return Err(Error::Undefined);
        }

        let intercept = x[0];
        Ok(LinearFit { slope, intercept })
    }
    */
}

struct RegularGrid {
    dataset: netcdf::File,
    x_size: usize,
    x_map: LinearFit<f64>,
    y_size: usize,
    y_map: LinearFit<f64>,
    // Save the dimensions order: ij or ji
    dimension_order: HashMap<String, String>,
}

impl RegularGrid {
    /*
        fn validate_as_regular_grid() {
            // Open a full 1D array for x and another for y
            // calculate delta and confirm that all are <0.01 diference (criteria for linear)
            //
        }
    */
    #[allow(dead_code)]
    fn open<P>(file: P, varname_x: &str, varname_y: &str) -> Self
    where
        P: AsRef<std::path::Path>,
    {
        // confirm it is linear
        // Define A & B coefficients
        // get i_size and j_size

        let dataset = netcdf::open(&file).unwrap();

        // Identify the variables that have the user defined dimensions
        // and create a map on the dimenson order
        let varnames = &dataset
            .variables()
            .into_iter()
            .filter(|v| {
                v.dimensions()
                    .into_iter()
                    .map(|v| v.name() == varname_x)
                    .any(|v| v)
            })
            .filter(|v| {
                v.dimensions()
                    .into_iter()
                    .map(|v| v.name() == varname_y)
                    .any(|v| v)
            })
            .filter_map(|v| {
                match &v
                    .dimensions()
                    .into_iter()
                    .map(|v| v.name())
                    .collect::<Vec<_>>()[..]
                {
                    [varname_x, varname_y] => Some((v.name(), "xy".to_string())),
                    [varname_y, varname_x] => Some((v.name(), "yx".to_string())),
                    _ => None,
                }
            })
            .collect::<Vec<_>>();

        let dimension_order: HashMap<String, String> = HashMap::from_iter(varnames.iter().cloned());

        let x_size = dataset.dimension_len(varname_x).unwrap();

        let x_map = LinearFit::from_fit(dataset.values(varname_x).unwrap()).unwrap();

        let y_size = dataset.dimension_len(varname_y).unwrap();
        let y_map = LinearFit::from_fit(dataset.values(varname_x).unwrap()).unwrap();

        Self {
            dataset,
            x_size,
            x_map,
            y_size,
            y_map,
            dimension_order,
        }
    }

    #[allow(dead_code)]
    /// Get the nearest `varname` value to the given `x` and `y` coordinates
    fn nearest(&self, varname: &str, point: Point<f64>) -> Result<f32> {
        let i = self.x_map.predict(*point.x()).round() as usize;
        if i >= self.x_size {
            return Err(Error::IndexOutOfBounds);
        }
        let j = self.y_map.predict(*point.y()).round() as usize;
        if j >= self.y_size {
            return Err(Error::IndexOutOfBounds);
        }

        match self.dimension_order.get(varname) {
            Some(v) => match v.as_str() {
                "xy" => {
                    trace!("Assuming dimension order is 'xy'");
                    self.dataset.get_variable(varname, i, j)
                }
                "yx" => {
                    trace!("Assuming dimension order is 'yx'");
                    self.dataset.get_variable(varname, j, i)
                }
                _ => return Err(Error::Undefined),
            },
            _ => return Err(Error::Undefined),
        }
    }
}
