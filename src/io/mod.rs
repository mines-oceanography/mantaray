//! Structures and functions to assist with reading and writing input and output
//!
//! Data types supported:
//! - netcdf4: reading bathymetry file
//! - netcdf3: creating files

mod netcdf;
pub mod utility;

use crate::error::{Error, Result};

trait Dataset {
    fn dimension_len(&self, name: &str) -> Result<usize>;
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
